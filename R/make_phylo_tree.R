#' Build a `phyloTree` object for geneplast
#'
#'
#'
#' @param sspids a vector or data frame containing NCBI Taxon IDs from the species of interest.
#' @param newick a rooted phylogenetic tree in Newick format.
#' @param verbose a logical value specifying whether or not to display detailed messages.
#' @return An object of class "phylo".
#' @export
#' @author Danilo O Imparato
make.phyloTree <- function(sspids=NULL, newick=NULL, verbose=TRUE){
    check_args(sspids, newick)
    if(!is.null(newick)) {
        phylo <- treeio::read.newick(newick)
        return(phylo)
    }
    if(verbose)cat("-Building eukaryotes tree from sspids...\n")
    taxa_of_interest <- check_sspids(sspids)
    download_and_extract(verbose)
    bfc <- .get_cache()
    ncbi_merged_ids <- readr::read_delim(
      paste0(BiocFileCache::bfccache(bfc), "/merged.dmp"),
      delim = "|",
      trim_ws = TRUE,
      col_names = c("taxid","new_taxid"),
      col_types = readr::cols_only(
        taxid = "c",
        new_taxid = "c"
      )
    )
    ncbi_updated_ids_lookup <- tibble::deframe(ncbi_merged_ids)
    ncbi_edgelist <- readr::read_delim(
      paste0(BiocFileCache::bfccache(bfc), "/nodes.dmp"),
      skip = 1,
      delim = "|",
      trim_ws = TRUE,
      col_names = c("n1","n2","rank"),
      col_types = "ccc"
    )
    ncbi_taxon_names <- readr::read_delim(
      paste0(BiocFileCache::bfccache(bfc), "/names.dmp"),
      delim = "|",
      trim_ws = TRUE,
      col_names = c("name","ncbi_name","type"),
      col_types = "cc-c"
    )
    taxa_of_interest_updated <-
      taxa_of_interest |>
      dplyr::left_join(ncbi_merged_ids) |>
      dplyr::mutate(new_taxid = dplyr::coalesce(new_taxid, taxid))

    if(verbose)cat("-Loading full TimeTree species tree...\n")
    timetree_tree <- ape::read.tree(
      file = system.file("extdata", "TTOL v5.1 Species Tree - Internal NCBI ID Labelled.nwk",
                         package = "geneplast.data", mustWork = TRUE)
    )
    # We cannot build a working tree with multiple empty labels
    # Therefore we assign unlabelled nodes an arbitrary label. eg: "fix_123"
    unlabelled_nodes_indexes <- timetree_tree$node.label == ""
    number_of_unlabelled_nodes <- sum(unlabelled_nodes_indexes)
    timetree_tree$node.label[unlabelled_nodes_indexes] <- paste0("fix_", 1:number_of_unlabelled_nodes)

    # It is much easier to work with the tree as an edgelist dataframe
    # One way to obtain the df is by converting it into an igraph obj, then df
    timetree_graph <- ape::as.igraph.phylo(timetree_tree, directed = TRUE)
    timetree_edgelist <- igraph::as_data_frame(timetree_graph, what = "edges")
    timetree_edgelist$from <- gsub("'", '', timetree_edgelist$from)
    timetree_edgelist$to <- gsub("'", '', timetree_edgelist$to)
    temp_timetree_edgelist_from <- ncbi_updated_ids_lookup[timetree_edgelist$from]
    temp_timetree_edgelist_to <- ncbi_updated_ids_lookup[timetree_edgelist$to]
    timetree_edgelist_updated <- timetree_edgelist
    timetree_edgelist_updated$from <- dplyr::coalesce(temp_timetree_edgelist_from, timetree_edgelist$from)
    timetree_edgelist_updated$to <- dplyr::coalesce(temp_timetree_edgelist_from, timetree_edgelist$to)

    if(verbose)cat("-Searching for missing taxa...\n")
    # Find missing taxa in TimeTree
    # These will be further grafted from NCBI Taxonomy tree
    # The search includes both tips and internal nodes from the tree.
    taxa_of_interest_not_in_timetree <- ! taxa_of_interest_updated$new_taxid %in% c(timetree_edgelist_updated$from, timetree_edgelist_updated$to)
    taxa_of_interest_not_in_timetree <- taxa_of_interest_updated[taxa_of_interest_not_in_timetree,]

    if(verbose)cat("-Mapping taxon IDs and scientific names...\n")
    ncbi_scientific_names <-
      ncbi_taxon_names |>
      dplyr::filter(type == "scientific name") |>
      dplyr::select(name, ncbi_name)
    ncbi_scientific_names_lookup <- tibble::deframe(ncbi_scientific_names)
    eukaryota_taxon_id <- subset(ncbi_taxon_names, ncbi_name == "Eukaryota", "name", drop = TRUE)
    ncbi_graph <- igraph::graph_from_data_frame(
      d = ncbi_edgelist[,2:1],
      directed = TRUE
    )
    rm(ncbi_edgelist, ncbi_merged_ids)
    gc()
    ncbi_tree <- treeio::as.phylo(ncbi_graph)
    if(verbose)cat("-NCBI tree created.\n")

    # Vertices from taxa of interest not present in TimeTree must be grafted
    taxids_not_in_timetree <- taxa_of_interest_not_in_timetree$new_taxid

    # For each missing taxon, go through each ancestral vertex checking if it exists in TimeTree
    ncbi_graph_vids <- igraph::V(ncbi_graph) |> as.numeric() |> as.character()
    ncbi_graph_vnames <- igraph::V(ncbi_graph)$name
    ncbi_graph_vid_to_vname <- setNames(ncbi_graph_vnames, ncbi_graph_vids)

    if(verbose && length(taxids_not_in_timetree) > 0) {
        "%i missing taxa. Looking for the closest available taxa up in the hierarchy. This might take a while..." |>
            sprintf(length(taxids_not_in_timetree)) |>
            cat()
    }
    timetree_edgelist_patch <- purrr::map_df(
      .x = taxids_not_in_timetree,
      .f = function(current_taxid) {
        root_vertex_id <- as.numeric(igraph::V(ncbi_graph)[[current_taxid]])

        found_ancestor <- NULL

        igraph::dfs(
          graph = ncbi_graph,
          root = root_vertex_id,
          mode = "in",
          unreachable = FALSE,
          order = FALSE,
          in.callback = function(graph, current_vertex, extra) {

            vertex_id <- as.character(current_vertex[["vid"]])
            vertex_taxid <- ncbi_graph_vid_to_vname[vertex_id]

            ancestor_is_in_timetree <- vertex_taxid %in% timetree_edgelist_updated$from

            if(ancestor_is_in_timetree) {
              found_ancestor <<- vertex_taxid
              return(TRUE)
            }

            return(FALSE)
          }
        )

        data.frame(from = found_ancestor, to = current_taxid)
      },
      .progress = TRUE
    )

    if(verbose)cat("-Grafting missing vertices...\n")
    # Grafting missing vertices to the TimeTree edgelist
    timetree_edgelist_patched <- rbind(timetree_edgelist_updated, timetree_edgelist_patch)

    # It can happen that, after forceful grafting, some edges are duplicated or link to themselves
    # This is not allowed, therefore we remove such occurrences
    timetree_edgelist_patched <- unique(timetree_edgelist_patched) |>
      dplyr::filter(from != to)

    # Creating an igraph from the grafted tree
    timetree_graph_patched <- igraph::graph_from_data_frame(timetree_edgelist_patched, directed = TRUE)
    if(verbose)cat("-Grafted tree created.\n")

    # Extract the eukaryotic subtree from the complete tree for just the taxa of interest.
    # Any non-eukaryotes are discarded after it.
    eukaryote_vertex_in_patched_timetree <- igraph::V(timetree_graph_patched)[eukaryota_taxon_id]

    leaves_of_interest_in_patched_timetree <- igraph::V(timetree_graph_patched)[taxa_of_interest_updated$new_taxid]

    # Traversing the TimeTree graph from the Eukaryota vertex towards the taxa of interest
    eukaryote_paths_in_patched_timetree <- igraph::shortest_paths(
      graph = timetree_graph_patched,
      from = eukaryote_vertex_in_patched_timetree,
      to = leaves_of_interest_in_patched_timetree,
      mode = "out"
    )$vpath

    eukaryote_vertices_in_patched_timetree <- eukaryote_paths_in_patched_timetree |> unlist() |> unique()

    timetree_eukaryote_subgraph <- igraph::induced_subgraph(
      graph = timetree_graph_patched,
      vids = eukaryote_vertices_in_patched_timetree,
      impl = "create_from_scratch"
    )

    timetree_tree_patched <- treeio::as.phylo(timetree_eukaryote_subgraph)
    timetree_tree_patched$tip.alias <- ncbi_scientific_names_lookup[timetree_tree_patched$tip.label]
    timetree_tree_patched$edge.length <- NULL
    if(verbose)cat("-Eukaryotes tree created.\n")
    timetree_tree_patched
}
