#!/usr/bin/env Rscript
###############################################################################
#                                                                             #
# Script for creating geneplast objects from OMA Browser database             #
#                                                                             #
# Authors: Leonardo RS Campos (leofields@gmail.com)                           #
#          Danilo Imparato (imparatodanilo@gmail.com)                         #
#                                                                             #
###############################################################################

# Packages ####
library(ape)
library(phytools)
library(tidyverse)
library(igraph)
library(tidylog)
library(magrittr)
library(treeio)

# Downloading files ####
if (!file.exists("working_directory/taxonomy")) {
  tmp <- tempfile()
  download.file("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", tmp)
  untar(tmp, exdir = "working_directory/taxonomy")
  unlink(tmp)
}

if (!file.exists("working_directory/oma-species.txt")) {
  download.file(
    "https://omabrowser.org/All/oma-species.txt",
    "working_directory/oma-species.txt"
  )
}

if (!file.exists("working_directory/oma-groups.txt.gz")) {
  download.file(
    "https://omabrowser.org/All/oma-groups.txt.gz",
    "working_directory/oma-groups.txt.gz"
  )
}

if (!file.exists("working_directory/timetree/Eukaryota_species.nwk")) {
  download.file(
    "http://www.timetree.org/ajax/direct_download?direct-download-format=newick&direct-download-id=23070&direct-download-rank=species",
    "working_directory/timetree/Eukaryota_species.nwk"
  )
}

# Loading data ####
oma_species <- read_tsv(
  "working_directory/oma-species.txt",
  skip = 3,
  col_names = c("oma_code", "taxid", "oma_name"),
  col_types = cols_only("c", "c", "c", "-", "-")
)

ncbi_merged_ids <- read_delim(
  "working_directory/taxonomy/merged.dmp",
  delim = "|",
  trim_ws = TRUE,
  col_names = c("taxid","new_taxid"),
  col_types = cols_only("c", "c", "-")
)

ncbi_edgelist <- read_delim(
  "working_directory/taxonomy/nodes.dmp",
  skip = 1,
  delim = "|",
  trim_ws = TRUE,
  col_names = c("n1", "n2"),
  col_types = cols_only("c", "c", "-", "-", "-", "-", "-", "-", "-", "-", "-",
                        "-", "-", "-")
)

ncbi_taxon_names <- read_delim(
  "working_directory/taxonomy/names.dmp",
  delim = "|",
  trim_ws = TRUE,
  col_names = c("name","ncbi_name","type"),
  col_types = cols_only("c", "c", "-", "c", "-")
)

# Creating graph ####
# leaving only scientific names
ncbi_taxon_names %<>%
  filter(type == "scientific name") %>%
  select(name, ncbi_name)

# finding Eukaryota taxid
eukaryota_taxon_id <-
  subset(ncbi_taxon_names, ncbi_name == "Eukaryota", "name", drop = TRUE)

# updating taxids
oma_species %<>%
  left_join(ncbi_merged_ids) %>%
  mutate(new_taxid = coalesce(new_taxid, taxid))

# creating graph
g <- graph_from_data_frame(ncbi_edgelist[, 2:1], directed = TRUE,
                           vertices = ncbi_taxon_names)

# easing memory
rm(ncbi_edgelist, ncbi_merged_ids)

# Traversing graph from the "Eukaryota" node to its leaves ####
eukaryote_root <- V(g)[eukaryota_taxon_id]
eukaryote_leaves <- V(g)[oma_species[["new_taxid"]]]

eukaryote_paths <- shortest_paths(g, from = eukaryote_root, to = eukaryote_leaves, mode = "out")$vpath

eukaryote_vertices <- eukaryote_paths %>% unlist %>% unique

eukaryote_tree <- induced_subgraph(g, eukaryote_vertices, impl = "create_from_scratch")

# Saving ####
ncbi_tree <- treeio::as.phylo(eukaryote_tree)

oma_eukaryotes <- oma_species %>%
  filter(new_taxid %in% ncbi_tree$tip.label) %>%
  inner_join(ncbi_taxon_names, by = c("new_taxid" = "name"))

write(oma_eukaryotes[["ncbi_name"]],"working_directory/oma_eukaryotes.txt")
###############################################################################
#                 EXTERNAL PROCEDURE - MANUAL INPUT #1                     ####
###############################################################################
#                                                                             #
# Access http://www.timetree.org/ website, under                              #
# BUILD A TIMETREE section, Load a List of Species,                           #
# choose the file oma_eukaryotes.txt and click the 'Upload' button.           #
#                                                                             #
# After processing, on the left pane, under Export Tree section,              #
# click the 'To Newick File' button, and save in the working directory.       #
#                                                                             #
###############################################################################
readline(prompt="EXTERNAL PROCEDURE - MANUAL INPUT #1\nPress ENTER when ready.")
if (!file.exists("working_directory/oma_eukaryotes.nwk")) {
  stop("Missing manual input file obtained from external procedure #1")
}

# loading newick tree obtained from timetree after external procedure
timetree_newick <- read.tree("working_directory/oma_eukaryotes.nwk")

# replacing timetree's underscores with spaces
timetree_newick[["tip.label"]] %<>% str_replace_all("_", " ")

# which timetree species' names exactly match with ncbi's
taxid_indexes <- timetree_newick[["tip.label"]] %>% match(oma_eukaryotes[["ncbi_name"]])

# find out which timetree species names didn't exactly match ncbi's
unmatched_names <- timetree_newick[["tip.label"]] %>% magrittr::extract(taxid_indexes %>% is.na)

ncbi_to_timetree <-
  data.frame(
    timetree_name = unmatched_names,
    ncbi_name = c(""),
    stringsAsFactors = FALSE
  )
write.csv(ncbi_to_timetree,file="working_directory/oma_ncbi_to_timetree.csv", row.names = FALSE)
###############################################################################
#                 EXTERNAL PROCEDURE - MANUAL INPUT #2                     ####
###############################################################################
#                                                                             #
# Manually fill the ncbi_name column in ncbi_to_timetree.csv file.            #
#                                                                             #
# To help this procedure of matching names, it is useful to search NCBI       #
# Taxonomy Browser checking for synonyms or merged taxons.                    #
#                                                                             #
###############################################################################
readline(prompt="EXTERNAL PROCEDURE - MANUAL INPUT #2\nPress ENTER when ready.")
# loading manually created lookup table to be joined
ncbi_to_timetree <- read_csv("working_directory/oma_ncbi_to_timetree.csv")

# joining info
species_dictionary <- oma_eukaryotes %>% left_join(ncbi_to_timetree)

# coalescing NAs to ncbi_name
species_dictionary %<>% mutate(timetree_name = coalesce(timetree_name, ncbi_name))
species_dictionary %<>% mutate(timetree_name = ifelse(timetree_name %in% timetree_newick[["tip.label"]], timetree_name, NA))

# annotating genera
species_dictionary %<>% mutate(genus_search = coalesce(timetree_name, ncbi_name) %>% strsplit(" ") %>% sapply("[", 1))

# unique genera
selected_genera <- species_dictionary[["genus_search"]] %>% unique

# the following genera names are unreliable and should not be searched for
duplicated_genera <- scan("inst/extdata/duplicated_genera.txt", what = "character")

# these are unreliable selected_genera:
unreliable_genera <- intersect(selected_genera, duplicated_genera)

# loading newick tree obtained from timetree
tree_85k <- read.tree("working_directory/timetree/Eukaryota_species.nwk")

# ensuring a cleaner newick file with only necessary data
# this is actually really important
tree_85k[["node.label"]] <- NULL
tree_85k[["edge.length"]] <- NULL

# replacing timetree's underscores with spaces
tree_85k[["tip.label"]] %<>% str_replace_all("_", " ")

# storing genus
tree_85k[["tip.genus"]] <- sapply(strsplit(tree_85k[["tip.label"]]," "), "[", 1)
tree_85k_genera <- tree_85k[["tip.genus"]] %>% unique

# subtracting unreliable genera
tree_85k_genera %<>% setdiff(unreliable_genera)

# keeping only selected genera, including unreliable ones
tree_genus <- tree_85k %$% keep.tip(., tip.label[tip.genus %in% selected_genera])
tree_genus[["tip.genus"]] <- sapply(strsplit(tree_genus[["tip.label"]]," "), "[", 1)

# unfound species which genera are present in the 85k tree
unfound_species <- species_dictionary %>% filter(is.na(timetree_name) & genus_search %in% tree_85k_genera)

# for each unfound species which genus is present in the 85k tree,
for(i in 1:nrow(unfound_species)){
  # we search for all species of this genus ("sister species") in the 85k tree
  # this part is tricky because bind.tip rebuilds the tree from scratch
  # so we need to keep removing underscores. there are better ways to do this.
  tip_genus <- tree_genus[["tip.label"]] %>% strsplit("[_ ]") %>% sapply("[", 1)
  sister_species <- tree_genus[["tip.label"]][tip_genus == unfound_species[[i, "genus_search"]]]
  # we obtain the sister_species' most recent common ancestor (MRCA)
  # c(.[1]) is a hack because the MRCA function only works with at least 2 nodes
  where <- getMRCA(tree_genus, sister_species %>% c(.[1]))
  # and then add a leaf node linked to this MRCA
  tree_genus %<>% bind.tip(tip.label = unfound_species[[i, "ncbi_name"]], where = where)
}

# for some reason bind.tip adds underscores to species names
tree_genus[["tip.label"]] %<>% str_replace_all("_", " ")

# keeping track of found species
found_species <- species_dictionary %>% filter(!is.na(timetree_name) | genus_search %in% tree_85k_genera)
# forced_name means it either was found in timetree or we forced it by looking at genera names
found_species %<>% mutate(forced_name = coalesce(timetree_name, ncbi_name))

# so we keep only found species in this tree we are building (timetree + forced by genera)
tree_genus %<>% keep.tip(found_species[["forced_name"]])

# converting to ncbi taxids
tree_genus[["tip.label"]] <- found_species[["new_taxid"]][match(tree_genus[["tip.label"]], found_species[["forced_name"]])]


##################################################
##                UNFOUND GENERA                ##
##################################################
## the following part fills in missing species  ##
## by searching for its closest relatives in    ##
## in the NCBI tree                             ##
##################################################

# converting ncbi phylo to igraph
graph_ncbi <- as.igraph.phylo(ncbi_tree, directed = TRUE)

# converting phylo to igraph
graph_genus <- as.igraph.phylo(tree_genus, directed = TRUE)

# for each species which genus is not in timetree
# we'll look for its two closest species (in the NCBI tree) which are present in the tree_genus we just built
unfound_genera <- species_dictionary %>% filter(is.na(timetree_name) & !genus_search %in% tree_85k_genera)

# this is the igraph equivalent of "phylo_tree$tip.label"
tip_nodes <- V(graph_ncbi)[degree(graph_ncbi, mode = "out") == 0]

# undirected distances between all species nodes
tip_distances <- graph_ncbi %>%
  distances(v = tip_nodes, to = tip_nodes, mode = "all") %>%
  as_tibble(rownames = "from") %>%
  pivot_longer(-from, names_to = "to", values_to = "distance")

# removing self references (zero distances)
tip_distances %<>% filter(distance > 0)

# we only want to search for species of unfound genera
tip_distances %<>% inner_join(unfound_genera %>% select(from = new_taxid))

# we only want to find species already present in the genus_tree
tip_distances %<>% inner_join(found_species %>% select(to = new_taxid))

# we only want the two closest relatives
tip_distances %<>%
  group_by(from) %>%
  top_n(-2, distance) %>% # top 2 smallest distances
  top_n(2, to) # more than 2 species have the same smallest distance, so we get the first ones

# out distance matrix between all nodes in tree, needed to find MRCAs
out_distances <- graph_genus %>% distances(mode = "out")

# for each species of unfound genera,
# we find the MRCA for its two closest relatives
unfound_genera_mrca <- tip_distances %>% group_by(from) %>% summarise(mrca = {
  # which rows have no infinite distances? the last one represents the MRCA
  mrca_row_index <- max(which(rowSums(is.infinite(out_distances[, to])) == 0))
  rownames(out_distances)[mrca_row_index]
})

# adding unfound genera species nodes
graph_genus %<>% add_vertices(
  nrow(unfound_genera_mrca),
  color = "red",
  attr = list(name = unfound_genera_mrca[["from"]])
)

# defining unfound genera species edges
# edges_to_add[1] -> edges_to_add[2], edges_to_add[2] -> edges_to_add[3]...
edges_to_add <-
  V(graph_genus)[unfound_genera_mrca
                 %>% select(mrca, from)
                 %>% t
                 %>% as.vector]$name

# connecting species leafs to the supposed MRCA
graph_genus %<>% add_edges(V(graph_genus)[edges_to_add])

# finally converting to phylo format
phylo_graph_genus <- treeio::as.phylo(graph_genus)

# adding tip.alias (this is not exported with write.tree)
phylo_graph_genus[["tip.alias"]] <-
  species_dictionary[["oma_name"]][match(phylo_graph_genus[["tip.label"]],
                                            species_dictionary[["new_taxid"]])]

# converting back to oma ids
phylo_graph_genus[["tip.label"]] <-
  species_dictionary[["taxid"]][match(phylo_graph_genus[["tip.label"]],
                                      species_dictionary[["new_taxid"]])]

# ensuring a cleaner newick file with only necessary data
phylo_graph_genus[["node.label"]] <- NULL
phylo_graph_genus[["edge.length"]] <- NULL

phyloTree <- phylo_graph_genus


####### COGDATA ####

###############################################################################
#                 EXTERNAL PYTHON SCRIPT                                   ####
###############################################################################
#                                                                             #
# Execute rearrange-oma-groups.py python script                               #
#                                                                             #
###############################################################################

# loading orthologous groups info
cogdata <- read_tsv(
  "working_directory/oma-groups-rearranged.txt",
  col_names = c("cog_id", "protein_id"),
  col_types = cols_only("c", "-", "c")
)

cogdata %<>%
  mutate(oma_code = protein_id %>% substr(1, 5)) %>%
  inner_join(oma_eukaryotes) %>%
  select(cog_id, protein_id, taxid)

colnames(cogdata) <- c("cog_id", "protein_id", "ssp_id")

# Fix problem in ogr.preprocess
cogdata$cog_id <- paste0("OMA", cogdata$cog_id)

# Human ENTREZ mapping
human_entrez_2_odb <- read_delim("odb10v1_genes-human-entrez.tsv",
                                 delim = "\t", escape_double = FALSE,
                                 col_names = FALSE, trim_ws = TRUE)
names(human_entrez_2_odb) <- c("protein_id", "gene_id")
cogdata <- cogdata %>% left_join(human_entrez_2_odb)

####### SSPIDS ####
sspids <- oma_eukaryotes[,c("taxid", "oma_name")]
colnames(sspids) <- c("ssp_id", "ssp_name")
sspids$domain <- c("Eukaryota")

cogids <- data.frame(cog_id=unique(cogdata$cog_id))
row.names(cogids) <- cogids$cog_id

save(cogdata,
     sspids,
     cogids,
     phyloTree,
     file = "inst/extdata/gpdata_oma_All.Jan2020.RData",
     compress = "xz")
