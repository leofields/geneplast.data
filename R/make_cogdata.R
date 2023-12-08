#' Parse orthogroups tabular output from OrthoFinder into a `cogdata` data frame for geneplast
#'
#'
#' @param file OrthoFinder orthogroups tabular file
#' @return cogdata data frame
#' @export
#' @author Leonardo RS Campos
make.cogdata <- function(file){
    orthogroups <- readr::read_tsv(file)
    orthogroups |>
        tidyr::pivot_longer(cols = !1, names_to = "ssp_id", values_to = "protein_id") |>
        tidyr::separate_rows(protein_id) |>
        dplyr::rename(cog_id = 1)
}
