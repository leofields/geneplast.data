Title <- c("STRING database v11.0", "OMA Browser All.Jan2020", "OrthoDB v10.1")
Description <-
    c("Eukaryotes input data for geneplast package extracted from STRING v11.0",
      "Eukaryotes input data for geneplast package extracted from OMA All.Jan2020",
      "Eukaryotes input data for geneplast package extracted from OrthoDB v10.1")
BiocVersion <- "3.12"
Genome <- NA
SourceType <- "TXT"
SourceUrl <- c("https://string-db.org/", "https://omabrowser.org/", "https://www.orthodb.org/")
SourceVersion <- c("11.0", "All.Jan2020", "v10.1")
Species <- NA
TaxonomyId <- NA
Coordinate_1_based <- FALSE
DataProvider <- c("STRING", "OMA", "OrthoDB")
Maintainer <- "Leonardo RS Campos <leofields@gmail.com>"
RDataClass <- "Rda"
DispatchClass <- "Rda"
ResourceName <- c("gpdata_string_v11.RData", "gpdata_oma_All.Jan2020.RData", "gpdata_orthodb_v101.RData")
RDataPath <- c("geneplast.data.string.v11", "geneplast.data.oma.All.Jan2020", "geneplast.data.orthodb.v101")
Tags <- "geneplast, datasets, eukaryota, orthology"

meta <- data.frame(Title,
                   Description,
                   BiocVersion,
                   Genome,
                   SourceType,
                   SourceUrl,
                   SourceVersion,
                   Species,
                   TaxonomyId,
                   Coordinate_1_based,
                   DataProvider,
                   Maintainer,
                   RDataClass,
                   DispatchClass,
                   ResourceName,
                   RDataPath,
                   Tags)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)