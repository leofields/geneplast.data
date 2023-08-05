.get_cache <- function() {
    cache <- tools::R_user_dir("geneplastdata", which="cache")
    BiocFileCache::BiocFileCache(cache)
}

download_and_extract <- function( verbose = FALSE ) {
    fileURL <- "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

    bfc <- .get_cache()
    rid <- bfcquery(bfc, "taxdump", "rname")$rid
    if (!length(rid)) {
        if( verbose )
            message( "Downloading and extracting NCBI taxdump files" )
        rid <- names(bfcadd(bfc, "taxdump", fileURL ))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid))) {
        options(timeout=300)
        bfcdownload(bfc, rid)
    }

    rpath <- bfcrpath(bfc, rids = rid)
    untar(rpath, exdir = bfccache(bfc))
}
