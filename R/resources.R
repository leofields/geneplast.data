
.get_cache <- function() {
    cache <- tools::R_user_dir("geneplastdata", which="cache")
    BiocFileCache::BiocFileCache(cache)
}

download_and_extract <- function( verbose = FALSE ) {
    fileURL <- "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, "taxdump", "rname")$rid
    if (!length(rid)) {
        if( verbose )
            message( "Downloading and extracting NCBI taxdump files" )
        rid <- names(BiocFileCache::bfcadd(bfc, "taxdump", fileURL ))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
        options(timeout=300)
        BiocFileCache::bfcdownload(bfc, rid)
    }

    rpath <- BiocFileCache::bfcrpath(bfc, rids = rid)
    untar(rpath, exdir = BiocFileCache::bfccache(bfc))
}
