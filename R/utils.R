
download_if_missing <- function(url, filename = basename(url)) {
    filename <- here::here("data-raw/download", filename)
    if (!file.exists(filename)) {
        download.file(url, filename)
    }
}

check_args <- function(arg1, arg2) {
    b1<-is.character(arg1) || is.integer(arg1) || is.numeric(arg1)
    b2<-is.matrix(arg1) || is.data.frame(arg1)
    b3<-is.null(arg2)
    if(!b1 && !b2 && b3){
        stop("A valid 'sspids' argument or 'newick' file must be provided! \n",
             call.=FALSE)
    }
}

check_sspids <- function(arg) {
    b1<-is.character(arg) || is.integer(arg) || is.numeric(arg)
    b2<-is.matrix(arg) || is.data.frame(arg)
    if(!b1 && !b2){
        stop("'sspids' should be a vector of characters or dataframe! \n",
             call.=FALSE)
    }
    if(b1){
        arg<-unique(as.character(arg))
        arg<-arg[!is.na(arg)]
        arg<-arg[arg!='']
        arg<-sort(arg)
        arg<-data.frame(taxid=arg,stringsAsFactors=FALSE)
    } else {
        clpars<-c("taxid")
        clname<-tolower(colnames(arg))
        if(!all(clpars%in%clname)){
            stop("'sspids' colnames should include: ",paste(clpars,collapse=", "),
                 call.=FALSE)
        }
        colnames(arg)<-clname
        arg<-arg[,c(clpars,clname[which(!clname%in%clpars)]),drop=FALSE]
        arg[,1]<-as.character(arg[,1])
        arg[,2]<-as.character(arg[,2])
        for(i in 1:ncol(arg)){
            if(!is.numeric(arg[,i]) || !is.integer(arg[,i])){
                arg[,i]<-as.character(arg[,i])
            }
        }
        arg<-data.frame(arg,stringsAsFactors=FALSE)
        uni<-unique(arg[,1])
        uni<-uni[!is.na(uni)]
        uni<-uni[uni!='']
        uni<-sort(uni)
        arg<-arg[match(uni,arg[,1]), ,drop=FALSE]
    }
    rownames(arg)<-arg[,1]
    return(arg)
}
