
#######
#### functions hidden to the user
####################################

.findPerl <- function(perl, verbose = "FALSE")
{
  errorMsg <- "perl executable not found. Use perl= argument to specify the correct path." 
  if (missing(perl))
    {
      perl = "perl"
    }
  perl = Sys.which(perl)
  if (perl=="" || perl=="perl")
    stop(errorMsg)
  if (.Platform$OS == "windows") {
    if (length(grep("rtools", tolower(perl))) > 0) {
      perl.ftype <- shell("ftype perl", intern = TRUE)
      if (length(grep("^perl=", perl.ftype)) > 0) {
        perl <- sub('^perl="([^"]*)".*', "\\1", perl.ftype)
      }
    }
      }
  if (verbose) cat("Using perl at", perl, "\n")
  perl
}
#.findPerlDoc <- function(perldoc, verbose = "FALSE")
#{
#  errorMsg <- "perldoc executable not found. Use perldoc= argument to specify the correct path." 
#  if (missing(perldoc))
#    {
#      perldoc = "perldoc"
#    }
#  perldoc = Sys.which(perldoc)
#  if (perldoc=="" || perldoc=="perldoc")
#    stop(errorMsg)
#  if (.Platform$OS == "windows") {
#    if (length(grep("rtools", tolower(perldoc))) > 0) {
#      perldoc.ftype <- shell("ftype perldoc", intern = TRUE)
#      if (length(grep("^perldoc=", perldoc.ftype)) > 0) {
#        perldoc <- sub('^perldoc="([^"]*)".*', "\\1", perldoc.ftype)
#      }
#    }
#      }
#  if (verbose) cat("Using perldoc at", perldoc, "\n")
#  perldoc
#}
.ClustalChecks <- function(ClustalCommand="clustalo") {
    #Check for clustalo in the PATH
    message("Checking if clustalo is in the PATH...")
    ClustalCommand <- Sys.which(ClustalCommand)
    if(ClustalCommand=="") {
        warning("Clustal Omega is not in the PATH:\nYou can either change clustalo command using lmParams function or use the web service. See ?setup")
        return()
    }
    message("Checking clustalo Version...")
    ClustalVersion <- system(paste(ClustalCommand ,  "--version") , intern=TRUE)
    if(!grepl("^1.2" , ClustalVersion))
        warning("Clustal Omega version could be not compatible. LowMACA was tested on 1.2.x")
}
.PerlModuleChecks <- function(stop=FALSE , perl="perl"){
    if(stop)
        myFunc <- stop
    else
        myFunc <- warning
    message("Checking perl installation...")
    perl <- tryCatch(.findPerl(perl=perl) , error=function(e) "no perl")
    if(perl=="no perl"){
        warning("perl executable is not in the PATH. Remember to install perl and its modules XML::Simple and LWP if you want to use web service aligner.")
        return()
    }
    message("Checking perl modules XML::Simple and LWP...")
    checkXML <- system( paste(perl , "-MXML::Simple -e 1") , intern=TRUE)
    failed <- !is.null(attr(checkXML, 'status')) && attr(checkXML, 'status') != 0
    if(failed)
        myFunc(paste( "XML::Simple module for perl is not installed. 
            If you don't want to install a local clustal omega and use the web service, XML::Simple is required", checkXML , sep="\n"))
    perl <- .findPerl()
    checkLWP <- system( paste(perl , "-MLWP -e 1") , intern=TRUE)
    failed <- !is.null(attr(checkLWP, 'status')) && attr(checkLWP, 'status') != 0
    if(failed)
        myFunc(paste( "LWP module for perl is not installed. 
            If you don't want to install a local clustal omega and use the web service, LWP is required", checkLWP , sep="\n"))
}
.myTrimmer <- function (object, ...) 
{
   s <- sub("^[\t\n\f\r ]*", "", as.character(object))
   s <- sub("[\t\n\f\r ]*$", "", s)
   s
}

showTumorType <- function() {
    mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
    all_cancer_studies <- getCancerStudies(mycgds)[,c(1,2)]
    all_cancer_studies2 <- unique(
        data.frame(
            Code=sapply(all_cancer_studies$cancer_study_id
                , function(x) strsplit(x , "_")[[1]][1])
            , Full_Name=sapply(all_cancer_studies$name
                , function(x) .myTrimmer(strsplit(x , "\\(")[[1]][1]))
                ))
    all_cancer_studies3 <- aggregate(Full_Name~Code, all_cancer_studies2
        , FUN=function(x) {paste(x , collapse="|")})
    out <- as.character(all_cancer_studies3$Code)
    names(out) <- all_cancer_studies3$Full_Name
    return(out)
}

.checkGene_to_geneID <- function(genes, myUni, myAlias) {
# This function checks the gene ids provided by the user
# transforming eventual aliases to official HugoSymbols,
# removing duplicated itmes and returning the corresponding EntrezID
    myAliasUnmapped <- myAlias[ myAlias$MappedByLowMACA=="no" , ]
    myAlias <- myAlias[ myAlias$MappedByLowMACA=="yes" , ]
    # make all symbols uppercase in order to avoid
    # ambiguities
    genes <- toupper(genes)
    # good genes are considered when provided as EntrezID or HugoSymbols
    good_genes <- c(as.character(myUni$Gene_Symbol), 
        as.character(myUni$Entrez))
    if( all(genes %in% good_genes) ) {
        message("All Gene Symbols correct!")
        # collect the annotations of id provided as gene symbols
        isGeneSymbol <- myUni$Gene_Symbol %in% genes
        Official_gs <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
        # collect the annotations of id provided as EntrezID
        isEntrez <- as.character(myUni$Entrez) %in% genes
        Official_entrez <- myUni[isEntrez, c("Gene_Symbol", "Entrez")]
        # merge the two annotations
        Official <- rbind(Official_gs, Official_entrez)
        Official$Alias <- rep(NA, nrow(Official))
        # check for duplicated itmes, in case there are make a
        # warning with the duplicated items
        if( any(duplicated(Official$Entrez)) ) 
        {
            warning("Either there were duplicated Gene Symbols or Entrez IDs 
                or you put a Gene Symbol along with its Entrez ID:"
              , immediate.=TRUE)
            print(Official[duplicated(Official$Entrez), ])
        }
        # remove duplicated items, drop the factors' levels 
        # and return
        Official <- unique(Official[, c("Gene_Symbol", "Entrez")])
        Official <- droplevels(Official)
        return(Official)
    } else {
        # collect the annotations of id provided as gene symbols
        isGeneSymbol <- myUni$Gene_Symbol %in% genes
        Official_gs <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
        # collect the annotations of id provided as EntrezID
        isEntrez <- as.character(myUni$Entrez) %in% genes
        Official_entrez <- myUni[isEntrez, c("Gene_Symbol", "Entrez") ]
        # merge the two annotations
        Official <- rbind(Official_gs, Official_entrez)
        Official$Alias <- rep(NA, nrow(Official))
        # find genes who were not provided as either EntrezID or HugoSymbol
        notOfficial <- setdiff(genes, good_genes)
        # try to assign them through alias
        bad_alias <- setdiff(notOfficial, myAlias$Alias)
        if( length(bad_alias)==0 ) 
        {
            notOfficial <- lapply(notOfficial, function(x) {
                gs <- unique(myAlias[myAlias$Alias==x,"Official_Gene_Symbol"])
                isGeneSymbol <- myUni$Gene_Symbol %in% gs
                official <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
                official$Alias <- rep(x, nrow(official))
                return(official)
                })
            if( all(sapply(notOfficial, nrow)==1)) {
                out <- do.call("rbind", notOfficial)
                message("These Genes were reverted to their official Gene Symbol:")
                print(out)
                out <- rbind(Official, out)
                if( any( duplicated(out$Entrez) ) ) {
                    warning("There were duplicated Gene Symbols or Entrez IDs 
                        or you put a Gene Symbol along with its Entrez ID:"
                      , immediate.=TRUE)
                    print(out[duplicated(out$Entrez), ])
                }
                out <- unique(out[, c("Gene_Symbol", "Entrez")])
                out <- droplevels(out)
                return(out)
                # return(droplevels( unique(out[, c(1, 2)]) ) )
            } else {
                message("There is an ambiguity with some aliases:")
                bad_alias_2 <- sapply(notOfficial, length)!=1
                bad_alias_3 <- do.call("rbind", notOfficial[bad_alias_2])
                print( bad_alias_3 )
                message("Choose a correct Gene Symbol and start over :(")
                return( bad_alias_3 )
            }
        } else {
            wrongGenes <- !( bad_alias %in% c(myAliasUnmapped$Official_Gene_Symbol ,
                                                myAliasUnmapped$Alias) )
            if( any(wrongGenes) ) {
                wrongGenes <- bad_alias[wrongGenes]
                message("There are invalid Gene Symbol or Entrez IDs:")
                print(wrongGenes)
                stop("Check manually and start over :(")                
            }
            unmappedGenes <- ( bad_alias %in% c(myAliasUnmapped$Official_Gene_Symbol ,
                                                myAliasUnmapped$Alias) )
            if( any(unmappedGenes) ) {
                unmappedGenes <- bad_alias[unmappedGenes]
                message("There are valid genes that have not been mapped by LowMACA:")
                print(unmappedGenes)
                stop("We are sorry, remove these genes and start over :(")
            }
        }
    }
}

.Trident_Score <- function(origMAlign , cons_mat="BLOSUM62", param=c(1 , 0.5 , 2) ) {
            m <- consensusMatrix(origMAlign)
            freq_mat <- pcm2pfm(m)
            aminos <- rownames(m)[rownames(m)!='-']
            if(cons_mat=="BLOSUM62") {
                data("BLOSUM62", envir = environment())
                BLOSUM62 <- get("BLOSUM62", envir  = environment())
                myBLOSUM <- BLOSUM62[aminos , aminos]
            } else {
                myBLOSUM <- cons_mat[aminos , aminos]
            }
            ### force negative elements on the diagonal to be 0
            if(any(diag(myBLOSUM)<0)) {
                message('There are some negative elements in the diagonal elements if consensus matrix')
                print(diag(myBLOSUM)[diag(myBLOSUM)<0])
                message('Forcing them to be zero')
                diag(myBLOSUM)[diag(myBLOSUM)<0] <- 0
            }
            myBLOSUM_valdar <- sapply(aminos , function(col) {
                sapply(aminos , function(row) {
                        myBLOSUM[col , row]/sqrt(myBLOSUM[col , col]*myBLOSUM[row , row])
                    })
                })
            #first member
            lambda_t <- 1/log2(min(21 , nrow(origMAlign)))
            log2freq_mat <- ifelse( is.infinite(log2(freq_mat)) , 0 , log2(freq_mat))
            t <- lambda_t*(-colSums(freq_mat*log2freq_mat))
            #second member
            MyNorm <- function(x) sqrt(sum(x^2))
            lambda_r <- 1/( sqrt(20)*(max(myBLOSUM_valdar) - min(myBLOSUM_valdar)) )
            r <- sapply(1:ncol(freq_mat) , function(i) {
                pos <- freq_mat[ aminos , i]
                amin <- names(pos[pos!=0])
                if(length(amin)==1) {
                    return(0)
                } else {
                    myBLOSUM_cut <- myBLOSUM_valdar[ , amin]
                    x_mean <- (1/length(amin))*rowSums(myBLOSUM_cut)
                    distance <- sum(sapply(1:ncol(myBLOSUM_cut) , function(x) {
                                                eucl <- MyNorm(x_mean-myBLOSUM_cut[ , x])
                                    }))
                    out <- distance/length(amin)
                    out <- out*lambda_r
                    return(out)
                }
            })
            #third member
            if( "-" %in% rownames(freq_mat) )
                g <- freq_mat["-" , ]
            else
                g <- 0
            #Finally
            result <- ((1-t)^param[1])*((1-r)^param[2])*((1-g)^param[3])
            return(result)
}

.clustalOAlign <- function(genesData, clustal_cmd, filename , mail , perlCommand)
{
    # # get the protein sequences corresponding to the selected genes
    # # from the uniprot_file
    # arrange the protein sequences to make a fasta format file
    fasta <- c()
    seq_names <- rownames(genesData)
    if(length(seq_names)>2000 && !is.null(mail))
        stop("You cannot evaluate more than 2000 sequences in web mode. Install a local clustalo")
    for (i in 1:nrow(genesData)) {
        seq_name <- paste( ">", seq_names[i], sep='')
        seq <- as.character(genesData[i, 'AMINO_SEQ'])
        fasta <- c(fasta, seq_name, seq)
    }

    # write the fasta file on a temporary file
    fastafile <- tempfile()
    write.table(fasta, fastafile, 
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    clustalout_clu <- tempfile()
    # run clustaomega output in Clustal format
    if(nrow(genesData)>2) {
        if(is.null(mail)) {
            #clustalout_clu <- tempfile()
            ClustalOmega <- Sys.which(clustal_cmd)
            if(ClustalOmega=="")
                stop("Clustal Omega command not found. clustalo is not in your PATH or it was not installed")
                #clustal omega pairwise distance matrix
            dist_mat <- tempfile()
            dist_cmd <- paste("--distmat-out" , dist_mat , sep="=")
            exec <- paste(ClustalOmega, '--outfmt=clu' , '--percent-id' 
                , dist_cmd , "--full --force" , '-i', fastafile, '>', clustalout_clu)
            #Windows doesn't accept the redirection >, so we must use shell
            if( Sys.info()['sysname']=="Windows" ) {
                exec <- gsub("\\\\" , "/" , exec)
                cmdCheck <- shell(exec , intern=TRUE)
            } else {
                cmdCheck <- system(exec , intern=TRUE)
            }
            #Check if system call to clustalo had 0 exit status
            failed <- !is.null(attr(cmdCheck, 'status')) && attr(cmdCheck, 'status') != 0
            if(failed)
                stop(paste("Alignment with ClustalOmega had non 0 exit status:",cmdCheck , sep="\n"))
            score <- .scoreMatrix(dist_mat , mail=mail)
            aln <- .clustalMatrix(clustalout_clu)
            see_aln <- readAAMultipleAlignment(filepath = clustalout_clu
                , format="clustal")
            if( is.null(filename) ) {
                unlink(clustalout_clu)
            } else {
                file.rename(clustalout_clu, filename)
            }
            unlink(dist_mat)
            unlink(fastafile)
        } else {
            perl <- .findPerl(perl=perlCommand)
            package.dir <- system.file(package='LowMACA')
            script <- file.path(package.dir,'clustalo_lwp.pl')
            mailArgument <- paste("--email" , mail)
            webOut <- file.path(tempdir() , "webClustal")
            webOutArgument <- paste("--outfile" , shQuote(webOut))
            exec <- paste(shQuote(perl)
                            , shQuote(script)
                            , "--noguidetreeout --stype protein --dismatout --outfmt clustal"
                            , mailArgument
                            , webOutArgument
                            , shQuote(fastafile)
                    )
            cmdCheck <- system(exec , intern=TRUE)
            failed <- !is.null(attr(cmdCheck, 'status')) && attr(cmdCheck, 'status') != 0
            if(failed)
                stop(paste("Alignment with ClustalOmega had non 0 exit status:",cmdCheck , sep="\n"))
            dist_mat <- paste0(webOut , ".pim.pim")
            clustalout_clu <- paste0(webOut , ".aln-clustal.clustal")
            score <- .scoreMatrix(dist_mat , mail=mail)
            aln <- .clustalMatrix(clustalout_clu)
            see_aln <- readAAMultipleAlignment(filepath = clustalout_clu
                , format="clustal")
            if( is.null(filename) ) {
                unlink(clustalout_clu)
            } else {
                file.rename(clustalout_clu, filename)
            }
            unlink(dist_mat)
            unlink(fastafile)
        }
    } else {
        if(nrow(genesData)==1) {
                # Fasta length
            len <- nchar(as.character(genesData[1, 'AMINO_SEQ']))
            aln <- data.frame(Gene_Symbol=rep(genesData[, 'Gene_Symbol'], len)
                           , Protein=rep(genesData[, 'UNIPROT'])
                           , Entrez=rep(genesData[, 'Entrez'])
                           , Align=1:len
                           , Ref=1:len
                        )
            see_aln <- AAMultipleAlignment(genesData[1, 'AMINO_SEQ'])
            score <- NA
                # Fasta length
            aln <- data.frame(Gene_Symbol=rep(genesData[ , 'Gene_Symbol'] , len)
                , domainID=rep(genesData[ , 'Pfam_ID'])
                , Entrez=rep(genesData[ , 'Entrez'])
                , Envelope_Start=rep(genesData[ , 'Envelope_Start'])
                , Envelope_End=rep(genesData[ , 'Envelope_End'] , len)
                , Align=1:len
                , Ref=1:len
                )
            score <- NA
        } else {
            if(is.null(mail)) {
            # With a 2 sequence alignment we do not create a distance matrix
                ClustalOmega <- Sys.which(clustal_cmd)
                if(ClustalOmega=="")
                    stop("Clustal Omega command not found. clustalo is not in your PATH or it was not installed")
                warning("There are less than 3 sequences aligned, so no distance matrix can be calculated" 
                    , immediate.=TRUE)
                exec <- paste(ClustalOmega, '--outfmt=clu' , '-i', fastafile, '>', clustalout_clu)
                #Windows doesn't accept the redirection with system >, so we must use shell
                if( Sys.info()['sysname']=="Windows" ) {
                    exec <- gsub("\\\\" , "/" , exec)
                    cmdCheck <- shell(exec , intern=TRUE)
                } else {
                    cmdCheck <- system(exec , intern=TRUE)
                }
                failed <- !is.null(attr(cmdCheck, 'status')) && attr(cmdCheck, 'status') != 0
                if(failed)
                    stop(paste("Alignment with ClustalOmega had non 0 exit status:",cmdCheck , sep="\n"))                
                score <- "It is not possible to calculate distance matrix with less than 3 sequences"
                aln <- .clustalMatrix(clustalout_clu)
                see_aln <- readAAMultipleAlignment(filepath = clustalout_clu
                    , format="clustal")
                if( is.null(filename) ) {
                    unlink(clustalout_clu)
                } else {
                    file.rename(clustalout_clu, filename)
                }
                unlink(fastafile)
            } else {
                perl <- .findPerl(perl=perlCommand)
                package.dir <- system.file(package='LowMACA')
                script <- file.path(package.dir,'clustalo_lwp.pl')
                mailArgument <- paste("--email" , mail)
                webOut <- file.path(tempdir() , "webClustal")
                webOutArgument <- paste("--outfile" , shQuote(webOut))
                exec <- paste(shQuote(perl)
                            , shQuote(script)
                            , "--noguidetreeout --stype protein --nodismatout --outfmt clustal"
                            , mailArgument
                            , webOutArgument
                            , shQuote(fastafile)
                    )
                cmdCheck <- system(exec , intern=TRUE)
                failed <- !is.null(attr(cmdCheck, 'status')) && attr(cmdCheck, 'status') != 0
                if(failed)
                    stop(paste("Alignment with ClustalOmega had non 0 exit status:",cmdCheck , sep="\n"))
                clustalout_clu <- paste0(webOut , ".aln-clustal.clustal")
                score <- "It is not possible to calculate distance matrix with less than 3 sequences"
                aln <- .clustalMatrix(clustalout_clu)
                see_aln <- readAAMultipleAlignment(filepath = clustalout_clu
                    , format="clustal")
                if( is.null(filename) ) {
                    unlink(clustalout_clu)
                } else {
                    file.rename(clustalout_clu, filename)
                }
                unlink(fastafile)
            }
        }
    }
    return(list(ALIGNMENT=aln, SCORE=score, CLUSTAL=see_aln))
}

.clustalMatrix <- function(filename)
{
    origMAlign <- readAAMultipleAlignment(filepath = filename , format="clustal")
    origMAlign_mat <- Biostrings::as.matrix(origMAlign)
    out <- apply(origMAlign_mat, 1 , function(seq) {
        match <- rep(NA, length(seq))
        ix <- which(seq != '-')
        match[ix] <- 1:length(ix)
        return(match)
        })
    colnames(out) <- rownames(origMAlign_mat)
    out <- as.data.frame(out)
    out$Align <- 1:nrow(out)
    out <- melt(out, id.vars='Align')
    colnames(out)[2:3] <- c('Seq','Ref')
    domainID <- out$Seq
    lines_name <- do.call('rbind', 
        sapply(as.character(out$Seq), strsplit, '\\|'))
    rownames(lines_name) <- NULL
    envelope <- do.call('rbind', sapply(lines_name[,4], strsplit, ';'))
    envelope <- apply(envelope, 2, as.numeric)
    lines_name <- data.frame(domainID, lines_name[,1:3], envelope)
    colnames(lines_name) <- c('domainID', 'Gene_Symbol', 'Pfam', 'Entrez', 
        'Envelope_Start', 'Envelope_End')
    out <- cbind(lines_name, out[, c('Align', 'Ref')])
    out$Ref <- out$Ref + out$Envelope_Start - 1
    return(out) 
}

    #Calculate similarity matrix and some summary statistics
.scoreMatrix <- function(filename , mail)
{
    dist_seq <- read.table(filename , row.names=1 , skip=1)
    if(!is.null(mail)) {
        rownames(dist_seq) <- dist_seq$V2
        dist_seq$V2 <- NULL
    }
    colnames(dist_seq) <- rownames(dist_seq)
    dist_seq <- as.matrix(dist_seq)
    diag(dist_seq) <- NA
    mean_dist <- rowMeans(dist_seq , na.rm=TRUE)
    median_dist <- apply(dist_seq , 1 , function(x) median(x , na.rm=TRUE) )
    max_dist <- apply(dist_seq , 1 , function(x) max(x , na.rm=TRUE) )
    min_dist <- apply(dist_seq , 1 , function(x) min(x , na.rm=TRUE) )
    summary <- data.frame(MEAN_SIMILARITY=mean_dist
                            , MEDIAN_SIMILARITY=median_dist
                            , MAX_SIMILARITY=max_dist
                            , MIN_SIMILARITY=min_dist
                            )
    rownames(summary) <- rownames(dist_seq)
    if( any(summary$MAX_SIMILARITY <= 20)) {
        warning("There are sequences very dissimilar to the others! 
            Consider to exclude them because the maximum similarity with any other sequence is less than 20%" , immediate.=TRUE)
        print(summary[ summary$MAX_SIMILARITY <= 20 , ])
    }
    return(list(DIST_MAT=dist_seq , SUMMARY_SCORE=summary))
}

.getGeneMutations <- function(myGenes=myGenes
                            ,parallelize=FALSE  
                            ,mutation_type=c("missense", "all", "truncating" , "silent") 
                            ,NoSilent=TRUE 
                            ,Exonic=TRUE
                            ,tumor_type="all"
                            )
{
    mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
    all_cancer_studies <- getCancerStudies(mycgds)[,c(1,2)]
    mutation_type <- mutation_type[1]
        # If I want just silent mutation, this overwrite NoSilent mode
    if(mutation_type=="silent") NoSilent=FALSE
    if(tumor_type[1]=="all") {
        chosenTumors <- all_cancer_studies[,1]
    } else {
        chosenTumors <- all_cancer_studies[grepl( paste(tumor_type , collapse="|") , all_cancer_studies[,1] , ignore.case=TRUE) , 1]
    }
        #Remove cell_lines and npc tumor that generate problems
    #chosenTumors <- chosenTumors[ !(chosenTumors %in% c("cellline_nci60" , "cellline_ccle_broad" , "npc_nusingapore")) ]
    if(parallelize) {
        applyfun <- mclapply
        options('mc.cores'=detectCores())
    } else applyfun <- lapply
    out_double <- applyfun(chosenTumors , function(i)
        {
            geneticProfile <- getGeneticProfiles(mycgds, i)[ ,c(1:2)]
            sel <- geneticProfile$genetic_profile_name=="Mutations"
            geneticProfile <- geneticProfile[sel, 1]
            caseList <- getCaseLists(mycgds, i)
            sel <- caseList$case_list_name=="Sequenced Tumors"
            if(any(sel)) {
                caseListID <- caseList[sel, 1]
                tryCatch(
                    muts <- getMutationData( mycgds 
                        , caseList=caseListID 
                        , geneticProfile=geneticProfile 
                        , genes=myGenes[ , 'Entrez'])
                    , error=function(e) message(paste("Impossible to retrive mutations from" , i , "study"))
                )
                if(!exists("muts")) {
                    muts <- NULL
                    patients <- NULL
                } else {
                    if(ncol(muts)!=22) {
                        muts <- NULL
                        patients <- NULL
                    } else {
                        patients <- strsplit(caseList[sel, 'case_ids'] , split=" ")[[1]]
                    }
                }
            } else {
                muts <- NULL
                patients <- NULL
            }
            return( list( out=muts , patients=patients) )
        })
    if(all(sapply(out_double , function(x) is.null(x$out)))) {
        message("There are no mutations available for the selected tumor types and genes")
        return(list( Mutations=NULL , AbsFreq=NULL ))
    }
        #These are the number of samples
    samples_out <- lapply(1:length(out_double) , function(x) out_double[[x]][["patients"]])
    names(samples_out) <- chosenTumors
        #Create a set of samples per tumor type and not per study type (ex. brca_tcga and brca_pub will be joined)
    chosenTumors_type <- unique(sapply(chosenTumors, function(x) strsplit(x , split="_")[[1]][1]))
    npat_type <- c()
    for (i in chosenTumors_type) {
        selected_tum <- samples_out[grep(i , names(samples_out) , value=TRUE)]
        selected_pat <- unique(unlist(selected_tum))
        npat_type <- c(npat_type , length(selected_pat) )
    }
    names(npat_type) <- chosenTumors_type  
        #Mutation Dataset
    mut <- as.data.frame( 
                rbindlist(
                    lapply(1:length(out_double) , function(x) out_double[[x]][["out"]]) 
                    ))
    mut$Tumor_Type <- sapply(mut$genetic_profile_id , function(x) strsplit(x , split="_")[[1]][1])
    goodCols <- c("entrez_gene_id"
                ,"gene_symbol"
                ,"case_id"
                ,"mutation_type"
                ,"amino_acid_change"
                ,"Tumor_Type"
                )
    mut <- mut[ , goodCols]
        #Delete all the splice sites outside the coding region (reported as e1-2 or similar notations)
    mut <- mut[ !grepl("^e" , mut$amino_acid_change) , ]
    mut$letter <- substr(mut$amino_acid_change,1,1)
    mut$position <- as.numeric(as.character(str_extract(
                    string=mut$amino_acid_change,pattern="\\d+")))
    mut <- data.frame(
        Entrez=mut$entrez_gene_id
        , Gene_Symbol=mut$gene_symbol
        , Amino_Acid_Letter=mut$letter
        , Amino_Acid_Position=mut$position
        , Amino_Acid_Change=mut$amino_acid_change
        , Mutation_Type=mut$mutation_type
        , Sample=mut$case_id
        , Tumor_Type=mut$Tumor_Type)
    for (i in colnames(mut)){
        if(class(mut[,i])=="factor")
            mut[,i] <- as.character(mut[,i])
    }
    if(NoSilent) {
        mut <- mut[ mut$Mutation_Type!="Silent" , ]
    }
    if(Exonic) {
        notTransc <- c("3'UTR"
                    ,"3'Flank"
                    ,"5'UTR"
                    ,"5'Flank"
                    ,"IGR1"
                    ,"IGR"
                    ,"Intron"
                    ,"RNA"
                    ,"Targeted_Region"
                    )
        mut <- mut[ !(mut$Mutation_Type %in% notTransc) , ] 
    }
    if( mutation_type=="missense" ) {
        missense <- c("Missense_Mutation"
                    ,"In_Frame_Del"
                    ,"In_Frame_Ins"
                    )
        mut <- mut[ mut$Mutation_Type %in% missense , ]
    }
    if( mutation_type=="silent" ) {
        mut <- mut[ mut$Mutation_Type=="Silent" , ]
    }
    if( mutation_type=="truncating" ) {
        truncating <- c("Frame_Shift_Del"
                        ,"Nonsense_Mutation"
                        ,"Translation_Start_Site"
                        ,"Frame_Shift_Ins"
                        ,"Nonstop_Mutation"
                        ,"Splice_Site"
                        ,"Indel"
                        )
        mut <- mut[ mut$Mutation_Type %in% truncating , ]
    }
    mut <- unique(mut)
        #Create frequency table by gene and tumor type
    tum_type <- factor(mut$Tumor_Type , levels=chosenTumors_type)
    gene <- tryCatch( factor(mut$Gene_Symbol , levels=levels(myGenes[ , "Gene_Symbol"])) 
                        , error=function(e) factor(mut$Gene_Symbol , levels=unique(mut$Gene_Symbol))
                        , warning=function(w) factor(mut$Gene_Symbol , levels=unique(mut$Gene_Symbol))
                        )
    sample_table <- as.data.frame.matrix(table(tum_type , gene ))
    sample_table2 <- data.frame(StudyID = chosenTumors_type
                                    ,Total_Sequenced_Samples = unname(npat_type)
                                    )
    sample_table_absFreq <- merge( sample_table2 , sample_table , by.x="StudyID" , by.y="row.names")
    return(list( Mutations=mut , AbsFreq=sample_table_absFreq ))
}

.getLocalGeneMutations <- function(myGenes=myGenes
                            ,localData=NULL  
                            ,mutation_type=c("missense", "all", "truncating" , "silent") 
                            ,NoSilent=TRUE 
                            ,Exonic=TRUE
                            ,tumor_type="all"
                            )
{
    if( is.null(localData) ) stop('no local file provided')
    ## all mutatios from local data
    mut <- localData
    ## filter: genes
    chosenGenes <- myGenes$Gene_Symbol
    mut <- mut[mut$Gene_Symbol %in% chosenGenes, ]
    ## filter: mutation type
    mutation_user_choiche <- mutation_type[1]
    chosenMutations <- unique(mut$Mutation_Type)
    ## check flags
    if(mutation_user_choiche=="silent") NoSilent=FALSE
    if(NoSilent)
        chosenMutations <- chosenMutations[chosenMutations != "Silent"]
    if(Exonic) {
        notTransc <- c("3'UTR", "3'Flank", "5'UTR", "5'Flank"
                    ,"IGR1", "IGR", "Intron", "RNA", "Targeted_Region")
        chosenMutations <- chosenMutations[!chosenMutations %in% notTransc]
    }
    ## check user choiche
    if( mutation_user_choiche=="missense" ) {
        chosenMutations <- chosenMutations[chosenMutations %in% 
            c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")]
    } else if ( mutation_user_choiche=="silent" ) {
        chosenMutations <- chosenMutations[chosenMutations == 'Silent']
    } else if( mutation_user_choiche=="truncating" ) {
        truncating <- c("Frame_Shift_Del"
                        ,"Nonsense_Mutation"
                        ,"Translation_Start_Site"
                        ,"Frame_Shift_Ins"
                        ,"Nonstop_Mutation"
                        ,"Splice_Site"
                        ,"Indel"
                        )
        chosenMutations <- chosenMutations[chosenMutations %in% truncating]
    }
    mut <- mut[mut$Mutation_Type %in% chosenMutations, ]
    ## filter: tumor types
    chosenTumors <- unique(mut$Tumor_Type)
    if( tumor_type[1]!="all" )
        chosenTumors <- chosenTumors[chosenTumors %in% tumor_type]
    mut <- mut[mut$Tumor_Type %in% chosenTumors, ]
    mut <- unique(mut)
    return(list( Mutations=mut , AbsFreq=NA ))
}

.alnWeights <- function(aln)
{
    aln_agg <- aln$ALIGNMENT
    aln_agg$pos_existance <- ifelse(is.na(aln_agg$Ref) , 0 , 1)
    aln_agg2 <- aggregate(pos_existance~Align , data=aln_agg , FUN=sum , simplify=TRUE)
    return(aln_agg2$pos_existance/sum(aln_agg2$pos_existance))
}

.MAD <- function(x) {
    med <- median(x , na.rm=TRUE)
    MAD <- median( abs(x - med) )
    return(1.4826 * MAD)
}

######################
###### calculate entropy
############################

.shannon <- function(q) {
    diff <- diff(q$x)[1]
    p <- q$y[q$y != 0]
    shan <- -sum(p*log(p))*diff
    return(shan)
}

.profileDensity <- function(profile, bw=NULL)
{
    nPos <- length(profile)
    positions <- which(as.logical(profile))
    positions <- rep(positions, times=profile[positions])
    if( is.null(bw) ) {
        d <- density(positions, from=1, to=nPos, n=nPos)
    } else {
        if( bw==0 ) {
            d <- list(x=1:nPos, y=profile, bw=0)
        } else {
            d <- density(positions, bw=bw, from=1, to=nPos, n=nPos)
        }
    }
    # normalize before output
    d$y <- d$y/sum(d$y)
    return(d)
}

.profileEntropy <- function(profile, bw=NULL, norm=TRUE, model=NULL, weights=NULL, ...) 
{
    d <- .profileDensity(profile, bw=bw, ...)
    ent <- .shannon(d)
    if( is.null(bw) ) bw <- d$bw
    if( norm ) {
        if( !is.null(model) ) {
            unif <- model(sum(profile))
        } else {
            len <- length(profile)
            nmut <- sum(profile)
            unif <- .sampleUnifEntropyL(len, nmut, bw=bw, weights=weights)
        }
        mean <- unif[[3]]-unif[[1]]
        var <- unif[[2]]^2
        ## check: if variance is 0, put the
        ## profile pvalue to zero
        if( var==0 ) {
            pval <- 1
        } else {
            shape <- mean^2/var
            scale <- var/mean
            pval <- pgamma(unif[[3]]-ent, shape=shape, scale=scale, lower.tail=FALSE)            
        }
        return(log10(pval))
    } else return(ent)
}

.sampleUnifEntropyL.old <- function(len, nmut, bw, nboot=1000, weights=NULL, 
    center=median, variability=.MAD)
{
    if(is.null(weights)) weights <- rep(1/len , len)
    boots <- sapply(1:nboot, function(i) {
        d <- density(sample(1:len, nmut, replace=TRUE , prob=weights), bw=bw, from=1, to=len, n=len)
        .shannon(d)
        })
    return(list(center=center(boots), variability=variability(boots), max=max(boots)))
}

.sampleUnifEntropyL <- function(len, nmut, bw, nboot=1000, weights=NULL, 
    center=median, variability=.MAD)
{
    if(is.null(weights)) weights <- rep(1/len , len)
    boots <- sapply(1:nboot, function(i) {
        positions <- sample(1:len, nmut, replace=TRUE , prob=weights)
        t <- table(positions)
        profile <- rep(0, len)
        profile[as.numeric(names(t))] <- t
        .profileEntropy(profile, bw=bw, norm=FALSE)
        })
    return(list(center=center(boots), variability=variability(boots), max=max(boots)))
}

.makeUniformModel <- function(mat, bw, nboot=1000, plotOUT=TRUE, 
    weights=NULL, center=median, variability=.MAD, parallelize=TRUE ) 
{
    if( parallelize ) {
        applyfun <- mclapply
    } else applyfun <- lapply
    geneLen <- ncol(mat)
    if( is.null(weights) ) weights <- rep(1/geneLen , geneLen)
    # minNMut <- min(apply(mat, 1, sum))
    # maxNMut <- sum(mat)
    minNMut <- floor(sum(mat)/10)*10 #round to the upper ten
    maxNMut <- ceiling(sum(mat)/10)*10 #round to the upper ten
    # nMutInt <- round(seq(from=minNMut,to=maxNMut, length.out=100))
    nMutInt <- unique(c(minNMut, maxNMut))
    if(length(nMutInt)==1) nMutInt <- c(nMutInt , nMutInt+1)
    outReal <- applyfun(nMutInt, function(i) 
        .sampleUnifEntropyL(geneLen, i, bw=bw, nboot=nboot , weights=weights,
            center=center, variability=variability))
    outReal <- do.call('cbind',outReal)
    polynomialModel <- function(x, par) {
        sapply( x
            , function(x_i)
                sum(sapply(1:length(par), function(i) 
                    if(is.na(par[i])) 0 else x_i^(i-1) * par[i] )))
    }
    pn.optim.aic <- function( tpts , experiment, variance=NULL ) {
        if( length(experiment) < 2 ) return(NA)
        polyOrderChisq <- function(i) {
            model <- lm( experiment~poly( tpts , i , raw=TRUE ) )
            return(list(par=model$coeff, value=AIC(model)))
        }
        sapply(1:min(30,length(tpts)-1), polyOrderChisq)
    }
    pnout <- pn.optim.aic(nMutInt, unlist(outReal[1,]), 1)
    degree <- min(which.min(unlist(pnout[2,])))
    par.mean <- pnout[1,degree]$par
    model.mean <- function(mut) polynomialModel(mut, par.mean)
    pnout <- pn.optim.aic(nMutInt, unlist(outReal[2,]), 1)
    degree <- min(which.min(unlist(pnout[2,])))
    par.sd <- pnout[1,degree]$par
    model.sd <- function(mut) polynomialModel(mut, par.sd)
    pnout <- pn.optim.aic(nMutInt, unlist(outReal[3,]), 1)
    degree <- min(which.min(unlist(pnout[2,])))
    par.max <- pnout[1,degree]$par
    model.max <- function(mut) polynomialModel(mut, par.max)
    modelUnif <- function(nmut) 
        list(mean=model.mean(nmut), sd=model.sd(nmut), max=model.max(nmut))
    if( plotOUT ) {
        par(mfrow=c(1,3))
        plot(nMutInt, unlist(outReal[1,]), xlab='n of mutations', ylab='',
            main='entropy center measure')
        lines(nMutInt, model.mean(nMutInt), col='red', lwd=3)
        plot(nMutInt, unlist(outReal[2,]), xlab='n of mutations', ylab='',
            main='entropy variability measure')
        lines(nMutInt, model.sd(nMutInt), col='red', lwd=3)     
        plot(nMutInt, unlist(outReal[3,]), xlab='n of mutations', ylab='',
            main='max entropy measure')
        lines(nMutInt, model.max(nMutInt), col='red', lwd=3)     
    }
    return(modelUnif)
}

.makeNullProfile <- function(mat, bw, nboot=1000, plotOUT=TRUE, 
    weights=NULL, center=median, variability=.MAD) 
{
    geneLen <- ncol(mat)
    if( is.null(weights) ) weights <- rep(1/geneLen , geneLen)
    nMut <- sum(mat)
    if( bw==0 ) {
        # to prevent having center and variability 
        # that both equals 0
        center <- mean
        variability <- sd
    }
    boots <- lapply(1:nboot, function(i) {
        # density(sample(1:geneLen, nMut, replace=TRUE , prob=weights), bw=bw, from=1, to=geneLen, n=geneLen)
        positions <- sample(1:geneLen, nMut, replace=TRUE , prob=weights)
        t <- table(positions)
        profile <- rep(0, geneLen)
        profile[as.numeric(names(t))] <- t
        .profileDensity(profile, bw=bw)
        })
    nullDensities <- sapply(boots, '[[', 'y')
    # calulate parameters for the gamma distribution
    mu <- apply(nullDensities, 1, center)
    s <- apply(nullDensities, 1, variability)
    s2 <- s^2
    # apply gamma distribution to find thresholds
    upperThreshold <- qgamma(.95, shape=mu^2/s2, scale=s2/mu)
    lowerThreshold <- qgamma(.05, shape=mu^2/s2, scale=s2/mu)

    # pvalue of every aa
    d <- .profileDensity(colSums(mat), bw=bw) #, from=1, to=geneLen, n=geneLen)
    pvals <- pgamma(d$y, shape=mu^2/s2, scale=s2/mu, lower.tail=FALSE)
    

    return(data.frame(
        mean=mu, lTsh=lowerThreshold, uTsh=upperThreshold, profile=d$y, pvalue=pvals#, qvalue=qvals
        ))

}
