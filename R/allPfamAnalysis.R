allPfamAnalysis <- function(repos
							, allLowMACAObjects=NULL
							, mutation_type=c("missense", "all", "truncating" , "silent") 
							, NoSilent=TRUE
							, mail=NULL
							, perlCommand="perl"
							, verbose=FALSE
							, conservation=0.1
							, use_hmm=FALSE
							, datum=FALSE
							, clustal_cmd="clustalo"
							, BPPARAM=bpparam("SerialParam")
                            ) {
	#library(LowMACA)
	# analyze all the pfams involved in the repos
	mutation_type <- mutation_type[1]
	if(mutation_type=="silent") NoSilent=FALSE
	#reading the file if it is not a variable
	if(is.character(repos))
		repos <- read.table(repos , sep="\t" , header=TRUE , as.is=TRUE)
	#Delete all non-SNVs mutation and all non-TCGA MutationType
    repos <- repos[ !is.na(repos$Mutation_Type) , ]
    bad_mut_types <- c("Fusion" , "COMPLEX_INDEL" , "vIII deletion" , "Splice_Site_SNP" , "Indel")
    repos <- repos[ !(repos$Mutation_Type %in% bad_mut_types) , ]
    repos <- repos[ !(repos$Amino_Acid_Change=="MUTATED") , ]
    #repos <- repos[ !grepl("^e" , repos$Amino_Acid_Change) , ]
    ########################
	# subsetting mutations
	######################
	if(NoSilent) {
        repos <- repos[ repos$Mutation_Type!="Silent" , ]
    }
	notTransc <- c("3'Flank"
                    ,"5'Flank"
                    ,"IGR1"
                    ,"IGR"
                    ,"Intron"
                    ,"RNA"
                    ,"Targeted_Region"
                    )
    repos <- repos[ !(repos$Mutation_Type %in% notTransc) , ]
    if( mutation_type=="missense" ) {
        missense <- c("Missense_Mutation"
                    ,"In_Frame_Del"
                    ,"In_Frame_Ins"
                    )
        repos <- repos[ repos$Mutation_Type %in% missense , ]
    }
    if( mutation_type=="silent" ) {
        repos <- repos[ repos$Mutation_Type=="Silent" , ]
    }
    if( mutation_type=="truncating" ) {
        truncating <- c("Frame_Shift_Del"
                        ,"Nonsense_Mutation"
                        ,"Translation_Start_Site"
                        ,"Frame_Shift_Ins"
                        ,"Nonstop_Mutation"
                        ,"Splice_Site"
                        ,"Indel"
                        ,"5'UTR"
                        ,"3'UTR"
                        )
        repos <- repos[ repos$Mutation_Type %in% truncating , ]
    }
	myPfam <- getMyPfam()
	allReposGenes <- unique(repos$Gene_Symbol)
	myPfamRepos <- myPfam[myPfam$Gene_Symbol %in% allReposGenes,]
	if(verbose)
		message("Splitting dataset into Pfams...")
	myPfamReposSplit <- split(myPfamRepos, myPfamRepos$Pfam_ID)
	####################
	# Parallel options (now exposed to user via BPPARAM)
	#####################
	# if(is.logical(parallel) && parallel==TRUE)
	# 	cores <- parallel::detectCores()
	# if(is.logical(parallel) && parallel==FALSE)
	# 	cores <- 1L
	# if(is.numeric(parallel))
	# 	cores <- min(parallel , parallel::detectCores())
	# if( cores > 1 ) {
	# 	applyfun <- bplapply
	# 	if( Sys.info()[['sysname']] == 'Windows' ){
	# 		#applyfun <- lapply
	# 		options(MulticoreParam=SnowParam(workers=cores))
	# 	} else {
	# 		options(MulticoreParam=MulticoreParam(workers=cores))
	# 	}
	# } else {
	# 	applyfun <- lapply
	# }
	##########################################
	# Create a LowMACA objects for every Pfam
	###########################################
	if(verbose)
		message("Creating LowMACA objects for every Pfam...")
	allPfamsLM <- bplapply(myPfamReposSplit, function(x) {
		if(Sys.info()[['sysname']] == 'Windows') library(LowMACA)
		#attach(loadNamespace(sub('package:' , '' , search()[search()!=".GlobalEnv"])))
		pfam <- unique(x$Pfam_ID)
		genes <- unique(x$Gene_Symbol)
		if(verbose)
		    message(paste('Analyzing PFAM ', pfam, '...', sep=''))
	    # if( cores>1 ) {
	    invisible(capture.output(suppressWarnings(suppressMessages({
				lmDomain <- newLowMACA(pfam=pfam, genes=genes)
				lmParams(lmDomain)$clustal_cmd <- clustal_cmd
				lmDomain <- setup(lmDomain, repos=repos 
					, mail=mail , perlCommand=perlCommand
					, use_hmm=use_hmm, datum=datum)
				lmDomain <- entropy(lmDomain)
				}))))
	    	return(lmDomain)
    # 	} else {
	   #  	suppressWarnings(suppressMessages({
				# lmDomain <- newLowMACA(pfam=pfam, genes=genes)
				# lmParams(lmDomain)$clustal_cmd <- clustal_cmd
				# lmDomain <- setup(lmDomain, repos=repos 
				# 	, mail=mail , perlCommand=perlCommand
				# 	, use_hmm=use_hmm, datum=datum)
				# lmDomain <- entropy(lmDomain)
				# }))
	   #  	return(lmDomain)
    	# }
    	} , BPPARAM=BPPARAM)
	##########################################################
	# extract gene level information on no alignment analysis
	###########################################################
	if(verbose)
		message("Performing Single Gene analysis...")
	singleSequence <- lapply(allPfamsLM , function(object) {
		if(Sys.info()[['sysname']] == 'Windows') library(LowMACA)
		#attach(loadNamespace(sub('package:' , '' , search()[search()!=".GlobalEnv"])))
		if(!verbose)
			suppressMessages(lfmSingleSequence(object, metric='qvalue', threshold=.05
				, conservation=conservation , mail=NULL , perlCommand="perl" , BPPARAM=BPPARAM))
		else
			lfmSingleSequence(object, metric='qvalue', threshold=.05
				, conservation=conservation , mail=NULL , perlCommand="perl" , BPPARAM=BPPARAM
				, verbose=TRUE)
		# } , BPPARAM=BPPARAM)
		})
	singleSequence <- do.call("rbind" , singleSequence)
	if( !is.null(allLowMACAObjects) )
		save(allPfamsLM, file=allLowMACAObjects)
	if(verbose)
		message("Performing binomial correction...")
	#########################################################
	# extract gene level information on alignment analysis
	########################################################
	pfamAnalysis <- bplapply(allPfamsLM, function(lmDomain) {
		if(Sys.info()[['sysname']] == 'Windows') library(LowMACA)
		## get the significant positions
		pfam <- unique(as.character(lmDomain@arguments$input$Pfam_ID))
		signPos <- lfm(lmDomain)
		if( is.null(signPos) || nrow(signPos)==0 ) return(NULL)
		signPosGene <- split(signPos, signPos$Gene_Symbol)
		## get all mutations that fall on genes pfam
		mutAligned <- lmDomain@mutations$aligned
		splitfun <- sapply(
			strsplit(rownames(mutAligned),'\\|')
			, '[', 1)
		nAllMutGene <- tapply(rowSums(mutAligned), splitfun, sum)
		if( length(nAllMutGene)==1 ) names(nAllMutGene) <- names(signPosGene)
		nAllMut <- sum(nAllMutGene)
		## compute p to be used in the binom test for
		## each significant mutation
		nBackgroundMutDomain <- nAllMut-nrow(signPos)
		nForegroundMutDomain <- sapply(
			split(signPos, signPos$Multiple_Aln_pos)
			, nrow
			)
		pBinom <- nForegroundMutDomain/
			(nForegroundMutDomain+nBackgroundMutDomain)
		## genes to be evaluated
		genesToBeEvaluated <- names(signPosGene)
		## first gene
		pfamOut <- lapply(genesToBeEvaluated, function(i) {
			signPosGene_i <- signPosGene[[i]]
			nAllMutGene_i <- nAllMutGene[[i]]
			nBackgroundMutGene_i <- nAllMutGene_i-nrow(signPosGene_i)
			positionPval <- sapply(
				split(signPosGene_i, signPosGene_i$Multiple_Aln_pos)
				, function(signPosGene_i_MultipleAlnPos_j) {
					nForegroundMutGene_i <- nrow(signPosGene_i_MultipleAlnPos_j)
					Multiple_Aln_pos <- as.character(signPosGene_i_MultipleAlnPos_j$Multiple_Aln_pos[1])
					binomPval <- pbinom(
						q=nForegroundMutGene_i-1
						, size=nForegroundMutGene_i+nBackgroundMutGene_i
						, prob=pBinom[Multiple_Aln_pos]
						, lower.tail=FALSE
						)
					})
			data.frame(
				Gene_Symbol=i
				, Multiple_Aln_pos=as.numeric(names(positionPval))
				, Pfam_ID=pfam
				, binomialPvalue=positionPval
				)
		})
		pfamOut <- do.call('rbind', pfamOut)
		out <- merge(pfamOut, signPos)
	} , BPPARAM=BPPARAM)
	pfamAnalysis <- do.call('rbind', pfamAnalysis)
	tryError <- tryCatch(rownames(pfamAnalysis) <- 1:nrow(pfamAnalysis) , error=function(e) "No Significant Mutations")
	if(!identical(tryError , "No Significant Mutations"))
		pfamAnalysis <- pfamAnalysis[ order(pfamAnalysis$metric) , ]
	else 
		pfamAnalysis <- tryError
	final_output <- list(AlignedSequence=pfamAnalysis , SingleSequence=singleSequence)
	return(final_output)
}
