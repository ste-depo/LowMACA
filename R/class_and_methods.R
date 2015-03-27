# class definition
setClass('LowMACA', 
	slots=c(
		arguments='list'
		, alignment='list'
		, mutations='list'
		, entropy='list'
		)
	, prototype=list(
		arguments=list(
			genes=NULL
			, pfam=NULL
			, input=''
			, mode=''
			, params=list(
				mutation_type=c('missense', 'all', 'truncating', 'silent')[1]
				, tumor_type='all'
				, min_mutation_number=1L
				, density_bw=0
				, clustal_cmd='clustalo'
				)
			, parallelize=list(
				getMutations=FALSE
				, makeAlignment=TRUE
				)
			)
		)
	,validity=function(object) {
	    object@arguments$params
	    if( class(object@arguments$params) != "list" )
	        return("Parameter value should be a list")
	    if( length(object@arguments$params) != 5 )
	        return("Parameter value should be a list of length 5")
	    if( !identical(
	        names(object@arguments$params),
	        c("mutation_type", "tumor_type", "min_mutation_number", 
	            "density_bw", "clustal_cmd")
	        ))
	        return("Parameter value should be a list whose names are:
	            mutation_type, tumor_type, min_mutation_number, density_bw and clustal_cmd")
	    if( ! object@arguments$params$mutation_type %in% 
	        c('missense', 'all', 'truncating', 'silent'))
	        return("mutation_type should be one of the following:
	            missense, all, truncating, silent")
	    if( class(object@arguments$params$tumor_type) != 'character' )
	        return('tumor_type must be a character')
	    if( !is.numeric(object@arguments$params$min_mutation_number) )
	        return('min_mutation_number must be an integer')
	    if( !is.numeric(object@arguments$params$density_bw) 
	    	& object@arguments$params$density_bw != 'auto' )
	        return('density_bw must be a number or "auto"')
	    if( class(object@arguments$params$clustal_cmd) != 'character' )
	        return('clustal_cmd must be a character')
    	TRUE
		}
	)

# getters
setGeneric('lmParams', function(object) standardGeneric('lmParams'))
setMethod('lmParams', 'LowMACA', function(object) {
	return(object@arguments$params)
	})

setGeneric('lmAlignment', function(object) standardGeneric('lmAlignment'))
setMethod('lmAlignment', 'LowMACA', function(object) {
	return(object@alignment)
	})

setGeneric('lmMutations', function(object) standardGeneric('lmMutations'))
setMethod('lmMutations', 'LowMACA', function(object) {
	return(object@mutations)
	})

setGeneric('lmEntropy', function(object) standardGeneric('lmEntropy'))
setMethod('lmEntropy', 'LowMACA', function(object) {
	return(object@entropy)
	})


setGeneric('parallelize', function(object) standardGeneric('parallelize'))
setMethod('parallelize', 'LowMACA', function(object) {
	return(object@arguments$parallelize)
	})

# setters
setGeneric('lmParams<-', function(object, value) standardGeneric('lmParams<-'))
setReplaceMethod('lmParams', 'LowMACA', function(object, value) {
	object@arguments$params <- value
	validObject(object)
	    #Check for correct clustalo command
    ClustalCommand <- Sys.which(object@arguments$params$clustal_cmd)
    if(ClustalCommand=="") {
        warning("The path to clustal omega is not correct. Change it ore use the web service. See ?setup for details")
    } else {
    	ClustalVersion <- system(paste(ClustalCommand ,  "--version") , intern=TRUE)
    	if(!grepl("^1.2" , ClustalVersion))
        	warning("Clustal Omega version could be not compatible. LowMACA was tested on 1.2.x")
	}
	object
})

setGeneric('parallelize<-', function(object, value) standardGeneric('parallelize<-'))
setReplaceMethod('parallelize', 'LowMACA', function(object, value) {
	object@arguments$parallelize <- value
	object
	})

# constructor
newLowMACA <- function(genes=NULL, pfam=NULL)
{
	# check arguments
	if( is.null(genes) & is.null(pfam) )
		stop('Either a gene or a pfam should be specified.')
	if( length(pfam)>1 )
		stop('Only one Pfam ID can be evaluated')
	# initialize the LowMaca object
	object <- new('LowMACA')
	if( !is.null(pfam) ) {
		object@arguments$mode <- 'pfam'
		object@arguments$pfam <- pfam
	} else object@arguments$mode <- 'genes'
	mode <- object@arguments$mode
	if( mode == 'genes' )
	{
	    # load annotation files
	    myAlias <- getMyAlias()
	    myUni <- getMyUni()
	    #
	    selectedColumns <- c('Gene_Symbol', 'Entrez', 'UNIPROT', 'AMINO_SEQ')
		geneID <- .checkGene_to_geneID(genes, myUni, myAlias)$Gene_Symbol
		selectedRows <- myUni$Gene_Symbol %in% geneID
		genesData <- myUni[selectedRows, selectedColumns]
		genesData$Pfam_ID <- '-'
		genesData$Envelope_Start <- 1
		genesData$Envelope_End <- nchar(genesData$AMINO_SEQ)
        seq_names <- paste(genesData[, 'Gene_Symbol']
            , genesData[, 'Pfam_ID']
            , genesData[, 'Entrez'] 
            , paste(genesData[, 'Envelope_Start'], 
                genesData[, 'Envelope_End'], sep=";") 
            , sep="|")
		rownames(genesData) <- seq_names
	} else {
		# load annotation files
	    myPfam <- getMyPfam()
	    # select rows and colums to match the input
	    selectedColumns <- c('Gene_Symbol', 'Pfam_ID', 'Entrez'
	    	, 'Envelope_Start', 'Envelope_End', 'UNIPROT', 'Pfam_Fasta')
	    selectedRows <- myPfam$Pfam_ID %in% pfam
	    if( !any(selectedRows) ) stop('Pfam name is not correct or is not mapped by LowMACA')
	    if( !is.null(genes) ) {
	    	# load annotation files
		    myAlias <- getMyAlias()
		    myUni <- getMyUni()
		    geneID <- .checkGene_to_geneID(genes, myUni, myAlias)$Gene_Symbol
		    selectedRows <- selectedRows & (myPfam$Gene_Symbol %in% geneID)
	    }
	    genesData <- myPfam[selectedRows, selectedColumns]
        seq_names <- paste(genesData[, 'Gene_Symbol']
            , genesData[, 'Pfam_ID']
            , genesData[, 'Entrez'] 
            , paste(genesData[, 'Envelope_Start'], 
                genesData[, 'Envelope_End'], sep=";") 
            , sep="|")
        colnames(genesData)[colnames(genesData) == 'Pfam_Fasta'] <- 'AMINO_SEQ'
        rownames(genesData) <- seq_names
	}
	genesData[, 'Gene_Symbol'] <- factor(genesData[, 'Gene_Symbol'])
	genesData[, 'Pfam_ID'] <- factor(genesData[, 'Pfam_ID'])
	genesData[, 'Entrez'] <- factor(genesData[, 'Entrez'])
	genesData[, 'UNIPROT'] <- factor(genesData[, 'UNIPROT'])
	## convert non-canonical amnio acids into their
	## most similar and canonical one
	##		U (selenocysteine) -> A (alanine)
	genesData$AMINO_SEQ <- gsub('U','A', genesData$AMINO_SEQ)
	## return the updated object
	object@arguments$input <- genesData
	return(object)

}

# methods
setGeneric('setup', function(object, repos=NULL, clustalo_filename=NULL 
	, mail=NULL , perlCommand="perl") standardGeneric('setup'))
setMethod('setup', 'LowMACA', function(object, repos=NULL, clustalo_filename=NULL 
	, mail=NULL , perlCommand="perl") {
	object <- alignSequences(object, clustalo_filename , mail , perlCommand)
	object <- getMutations(object, repos=repos)
	object <- mapMutations(object)
	return(object)
	})

setGeneric('alignSequences', function(object, clustalo_filename=NULL 
	, mail=NULL , perlCommand="perl") standardGeneric('alignSequences'))
setMethod('alignSequences', 'LowMACA', function(object, clustalo_filename=NULL 
	, mail=NULL , perlCommand="perl") {
	message("Aligning sequences...")
	clustal_cmd <- object@arguments$params$clustal_cmd
	genesData <- object@arguments$input
	if(!is.null(mail))
		.PerlModuleChecks(stop=TRUE , perl=perlCommand)
	object@alignment <- .clustalOAlign(genesData, clustal_cmd, clustalo_filename , mail , perlCommand)
	if( nrow(genesData)>1 ) {
		m <- object@alignment$CLUSTAL
		cm <- consensusMatrix(m)[-1,] # -1 removes gaps
		object@alignment$df <- data.frame(
			consensus=apply(cm,2, function(x) rownames(cm)[which.max(x)])
			, conservation=.Trident_Score(object@alignment$CLUSTAL)
			)
	} else {
		object@alignment$df <- data.frame(
			consensus=strsplit(genesData$AMINO_SEQ,'')[[1]]
			, conservation=rep(1, length(genesData$AMINO_SEQ))
			)				
	}
	return(object)
	})

setGeneric('getMutations', function(object, repos=NULL) standardGeneric('getMutations'))
setMethod('getMutations', 'LowMACA', function(object, repos=NULL) {
	genes <- object@arguments$input
	mutation_type <- object@arguments$params$mutation_type
	tumor_type <- object@arguments$params$tumor_type
	parallelGetMut <- object@arguments$parallelize$getMutations
	# outputFolder <- object@arguments$paths$output_folder
	if( is.null(repos) ) {
		message("Getting mutations from cancers studies...")
		gmOut <- .getGeneMutations(genes, 
			mutation_type=mutation_type, tumor_type=tumor_type, 
			parallelize=parallelGetMut)
	} else {
		message("Filtering mutations from local repository...")
		gmOut <- .getLocalGeneMutations(genes, 
			mutation_type=mutation_type, tumor_type=tumor_type, 
			localData=repos)
	}
	object@mutations$data <- gmOut$Mutations
	object@mutations$freq <- gmOut$AbsFreq
	return(object)
	})

setGeneric('mapMutations', function(object) standardGeneric('mapMutations'))
setMethod('mapMutations', 'LowMACA', function(object) {
	mut <- object@mutations$data
	if( nrow(mut)==0 )
		return(object)
	alignment <- object@alignment$ALIGNMENT
	alignmentLength <- max(alignment$Align)
	# # elaborate the mutations and map them with the alignment
	mut_aligned <- merge( mut , alignment[!is.na(alignment$Ref) , ]
					, by.x=c("Entrez" , "Gene_Symbol" , "Amino_Acid_Position") 
					, by.y=c("Entrez" , "Gene_Symbol" , "Ref") 
					, all.x=TRUE)
	mut_aligned <- mut_aligned[!is.na(mut_aligned$Align), ]
	# parameters for the analysis [can't be set by the user, so far]
	splitCriterion <- c('Gene_Symbol', 'Tumor_Type', 'domainID')
	mode <- object@arguments$mode
	if( mode == 'genes' ) {
		splitCriterion <- splitCriterion[1]
	} else splitCriterion <- splitCriterion[3]
	
	# split data
	mut_aligned.split <- split(mut_aligned, mut_aligned[[splitCriterion]])
	mut_aligned.extended <- matrix(0, nrow=length(mut_aligned.split)
		, ncol=alignmentLength)
	for(i in 1:length(mut_aligned.split)) {
		tmp <- sapply(
			split(mut_aligned.split[[i]]$Align, mut_aligned.split[[i]]$Align)
			, length
			)
		if( length(tmp)>0 ) 
			mut_aligned.extended[i,as.numeric(names(tmp))] <- tmp
	}
	rownames(mut_aligned.extended) <- names(mut_aligned.split)
	# filter data based on number of mutations
	minNMut <- object@arguments$params$min_mutation_number
	tokeep <- rowSums(mut_aligned.extended) >= minNMut
	if(any(!tokeep))
	{
		warning(
			paste("We excluded these genes (or domains) because they have less than"
				, minNMut, "mutations")
			, immediate.=TRUE)
		excludedGenes <- rownames(mut_aligned.extended[!tokeep, ])
		print(excludedGenes)
		object@mutations$excluded <- excludedGenes
	}
	mut_aligned.extended <- mut_aligned.extended[tokeep, , drop=FALSE]
	object@mutations$aligned <- mut_aligned.extended

	## add to the alignment data frame the frquency between the 
	## number of mutations that are feasible to come from a misalignment
	## of reads coming from another domain taken into consideration
	freq <- sapply(1:ncol(mut_aligned.extended), function(pos) {

		## alignment data
		selectedRows <- object@alignment$ALIGNMENT$Align == pos
		alnData <- object@alignment$ALIGNMENT[selectedRows, ]
		selectedColumns <- c('Gene_Symbol','Ref','Align')
		alnData <- alnData[!is.na(alnData$Ref), selectedColumns]
		colnames(alnData) <- c('Gene_Symbol', 'Amino_Acid_Position'
			, 'Consensus_Position')
		selectedColumns <- c('Gene_Symbol', 'Amino_Acid_Position'
			, 'Amino_Acid_Change')
		mutData <- object@mutations$data[, selectedColumns]
		mutData$Amino_Acid_Change <- substr(
			mutData$Amino_Acid_Change
			, nchar(mutData$Amino_Acid_Change)
			, nchar(mutData$Amino_Acid_Change)
			)
		colnames(mutData) <- c('Gene_Symbol', 'Amino_Acid_Position'
			, 'Mutation_letter')
		mergedData <- merge(mutData, alnData)
		## amino acids from mutated proteins
		mutAA <- table(mergedData$Mutation_letter)
		## amino acids in the alignment
		consensusAA <- object@alignment$df$consensus[pos]
		## frequency
		f <- sum(mutAA[names(mutAA)%in%consensusAA])/sum(mutAA)
		return(f)
		})
	object@alignment$df$misalnFreq <- freq
	# update the object
	return(object)

	})


setGeneric('entropy', function(object, bw=NULL) standardGeneric('entropy'))
setMethod('entropy', 'LowMACA', function(object, bw=NULL) {
	if( nrow(object@mutations$data)==0 ) {
		object@entropy$absval <- NA
		object@entropy$log10pval <- NA
		object@entropy$pvalue <- NA
		return(object)
	}
	mut_extended <- object@mutations$aligned
	# outputFolder <- object@arguments$paths$output_folder
	alignment    <- object@alignment
	# Uniform variable
	message("Making uniform model...")
	if( is.null(bw) ) bw <- object@arguments$params$density_bw else 
		object@arguments$params$density_bw <- bw
	if( bw=='auto' ) 
	{
		# calculate the bw from the global profile
		bw <- .profileDensity(colSums(mut_extended))$bw
	}
	message(paste('Assigned bandwidth:', round(bw,2)))
	object@entropy$bw <- bw
	weights <- .alnWeights(alignment)
	# if( plotUniform ) pdf(file.path(outputFolder , "Uniform_Model.pdf"))
	object@entropy$uniform <- .makeUniformModel(mut_extended, bw=bw, nboot=1000, 
		weights=weights, plotOUT=FALSE)
	# if( plotUniform ) dev.off()

	# Calculate the entropy values
	globalProfile <- colSums(mut_extended)
	uniform <- object@entropy$uniform
	absval <- .profileEntropy(globalProfile, norm=FALSE, bw=bw)
	log10pval <- .profileEntropy(globalProfile, norm=TRUE, bw=bw, model=uniform)
	pvalue <- 10^log10pval

	object@entropy$absval <- absval
	object@entropy$log10pval <- log10pval
	object@entropy$pvalue <- pvalue

	# null profile
	# Null profile calculated on the global profile
	nullOut <- .makeNullProfile(mut_extended, bw=bw, 
		nboot=1000, weights=.alnWeights(alignment))

	pvals <- nullOut$pvalue
	filteredPvals <- pvals
	filteredPvals[object@alignment$df$conservation < .1] <- NA
	qvals <- p.adjust(filteredPvals, method='BH')

	object@alignment$df <- cbind(
		object@alignment$df[, c('consensus', 'conservation')]
		, nullOut, qvalue=qvals
		, misalnFreq=object@alignment$df[, 'misalnFreq']
		)
	return(object)

	})

setMethod('show', 'LowMACA', function(object) {

	nItems <- nrow(object@arguments$input)
	message(paste('\nLowMACA object with', nItems, 'items:'))

	items <- object@arguments$input$Gene_Symbol
	itemsText <- paste(head(items), collapse=', ')
	message(paste(itemsText, ifelse(nItems>6, ', ...\n', '\n'), sep=''))

	if( length(object@entropy)>0 ) {
		objectPvalue <- object@entropy$pvalue
		objectPvalue <- signif(objectPvalue, 5)
		message(paste('The p-value of the global aligned profile is:', objectPvalue, '\n'))
	}
	})

###################
###### SUMMARY METHODS
############################

setGeneric('lfm', function(object, metric='qvalue', threshold=.05, conservation=0.1) 
	standardGeneric('lfm'))
setMethod('lfm', 'LowMACA', function(object, metric='qvalue', threshold=.05, conservation=0.1) {
	if( nrow(object@mutations$data)==0 ) {
		message('No mutations available for this object.')
	} else {
		if( metric == 'qvalue' ) {
			if( conservation != 0.1 ) {
				## in case a new conservation threshold is given,
				## recalculate the qvalue
				filteredPvals <- object@alignment$df$pvalue
				filteredPvals[object@alignment$df$conservation < conservation] <- NA
				qvalue <- p.adjust(filteredPvals, method='BH')
			} else {
				## use qvalue already stored
				qvalue <- object@alignment$df$qvalue
			}
			signifPos <- which(qvalue < threshold)
		} else {
			## use pvalue
			pvalue <- object@alignment$df$pvalue
			signifPos <- which(pvalue < threshold)
		}
		out <- lapply(signifPos, function(pos) {
			selectedRows <- object@alignment$ALIGNMENT$Align == pos
			alnData <- object@alignment$ALIGNMENT[selectedRows, ]
			selectedColumns <- c('Gene_Symbol', 'Ref', 'Envelope_Start'
				, 'Envelope_End')
			alnData <- alnData[!is.na(alnData$Ref), selectedColumns]
			colnames(alnData) <- c('Gene_Symbol', 'Amino_Acid_Position'
				, 'Envelope_Start', 'Envelope_End')
			selectedColumns <- c('Gene_Symbol', 'Amino_Acid_Position'
				, 'Amino_Acid_Change', 'Sample', 'Tumor_Type')
			mutData <- object@mutations$data[, selectedColumns]
			lfm <- merge(mutData, alnData)
			if( nrow(lfm)>0 )
				lfm$Multiple_Aln_pos <- pos
				lfm$metric <- switch(metric, 'pvalue'=pvalue[pos], 'qvalue'=qvalue[pos])
			return(lfm)
			})
		out <- do.call('rbind', out)
		# load data
	    myUni <- getMyUni()
	    selectedColumns <- c('Gene_Symbol', 'Entrez', 'Entry', 'UNIPROT'
	    	, 'Chromosome', 'Protein.name')
		myUniSmall <- myUni[, selectedColumns]
		out <- merge(out, myUniSmall)
		return(out)
	}
	})

############
###### GLOBAL ANALYSIS
#############################



############
###### PLOTTING METHODS
###############################

setGeneric('nullProfile', function(object) standardGeneric('nullProfile'))
setMethod('nullProfile', 'LowMACA', function(object) {

	if( nrow(object@mutations$data)==0 ) {
		message('No mutations available for this object.')
	} else {

		mean <- object@alignment$df$mean
		lowerThreshold <- object@alignment$df$lTsh
		upperThreshold <- object@alignment$df$uTsh
		profile <- object@alignment$df$profile
		pvalue <- object@alignment$df$pvalue
		qvalue <- object@alignment$df$qvalue
		misalnFreq <- object@alignment$df$misalnFreq

		#
	    qvalSignif_x <- which(qvalue < 5e-2)
	    qvalSignif_y <- profile[qvalSignif_x] + max(profile)/20
	    over <- profile > upperThreshold
	    ylim <- range(c(profile, upperThreshold, lowerThreshold, qvalSignif_y))

	    # red bars of the resudues over the threshold
	    plot(profile, main='', type='h'
	        , ylim=ylim, bty='n',xaxt='n',xlab='', col='orange', lwd=2)
	    # black bars of the other resudues and of residues over, but only
	    # the part below the threshold
	    blackProfile <- sapply(1:length(profile), 
	    	function(i) min(profile[i], upperThreshold[i]))
	    lines(blackProfile, col='black', lwd=2, type='h')
	    lines(upperThreshold, lty=2, lwd=2, col='blue')
	    ### plot an asterisk above the residues which are significant in 
	    # terms of the pvalue
	    if( length(qvalSignif_x) > 0 ) {
	    	text(qvalSignif_x, qvalSignif_y, '*', cex=2, col='red')
	    	## if more than 70% of mutations can be mapped on
	    	## another protein of the domain, plot the BLUE ASTERISK
	    	## instead of RED
	    	misaligned <- misalnFreq[qvalSignif_x] > .7
	    	if(any(misaligned)) {
	    		text(qvalSignif_x[misaligned], qvalSignif_y[misaligned]
	    			, '*', cex=2, col='blue')
	    		legend('topright', pch='*', col='blue', legend='possibly misaligned', cex=2)
	    	}
	    }
	}

	})

setGeneric('lmPlot', function(object) standardGeneric('lmPlot'))
setMethod('lmPlot', 'LowMACA', function(object) {
	if( nrow(object@mutations$data)==0 ) {
		message('No mutations available for this object.')
	} else {

		mode <- object@arguments$mode
		nObj <- nrow(object@arguments$input)

		origMAlign <- object@alignment$CLUSTAL
		m <- consensusMatrix(origMAlign)
		motif <- pcm2pfm(m)
		motif <- new('pfm', mat=motif, name='', color=colorset(alphabet='AA'))
		log10pval <- object@entropy$log10pval
		pvalue <- object@entropy$pvalue

		par(mar=c(2,4,2,2))
		if( nObj > 1 ) {
			layoutMatrix <- as.matrix(c(1,1,2,2,3,4,4))
			plotType <- 4
		} else if( mode == 'genes' ) {
			layoutMatrix <- as.matrix(c(1,1,1,1,3,2,2,2,2))
			plotType <- 3
		} else {
			layoutMatrix <- as.matrix(c(1,1,2,2))
			plotType <- 2
		}
		layout(layoutMatrix)

		## plot 1
		if( plotType == 3) par(mar=c(0,4,2,2))
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
		lenAln <- ncol(object@mutations$aligned)
		mut_aligned <- object@mutations$aligned
		# if( plotType == 3 ) colnames(mut_aligned) <- 1:length(mut_aligned)
		barplot(
			mut_aligned
			, col=myPalette(nrow(object@mutations$aligned))
			, border=if(lenAln<300) 'black' else NA
			, main=paste( "Mutations by position\nLog10 P-Value:" , round(log10pval , 2) , 
					"P-Value:" , signif( pvalue , 2 ) , "Bw:" , signif(object@entropy$bw,3) ) 
			, ylab='Mutations'
			)
		## plot 2
		if( plotType == 3) par(mar=c(2,4,0,2))
		nullProfile(object)
		if( plotType == 3) {
			pSeq <- round(seq(1, nchar(object@arguments$input$AMINO_SEQ), length.out=20))
			axis(1, at=pSeq, labels=pSeq)
		}
		## if there is more than one object also plot 
		## conservation plots	
		if( nObj > 1 ) {
			## plot 3
			barplot(object@alignment$df$conservation
				, col='darkgoldenrod1', ylab='Conservation')
			## plot 4
			plot(motif, ylab='Logo')
		}

		#############
		## in case the analysis os on a single gene 
		## plot its domains
		############### 
		if( plotType == 3 ) {

			myPfam <- getMyPfam()
			gene <- as.character(object@arguments$input$Gene_Symbol)
			domains <- myPfam[myPfam$Gene_Symbol==gene
				, c("Envelope_Start" , "Envelope_End" , "Pfam_Name")]
			## create empty plot
			plot.new()
			plot.window(
				xlim=c(1,nchar(object@arguments$input$AMINO_SEQ))
				, ylim=c(0,0.05)
				)
			par(mar=c(0,0,0,0))
			## plot domains
			if(nrow(domains)>0) {
				for (i in 1:nrow(domains)) {
					xleft=domains[i , "Envelope_Start"]
					xright=domains[i , "Envelope_End"]
					ytop=0.05
					ybottom=0
					col=topo.colors(nrow(domains) , alpha=0.5)[i]
					characters <- nchar(domains[i , "Pfam_Name"])
					rect(xleft=xleft , xright=xright 
						, ytop=ytop , ybottom=ybottom 
						, col=col )
					if(characters<=3){
						text(x=(xright+xleft)/2 , y=0.025 
							, domains[i , "Pfam_Name"] 
							, font=2)
					} else {
						if((xright-xleft)<=100){
							text(x=(xright+xleft)/2 , y=0.025 
								, domains[i , "Pfam_Name"] 
								, font=2 , cex=0.8)
						} else {
							text(x=(xright+xleft)/2 , y=0.025 
								, domains[i , "Pfam_Name"] 
								, font=2)
						}
					}
				}
			} else {
				text(x=nchar(object@arguments$input$AMINO_SEQ)/2, y=0.025
					, 'no pfam domains within gene sequence', cex=1.5)
			}
		}

	}

	})


setGeneric('bpAll', function(object) standardGeneric('bpAll'))
setMethod('bpAll', 'LowMACA', function(object) {
	if( nrow(object@mutations$data)==0 ) {
		message('No mutations available for this object.')
	} else {

		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

		mut_aligned_extended <- object@mutations$aligned
		gene_subset <- object@mutations$aligned
		colnames(gene_subset) <- as.character(1:ncol(gene_subset))
		log10pval <- ifelse( is.numeric(object@entropy$log10pval) , object@entropy$log10pval , 0)
		pvalue <- ifelse( is.numeric(object@entropy$pvalue) , object@entropy$pvalue , 1)
		barplot(gene_subset 
			, col=myPalette(nrow(gene_subset))
			, border=NA
			, las=2
			, cex.names=0.3 
			, main=paste( "Mutations by position\nEntropy Z-Score:" , round(log10pval , 2) , 
					"P-Value:" , signif( pvalue , 2 ) ) 
			)
		legend("topright" 
			, legend=rownames(gene_subset) 
			, fill=myPalette(nrow(gene_subset)) 
			, col=myPalette(nrow(gene_subset)) 
			, border=myPalette(nrow(gene_subset)) 
			, cex=0.5)
	}
	})


setGeneric('protter', function(object, filename='protter.png', threshold=5e-2) standardGeneric('protter'))
setMethod('protter', 'LowMACA', function(object, filename='protter.png', threshold=5e-2) {
	pvalues <- object@alignment$df$pvalue
	qvalues <- object@alignment$df$qvalue
	message(paste('Writing', filename, '...'))
	Sequence <- paste(object@alignment$df$consensus, collapse='')
	Mutation_pvalues <- paste(which(pvalues < threshold), collapse=',')
	Mutation_qvalues <- paste(which(qvalues < threshold), collapse=',')
	WebQuery <- paste("http://wlab.ethz.ch/protter/create?seq=", Sequence, 
		"&tm=auto&mc=lightsalmon&lc=blue&tml=numcount&numbers&legend&n:signal%20peptide,cc:white,fc:blue,bc:blue=Phobius.SP&n:N-glyco%20motif,s:box,fc:forestgreen,bc:forestgreen=(N).[ST]&n:MUT_pvalue,bc:orange="
		, Mutation_pvalues, "&n:MUT_qvalue,bc:red=", Mutation_qvalues, "&format=png",sep="")
	if(Sys.info()['sysname']=="Windows")
		download.file(WebQuery, destfile=filename , mode="wb")
	else
		download.file(WebQuery, destfile=filename)
	message(paste("Protter plot saved as:" , normalizePath(filename)))
	})
