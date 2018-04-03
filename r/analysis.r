
# nourdine.bah@crick.ac.uk
# philip.east@crick.ac

library(R.utils)
library(optparse)
library(rtracklayer)
library(plyr)
library(dplyr)
library(reshape)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(DOSE)
library(clusterProfiler)
library(biomaRt)
library(openxlsx)
library(VennDiagram)
library(grid)
library(gridExtra)
library(grDevices)
library(RColorBrewer)
library(mixOmics)
library(vsn)


###############################################################################
##                                                                           ##
##                      SESSION AND CONTACT INFORMATION                      ##
##                                                                           ##
###############################################################################

EMAIL <- paste0(Sys.getenv("USER"), "@crick.ac.uk")
R_VERSION <- base::version[['version.string']]
ENSEMBL_VERSION <- 86

DESEQ2_VERSION <- installed.packages()[grep("DESeq2",
														  rownames(installed.packages())),
															"Version"]
LFC_INTERPRETATION <-
	"a positive LFC for A.vs.B means a larger expression in A"

information.names <- list(email="Contact:",
								  r.version="R version:",
								  ensembl.version="Ensembl version:",
								  deseq2.version="DESeq2 version:",
								  lfc.interpretation="LFC interpretation:")

information <- list(email=EMAIL,
						  r.version=R_VERSION,
						  ensembl.version=ENSEMBL_VERSION,
						  deseq2.version=DESEQ2_VERSION,
						  lfc.interpretation=LFC_INTERPRETATION)

information.df <- data.frame(names=unlist(information.names),
									  info=unlist(information))



###############################################################################
##                                                                           ##
##                            GGPLOT2 THEME                                  ##
##                                                                           ##
###############################################################################

common_theme <- theme(
							 plot.title=element_text(hjust=0.5)
							 )



###############################################################################
##                                                                           ##
##                               FUNCTIONS                                   ##
##                                                                           ##
###############################################################################

###
get_counts <-
	# take a file path and return a list of two data frames : the whole
	# result file ($df) and just the count ($count)
	function(filepath,
				count_colname="expected_count",
				geneid_colname="gene_id") {

		# extract the sample ID from the file name
		sample_id <- sub("\\..*$", "", basename(filepath))

		# load the result file and add the sample ID
		df <- read.table(filepath, header=T)
		df$sample_id <- sample_id

		# the count data frame is composed of just one column storing the counts.
		# This column is named with the sample ID and the row names are
		# corresponding to the gene IDs.
		count <- data.frame(df[,count_colname])
		rownames(count) <- df[,geneid_colname]
		names(count) <- sample_id

		return(list(df=df, count=count))
	}


###
check_path <-
	# assure paths existence
	function(args) {

		# need this for absolute path checking
		if ( ! "package::R.utils" %in% search() ) { library(R.utils) }


		path_args <- names(args)[grep("-(dir|file)$", names(args))]
		lapply(path_args,

				 function(arg) {

					 # the argument path
					 path <- args[[arg]]

					 if (grepl("-dir$", arg)) {

						 # each path as to be an absolute path
						 if (!isAbsolutePath(path)) {
							 #stop(paste("Error: the", arg,
							 #				"is not an absolute path.\n")
						 }
						
						 # the output need to created if doesn't exist
						 if (path == "output-dir") {

							 # creates directory if not exists
							 if (!dir.exists(path)) {
								 dir.create(path, recursive=T)
							 }

						 } else {

							 # creates directory if not exists
							 if (!dir.exists(path)) {
								 stop(paste("Error: the", arg,
												"argument directory does not exist.\n"))
							 }
						 }

					 } else if (grepl("-file$", arg)) {
						 if ( !isAbsolutePath(path) || !file.exists(path) ) {
							 #stop(paste("Error: the", arg, "is not valid.\n"))
						 }
					 }

				 })
	}


###
ma_df <-
	# create a data frame containing the M and the A values
	function(comb, df, label_colname, count_colname) {

		# progressing info
		message <- paste("MA plot for", paste(comb, collapse=" VS "))
		cat(paste0(message, "\n"))

		# each sample
		x <- log2( df[ df[,label_colname]==comb[1], count_colname ] + 1 )
		y <- log2( df[ df[,label_colname]==comb[2], count_colname ] + 1 )

		# M and A values
		m <- x - y
		a <- (x + y) / 2

		# the data frame and title
		d <- data.frame(m, a)
		d$c <- paste(comb[1], "\nVS\n", comb[2])

		return(d)
	}


###
maplot <-
	# produces the MA plots between replicates for each condition
	function(df,
				exp_condition,
				condition_colname="condition",
				label_colname="sample_id",
				count_colname="expected_count",
				sample_size=10000,
				pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# need this for sample_n
		if ( ! "package::dplyr" %in% search() ) { library(dplyr) }

		# need to plot
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		# need plyr's ddply function
		if ( ! "package::plyr" %in% search() ) { library(plyr) }

		## xxxxxxxx #############################################################

		# the current condition
		condition <- unique(df[,condition_colname])
		if (length(condition)!=1) {
			stop("Error: condition is not unique for MA plot creation.")
		}

		# the replicates for the current condition
		samples <- unique(df[,label_colname])
		if (length(samples)<2) {
			message <- paste("Warning: just less than two samples for", condition)
			cat(paste0(message, "\n"))
			return(1)
		}

		# all possible combinations between replicates
		combL <- combn(samples, 2, simplify=F)
		combL <- lapply(combL, sort)

		# built MA-values-containing data frame for each comparison and fuse them
		l <- lapply(combL, ma_df, df, label_colname, count_colname)
		madf <- do.call(rbind, l)

		# have to sample the data frame because plot printing takes too long
		# otherwise
		madf <- plyr::ddply(madf, "c", sample_n, sample_size, replace=T)

		# build the plot
		g <- ggplot(madf, aes_string(x="a", y="m"))
		g <- g + geom_point(size=1.5, alpha=1/5)
		g <- g + geom_hline(yintercept=0, color="blue3")
		g <- g + stat_smooth(se=FALSE, method="loess", color="red3")
		g <- g + ylab("M")
		g <- g + xlab("A")
		g <- g + facet_wrap(~ c, ncol=5)
		g <- g + ggtitle(paste0("MA plots (", condition, ")"))
		g <- g + common_theme + theme(aspect.ratio=1,
												strip.text=element_text(size=14))

		# never remember it it's necessary on HPC
		print(g)

		# save the plot
		dir <- file.path(pdf_dir, "maplot")
		if (!dir.exists(dir)) { dir.create(dir, recursive=T) }
		filename <- paste0(sub(CONDITION_SEPARATOR, "_", condition), ".pdf")
		filename <- tolower(filename)
		path <- file.path(dir, filename)
		ggsave(path)

		return(0)
	}


###
normalisation <-
	# plots the different normalisation methods
	function(count, sample_size=10000, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# need edgeR normalisation methods
		if ( ! "package::edgeR" %in% search() ) { library(edgeR) }

		# need reshape's melt function
		if ( ! "package::reshape" %in% search() ) { library(reshape) }

		# need plyr's ddply function
		if ( ! "package::plyr" %in% search() ) { library(plyr) }

		# need this for sample_n
		if ( ! "package::dplyr" %in% search() ) { library(dplyr) }

		# need to plot
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		## xxxxxxxx #############################################################

		

		# build edgeR matrix
		dgel <- edgeR::DGEList(count)

		## the 4 different methods ##############################################

		# factors for total read count normalisation method
		count.nmlzd.cpm <- edgeR::cpm(dgel)
		tot <- data.frame(1e6 / colSums(count))
		names(tot) <- c("tot")

		# factors for upperquartile normalisation method
		upq <- edgeR::calcNormFactors(dgel, method="upperquartile",
												p=.75)$samples$norm.factors
		names(upq) <- c("upq")

		# factors for relative log expression normalisation method
		rle <- edgeR::calcNormFactors(dgel, method="RLE")$samples$norm.factors
		names(rle) <- c("rle")
		
		# factors for trimmed mean of m-values normalisation method
		tmm <- edgeR::calcNormFactors(dgel, method="TMM")$samples$norm.factors
		names(tmm) <- c("tmm")

		## xxx x xxxxxxxxx xxxxxxx ##############################################


		# fuse all the 4 methods
		norm <- as.data.frame(t(cbind(tot, upq, rle, tmm)))
		norm <- rbind(norm, non=rep(1, 6))

		# apply the normalisation methods to the data
		l <- lapply(row.names(norm),
						function(x) {
							c <- count
							r <- norm[x,]
							for (n in names(r)) {
								c[,n] <- as.numeric(c[,n]) * as.numeric(r[n])
							}
							names(c) <- gsub("\\..*", "", names(c))
							c$method <- x
							return(c)
						})
		df <- do.call(rbind, l)

		# change the name so it can be properly display by ggplot2
		df$method <- factor(df$method,
								  levels=c("non", "tot", "tmm", "upq", "rle"))
		df$method <- gsub("non", "Non normalised", df$method)
		df$method <- gsub("tot", "Total read count", df$method)
		df$method <- gsub("tmm", "Trimmed M's mean", df$method)
		df$method <- gsub("upq", "Upper quantile", df$method)
		df$method <- gsub("rle", "Relative Log", df$method)

		# build the data frame for the rank-mean plot
		df <- plyr::ddply(df, "method",
						function(x) {

							# the for each gene
							rank <- cbind(x, mean=apply(x[,names(count)], 1, mean))

							# the rank of the mean
							rank$rank <- rank(as.numeric(rank$mean))

							# transform the data frame to the long format
							rank <- reshape::melt(rank,
														 id.vars=c("method", "mean", "rank"))
							names(rank) <-
								c("method", "mean", "rank", "sample", "count")
							rank$count <- as.numeric(rank$count)

							# because otherwise the plot takes ages to print
							rank <- dplyr::sample_n(rank, sample_size, replace=T)

							# remove 0s
							rank <- subset(rank, count > 0)

							return(rank)
						})
		df$sample <- factor(df$sample, levels=names(count))

		# local theme for the rank-mean plot
		theme <- theme(plot.title=element_text(size=22, hjust=.5),
							legend.text=element_text(size=15),
							legend.title=element_text(size=17),
							axis.title=element_text(size=17),
							strip.text=element_text(size=15),
							legend.position="top",
							aspect.ratio=1)

		# build the rank-mean plot
		g <- ggplot(df, aes_string(x="rank", y="log2(count+1)", color="sample"))
		g <- g + geom_point(size=.4)
		g <- g + facet_wrap(as.formula("~ method"), ncol=3, scales="free")
		g <- g + xlab("Rank of the mean")
		g <- g + ylab("Mean raw counts (log2(count+1))")
		g <- g + scale_color_discrete(name="Sample")
		g <- g + ggtitle("Normalisation methods")
		g <- g + common_theme + theme
		g <- g + guides(colour=guide_legend(nrow=4))
		print(g)

		# save the plot
		path <- file.path(pdf_dir, "/normalisation.pdf")
		ggsave(path)

		return(df)
	}


###
plot_pca <-
	# takes a DESeq2 data set object, plots the PCA from it and returns the
	# plot
	function(dds, design, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# need the plotPCA method from the DESeq2 package
		if ( ! "package::DESeq2" %in% search() ) { library(DESeq2) }

		# need to plot
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		## xxxxxxxx #############################################################

		
		# let DESeq do the PCA for us
		pca <- DESeq2::plotPCA(dds,
									  intgroup="condition",
									  ntop=1000,
									  returnData=T)

		# fuse the PCA data with the design data
		design <- design[,-grep("^condition$", names(design))]
		pca <- merge(pca, design, by.x="name", by.y="sample_id")

		# get the conditions of the experiment from the design matrix
		reserved_words <- c("sample", "sample_id", "condition")
		reserved_words <- paste0("^", reserved_words, "$")
		regex <- paste(reserved_words, collapse="|")
		conditions <- sort(names(design)[-grep(regex, names(design))])

		# create the aesthetics depending on the number of experimental
		# conditions
		if (length(conditions)==0) {

			stop("Error: Something very very weird happened :/")

		} else if (length(conditions)==1) {

			aes <- aes_string(
									x="PC1",
									y="PC2",
									label="name",
									color=conditions[1]
									)

		} else if (length(conditions)==2) {

			aes <- aes_string(
									x="PC1",
									y="PC2",
									label="name",
									color=conditions[1],
									shape=conditions[2]
									)

		} else {

			aes <- aes_string(
									x="PC1",
									y="PC2",
									label="name",
									color=conditions[1],
									shape=conditions[2],
									fill=conditions[3]
									)
		}

		# build the plot and save it
		g <- ggplot(pca, aes) + geom_point(size=3)
		g <- g + ggtitle(paste0("Principal component analysis"))
		g <- g + geom_label_repel(aes,
										  segment.colour="black",
										  fontface='bold',
										  box.padding=unit(0.35, "lines"),
										  point.padding=unit(0.5, "lines"))
		g <- g + scale_fill_brewer(palette="Set2")
		g <- g + scale_color_brewer(palette="Set1")
		g <- g + common_theme
		print(g)
		path <- file.path(pdf_dir, "pca_plot.pdf")
		ggsave(path)

		return(g)
	}


###
plot_heatmap <-
	# plots heatmap from a DESeq data set object
	function(dds, design, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# need if for colorRamPalette
		if ( ! "package::grDevices" %in% search() ) { library(grDevices) }

		# need if for brewer.pal
		if ( ! "package::RColorBrewer" %in% search() ) { library(RColorBrewer) }

		# need if for cim
		if ( ! "package::mixOmics" %in% search() ) { library(mixOmics) }

		## xxxxxxxx #############################################################

		# compute the distance
		dist <- as.matrix(dist(t(assay(dds))))
		dist <- dist/max(dist)

		# get the conditions of the experiment from the design matrix
		reserved_words <- c("sample", "sample_id", "condition")
		reserved_words <- paste0("^", reserved_words, "$")
		regex <- paste(reserved_words, collapse="|")
		conditions <- sort(names(design)[-grep(regex, names(design))])

		# three color sets for three conditions max
		palettes <- c("Set1", "Set2", "Set3")

		# map condition to color
		maps <-
			lapply(1:length(conditions),
					 function(n, d, p) {

						 c <- conditions[n]
						 factors <- sort(as.character(unique(d[,c])))
						 nfact <- length(factors)

						 # don't understand why the minimum is 3...
						 if (nfact<3) {
							 colors <- brewer.pal(3, p[n])
						 } else {
							 colors <- brewer.pal(nfact, p[n])
						 }

						 colors <- colors[1:nfact]
						 map <- data.frame(f=factors, c=colors)

						 return(map)

					 }, design, palettes)
		names(maps) <- conditions

		# add the color to the design data frame
		design_colors <-
			mapply(function(map, n) {
					 df <- merge(design, map, by.x=n, by.y="f")
					 return(df)
					 }, maps, names(maps), SIMPLIFY=F)

		# the color vector that will be added to the heatmap
		side_colors <-
			mapply(function(d, n) {
					 rows <- match(colnames(dist), d[,"sample_id"])
					 return(d[rows,"c"])
					 }, design_colors, names(design_colors))

		# generate the palette
		hmcol <- grDevices::colorRampPalette(brewer.pal(9, "GnBu"))(16)

		# the legend annotation for each color
		annot_maps <-
			mapply(function(d, n) {
					 d[,"f"] <- paste0(d[,"f"], " (", n, ")")
					 return(d)
					 }, maps, names(maps), SIMPLIFY=F)

		# all the colors and factors
		general_map <- do.call(rbind, annot_maps)

		# the plot legend
		legend <-
			list(
				  legend=general_map[,"f"],
				  # fucking cast :/ :/ :/ .....
				  col=as.character(general_map[,"c"]),
				  title="Condition",
				  cex=.7
				  )

		# build the plot and save it
		path <- file.path(pdf_dir, "distance.pdf")
		pdf(path)
		cim(
			 dist,
			 color=rev(hmcol),
			 col.sideColors=side_colors,
			 row.sideColors=side_colors,
			 symkey=F,
			 margins=c(9, 9),
			 legend=legend
			 )
		dev.off()
		return(maps)
	}


###
meansdplot <-
	# builds a mean sd plot, saves it and returns it
	function(dds, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# need vsn's meanSdPlot method
		if ( ! "package::gridExtra" %in% search() ) { library(gridExtra) }

		# need DESeq for the transformation methods
		if ( ! "package::DESeq2" %in% search() ) { library(DESeq2) }

		# need to plot
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		# need gridExtra's arrangeGrob method
		if ( ! "package::gridExtra" %in% search() ) { library(gridExtra) }

		## xxxxxxxx #############################################################

		# non zero counts
		notAllZero <- (rowSums(counts(dds))>0)

		# different types of transformation
		ntd <- DESeq2::normTransform(dds)
		rld <- DESeq2::rlog(dds, blind=F)
		vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=F)

		# theme
		t <- theme(
					  plot.title=element_text(hjust=.5),
					  aspect.ratio=1,
					  legend.position="top"
					  )

		# plots
		msdplot.raw <- vsn::meanSdPlot(as.matrix(count)[notAllZero,])$gg
		msdplot.ntd <- vsn::meanSdPlot(assay(ntd)[notAllZero,])$gg
		msdplot.vsd <- vsn::meanSdPlot(assay(vsd)[notAllZero,])$gg
		msdplot.rld <- vsn::meanSdPlot(assay(rld)[notAllZero,])$gg

		# add title and theme
		msdplot.raw <- msdplot.raw + t + ggtitle("Raw counts")
		msdplot.ntd <- msdplot.ntd + t + ggtitle("Log")
		msdplot.rld <- msdplot.rld + t + ggtitle("Regularised log")
		msdplot.vsd <- msdplot.vsd + t + ggtitle("Variance stabilizing")
		
		# diplay all the plots
		grob <- gridExtra::arrangeGrob(
												 msdplot.raw,
												 msdplot.ntd,
												 msdplot.vsd,
												 msdplot.rld,
												 nrow=2
												 )
		
		# saves the plot
		path <- file.path(pdf_dir, "meansdplot.pdf")
		ggsave(path, grob)

		return(grob)
	}


###
create_dds <-
	# creates a deseq object containing just the contrasts of interest
	function(contrast,
				count,
				design,
				condition_colname="condition",
				label_colname="sample_id") {

		# need DESeq2
		if ( ! "package::DESeq2" %in% search() ) { library(DESeq2) }

		# progressing info
		message <- paste("Create DESeq2 object for",
							  paste(contrast, collapse=" VS "))
		cat(paste0("\n", message, "\n"))

		# get the samples associated to the current contrasts
		sample <- design[design[,condition_colname] %in% contrast,
							  c(label_colname, condition_colname)]

		# get the conditions corresponding to the samples
		id <- sample[,label_colname]
		condition <- as.factor(sample[,condition_colname])

		# build the deseq object
		dds <- DESeq2::DESeqDataSetFromMatrix(
														  countData= round(count[,id]),
														  colData = DataFrame(condition),
														  design = as.formula("~ condition")
														  )
		dds <- DESeq2::DESeq(dds)

		return(dds)
	}


###
get_results_from_dds <-
	# returns results of a conditions comparison and saves the corresponding
	# MA plot
	function(dds, conditions, alpha=.05, coefficient=2, pdf_dir=PDF_DIRPATH) {

		# need DESeq2
		if ( ! "package::DESeq2" %in% search() ) { library(DESeq2) }

		# compute results
		contrast <- c("condition", as.character(conditions))
		results <- DESeq2::results(dds, alpha=alpha, contrast=contrast)

		# pdf path and file name
		dir <- file.path(pdf_dir, "deseq2_maplot")
		title <- paste(conditions, collapse=" VS ")
		name <- paste(conditions, collapse="_vs_")
		if (!dir.exists(dir)) {
			dir.create(dir, recursive=T)
		}

		# progressing info
		message <- paste("Get results for", title)
		cat(paste0("\n", message, "\n"))

		# save the MA plot
		filename <- paste0(tolower(name), ".pdf")
		pdf(file=file.path(dir, filename))
		DESeq2::plotMA(results, main=title, alpha=alpha)
		dev.off()

		# save the shrunk MA plot
		shrunk_results <- lfcShrink(dds, coef=coefficient, res=results)
		filename <- paste0("shrunk_", tolower(name), ".pdf")
		pdf(file=file.path(dir, filename))
		DESeq2::plotMA(shrunk_results,
							main=paste(title, "(shrunk)"),
							alpha=alpha)
		dev.off()

		return(results)
	}


###
get_significant_results <-
	# applies a threshold on the adjusted p value of a DESeq2 results object
	function(results, threshold) {
		results <- results[!is.na(results$padj),]
		results <- results[results$padj<threshold,]
		return(results)
	}


###
ontology_dotplot <-
	# changes the way DOSE plots the data in order to have shorter and unique
	# y-axis labels
	function(g) {

		# need DOSE's doplot function
		if ( ! "package::DOSE" %in% search() ) { library(DOSE) }

		# the dotplot function from DOSE
		dplot <- DOSE::dotplot(g)

		# this can create duplicates
		dplot$data$Description <- gsub(",.*", "", dplot$data$Description)

		# number duplicates
		dplot$data$Description <- make.unique(dplot$data$Description, sep = "_")

		# order depending on the gene ratio rank
		idx <- order(dplot$data$GeneRatio, decreasing=F)
		dplot$data$Description <-
			factor(dplot$data$Description, levels=dplot$data$Description[idx])

		# we want the plot itself
		return(dplot)
	}


###
ontology <-
	# performs the three types of gene ontology and the pathways analysis,
	# plots them and returns them as a list
	function(results,
				analysis,
				name,
				binomial,
				annotation,
				pvalue_threshold=.01,
				qvalue_threshold=.05,
				adj_method="BH",
				pdf_dir=PDF_DIRPATH) {

		# progressing info
		message <- paste("Perform gene ontology analysis on", name)
		cat(paste0("\n", message, "\n"))

		## packages #############################################################

		# need clusterProfiler for GO analysis
		if ( ! "package::clusterProfiler" %in% search() ) {
			library(clusterProfiler)
		}

		# need to plot
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		# need gridExtra's arrangeGrob method
		if ( ! "package::gridExtra" %in% search() ) { library(gridExtra) }

		## xxxxxxxx #############################################################


		## organism choosing ####################################################

		if (tolower(binomial)=="hs") {

			# load the map and assign it
			if ( ! "package::org.Hs.eg.db" %in% search() ) {
				library(org.Hs.eg.db)
			}
			db = org.Hs.eg.db

			# the annotation
			if (missing(annotation)) { annotation = "ENSEMBL" }


		} else if (tolower(binomial)=="mm") {

			# load the map and assign it
			if ( ! "package::org.Mm.eg.db" %in% search() ) {
				library(org.Mm.eg.db)
			}
			db = org.Mm.eg.db

			# the annotation
			if (missing(annotation)) { annotation = "ENSEMBL" }


		} else if (tolower(binomial)=="dm") {

			# load the map and assign it
			if ( ! "package::org.Dm.eg.db" %in% search() ) {
				library(org.Dm.eg.db)
			}
			db = org.Dm.eg.db

			# the annotation
			if (missing(annotation)) { annotation = "FLYBASE" }


		} else if (tolower(binomial)=="sc") {

			# load the map and assign it
			if ( ! "package::org.Sc.sgd.db" %in% search() ) {
				library(org.Sc.sgd.db)
			}
			db = org.Sc.sgd.db

			# the annotation
			if (missing(annotation)) { annotation = "ORF" }


		} else if (tolower(binomial)=="xl") {

			return(1)

			# load the map and assign it
			if ( ! "package::org.Xl.eg.db" %in% search() ) {
				library(org.Xl.eg.db)
			}
			db = org.Xl.eg.db

			# the annotation
			if (missing(annotation)) { annotation = "" }


		} else if (tolower(binomial)=="pf") {

			return(1)

			# load the map and assign it
			if ( ! "package::org.Pf.plasmo.db" %in% search() ) {
				library(org.Pf.plasmo.db)
			}
			db = org.Pf.plasmo.db

			# the annotation
			if (missing(annotation)) { annotation = "" }

		} else {
			return(1)
		}


		## output directory #####################################################
		dir <- file.path(pdf_dir, "ontology", tolower(analysis))
		if (!dir.exists(dir)) {
			dir.create(dir, recursive=T)
		}


		## presence of differentially-expressed genes checking ##################

		# print that if the GO analysis fails
		dfail <- data.frame(x=0, y=0, label="There is not enough genes...")
		gfail <- ggplot(dfail, aes(x=x, y=y, label=label)) + geom_label()

		# the differentially expressed genes
		genes <- rownames(results)

		# the code cannot be run if there is no gene
		if (length(genes)==0) {
			lab <- "There is no differentially epxressed genees.."
			dnogene <- data.frame(x=0, y=0, label=lab)
			gnogene <- ggplot(dnogene, aes(x=x, y=y, label=label)) + geom_label()
			print(gnogene)
			filename <- paste0(tolower(name), ".pdf")
			ggsave(file.path(dir, filename), plot=gnogene, width=19.2, height=12)
			return(NULL)
		}


		## plot parameters ######################################################

		# the main title of the final plot
		title <- paste0(" (", name, ", ", length(genes), " genes)")

		# the theme of the final plot
		theme <- theme(plot.title=element_text(hjust=.5), legend.position="top")


		## gene ontology ########################################################

		# cellular component GO analysis
		egocc <- clusterProfiler::enrichGO(genes,
													  OrgDb=db,
													  keytype=annotation,
													  ont="CC",
													  pAdjustMetho=adj_method,
													  pvalueCutoff=pvalue_threshold,
													  qvalueCutoff=qvalue_threshold,
													  readable=T)
		if (!is.null(egocc)) {
			cc <-
				ontology_dotplot(egocc) +
					ggtitle(paste0("Cellular Component", title)) + theme
		} else {
			cc <- gfail
		}

		# molecular function GO analysis
		egomf <- clusterProfiler::enrichGO(genes,
													  OrgDb=db,
													  keytype=annotation,
													  ont="MF",
													  pAdjustMetho=adj_method,
													  pvalueCutoff=pvalue_threshold,
													  qvalueCutoff=qvalue_threshold,
													  readable=T)
		if (!is.null(egomf)) {
			mf <-
				ontology_dotplot(egomf) +
					ggtitle(paste0("Molecular Function", title)) + theme
		} else {
			mf <- gfail
		}

		# biological process GO analysis
		egobp <- clusterProfiler::enrichGO(genes,
													  OrgDb=db,
													  keytype=annotation,
													  ont="BP",
													  pAdjustMetho=adj_method,
													  pvalueCutoff=pvalue_threshold,
													  qvalueCutoff=qvalue_threshold,
													  readable=T)
		if (!is.null(egobp)) {
			bp <-
				ontology_dotplot(egobp) +
					ggtitle(paste0("Biological Process", title)) + theme
		} else {
			bp <- gfail
		}


		## pathways analysis ####################################################

		# converts to ENTREZID because it's the only what is accepted by KEGG
		entrez <- clusterProfiler::bitr(genes,
												  fromType=annotation,
												  toType="ENTREZID",
												  OrgDb=db)
		entrez <- entrez[,"ENTREZID"]

		# sometimes need to add the key type depends on the organism
		if (tolower(binomial)=="hs") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 organism="hsa",
														 pvalueCutoff=qvalue_threshold)

		} else if (tolower(binomial)=="mm") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 organism="mmu",
														 pvalueCutoff=qvalue_threshold)

		} else if (tolower(binomial)=="dm") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 organism="dme",
														 keyType="ncbi-geneid",
														 pvalueCutoff=qvalue_threshold)

		} else if (tolower(binomial)=="sc") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 organism="sce",
														 pvalueCutoff=qvalue_threshold)

		} else if (tolower(binomial)=="xl") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 organism="xla",
														 pvalueCutoff=qvalue_threshold)

		} else if (tolower(binomial)=="pf") {

			kk <- clusterProfiler::enrichKEGG(gene=entrez,
														 #organism="pfd",
														 organism="pfa",
														 pvalueCutoff=qvalue_threshold)
		}

		if (!is.null(kk)) {
			kegg <-
				ontology_dotplot(kk) + ggtitle(paste0("KEGG Pathways", title)) +
					theme
		} else {
			kegg <- gfail
		}
		

		## PLOTTING #############################################################

		g <- gridExtra::arrangeGrob(bp, mf, cc, kegg, nrow=2)
		print(g)
		filename <- paste0(tolower(name), ".pdf")
		ggsave(file.path(dir, filename), plot=g, width=19.2, height=12)

		return(list(egocc, egomf, egobp, kk))
	}


###
get_rgm <-
	function(condition, dds) {
		rgm <- do.call(cbind,
							by(t(assay(dds)),
								colData(dds)[,condition],
								function (mat) {
									return(colMeans(mat))
								})
							)
		return(rgm)
	}


###
create_output_data_frame <-
	function(deglist, annot, mart) {

		# progressing info
		message <- "Create output data frame"
		cat(paste0("\n", message, "\n"))

		# all the differencially expressed genes whatever the contrast is
		genes <- unique( unlist( lapply( deglist, rownames ) ) )

		# the symbols of these genes
		symbols <- annot[match(genes, annot[,"gene_id"]), "gene_name"]

		# the left-sided base of the data frame
		dbase <- data.frame(gene=genes, symbol=symbols, round(rgm, 3)[genes,])

		# add a column for each contrast
		for (i in names(deglist)) {

			# log fold change
			dbase[rownames(deglist[[i]]), paste(i, "LFC", sep=".")] <-
				round(deglist[[i]][,"log2FoldChange"], 3)

			# p value
			dbase[rownames(deglist[[i]]), paste(i, "p.value", sep=".")] <-
				signif(deglist[[i]][,"pvalue"], 5)

			# q value
			dbase[rownames(deglist[[i]]), paste(i, "q.value", sep=".")] <-
				signif(deglist[[i]][,"padj"], 5)
		}

		# add the description for each as the last column
		if (length(genes)>0) {
			ensemblRes <- getBM(attributes=c("ensembl_gene_id", "description" ),
									  filters="ensembl_gene_id",
									  values=genes,
									  mart=mart)
			dbase[ensemblRes[,1],"description"] <- ensemblRes[,2]
		}

		# empty string the gene is not differentially expressed
		dbase[is.na(dbase)] <- ""

		return(dbase)
	}


export_xlsx <-
	# saves as an xlsx file
	function(df, analysis, contrast_count, summary, xlsx_dir) {

		# need to manipulate xlsx files
		if ( ! "package::openxlsx" %in% search() ) { library(openxlsx) }

		# the workbook
		wb <- createWorkbook("dgea")

		# the tabs
		addWorksheet(wb, "Differential genes")
		addWorksheet(wb, "Summary")
		addWorksheet(wb, "Session info")

		# the number of each worksheet
		diff_ws <- 1
		summ_ws <- 2
		sess_ws <- 3

		## the summary worksheet ################################################

		# just in case
		comparison_col <- 1
		genes_col <- 2

		# write the summary
		writeData(
					 wb			=	wb,
					 sheet		=	summ_ws,
					 x				=	summary,
					 startRow	=	1,
					 startCol	=	1,
					 colNames	=	T,
					 rowNames	=	F
					 )

		# header summary style
		headerSummaryStyle <-
			createStyle(
							halign			=	"center",
							textDecoration	=	c("italic", "bold", "underline")
							)
		addStyle(
					wb		=	wb,
					sheet	=	summ_ws,
					style	=	headerSummaryStyle,
					cols	=	1:ncol(summary),
					rows	=	1,
					stack	=	T
					)

		# comparison summary style
		comparisonSummaryStyle <-
			createStyle(
							halign			=	"left",
							textDecoration	=	c("italic")
							)
		addStyle(
					wb		=	wb,
					sheet	=	summ_ws,
					style	=	comparisonSummaryStyle,
					cols	=	comparison_col,
					rows	=	2:(nrow(summary)+1),
					stack	=	T
					)

		# genes summary style
		genesSummaryStyle <-
			createStyle(
							halign			=	"right",
							textDecoration	=	c("italic")
							)
		addStyle(
					wb		=	wb,
					sheet	=	summ_ws,
					style	=	genesSummaryStyle,
					cols	=	genes_col,
					rows	=	2:(nrow(summary)+1),
					stack	=	T
					)

		# adjust columns widths
		setColWidths(
						 wb		=	wb,
						 sheet	=	summ_ws,
						 cols		=	1:ncol(summary),
						 widths	=	rep("auto", ncol(summary))
						 )
		
		
		## the session worksheet ################################################
		
		# the session info as a two columns data frame
		writeData(
					 wb			=	wb,
					 sheet		=	sess_ws,
					 x				=	information.df,
					 startRow	=	1,
					 startCol	=	1,
					 colNames	=	F,
					 rowNames	=	F
					 )

		# the left column is the name of the information, this is its style
		infoNameStyle <-
			createStyle(
							halign			=	"left",
							textDecoration	=	c("italic", "bold", "underline")
							)

		# the right column is the the information it self, this is its style
		infoStyle <-
			createStyle(
							textDecoration	=	c("italic")
							)
		
		# add the style of the left column
		addStyle(
					wb		=	wb,
					sheet	=	sess_ws,
					style	=	infoNameStyle,
					cols	=	1,
					rows	=	1:nrow(information.df)
					)


		# add the style of the right column
		addStyle(
					wb		=	wb,
					sheet	=	sess_ws,
					style	=	infoStyle,
					cols	=	2,
					rows	=	1:nrow(information.df)
					)

		# everything is written in red
		lfcStyle <- createStyle(
										fontColour	=	"#ff0000"
										)
		addStyle(
					wb		=	wb,
					sheet	=	sess_ws,
					style	=	lfcStyle,
					cols	=	1:ncol(information.df),
					rows	=	nrow(information.df),
					stack	=	T
					)
		
		# adjust columns widths
		setColWidths(
						 wb		=	wb,
						 sheet	=	sess_ws,
						 cols		=	1:ncol(information.df),
						 widths	=	rep("auto", ncol(information.df))
						 )
		

		## the differentially expressed genes worksheet #########################
		
		# the worksheet at the extreme left top corner
		start_row <- 1
		start_col <- 1

		# there is three statistics for each condition (LFC, pvalue and qvalue)
		stat_count <- 3


		# the two first left columns are just gene id and the symbol of the gene
		idcount <- 2

		# the row range for the values
		values_row_range <- (start_row + 2) : ( nrow(df) + start_row + 1 )

		# the "Mean of reads" column range
		mean_col_range <- (idcount+1) : ( ncol(rgm) + (idcount+start_col-1) )
		
		# the first row will just contain the name of each contrast
		upper_row <- gsub("^(.*)\\.(LFC|p\\.value|q\\.value)$", "\\1", names(df))
		upper_row <- gsub("gene|symbol|description", "", upper_row)

		# in the first row, the cells cells for the mean of reads for each
		# condition will be fused
		upper_row[mean_col_range] <- "Mean of reads"
		upper_row <- t(upper_row)
		
		# the second row will contain the type of the information of column
		lower_row <- gsub("^.*\\.(LFC|p\\.value|q\\.value)$", "\\1", names(df))
		lower_row <- t(lower_row)

		# mean of reads colors
		mean_bg_color <- "#66a3ff"
		mean_fg_color <- "#1a75ff"

		# statistics colors
		diff_fg_color_lfc <- "#9999ff"
		diff_fg_color_pvalue <- "#4d4dff"
		diff_fg_color_qvalue <- "#0000b3"
		diff_bg_color <- "#b3b3ff"
		
		# write the upper row
		writeData(
					 wb			=	wb,
					 sheet		=	diff_ws,
					 x				=	upper_row,
					 startRow	=	start_row,
					 startCol	=	start_col,
					 colNames	=	F,
					 rowNames	=	F
					 )

		# write the lower row
		writeData(
					 wb			=	wb,
					 sheet		=	diff_ws,
					 x				=	lower_row,
					 startRow	=	start_row+1,
					 startCol	=	start_col,
					 colNames	=	F,
					 rowNames	=	F
					 )

		# write all the values
		writeData(
					 wb			=	wb,
					 sheet		=	diff_ws,
					 x				=	df,
					 startRow	=	start_row+2,
					 startCol	=	start_col,
					 colNames	=	F,
					 rowNames	=	F
					 )

		# fuse the "Mean of reads" cells in the first row
		mergeCells(
					  wb		=	wb,
					  sheet	=	diff_ws,
					  cols	=	mean_col_range,
					  row		=	start_row
					  )
		
		# the style of the fused cell "Mean of reads" in the first row
		meanUpperStyle <-
			createStyle(
							halign			=	"center",
							fontColour		=	"#ffffff",
							fgFill			=	mean_bg_color,
							textDecoration	=	"bold",
							border			=	"TopBottomLeftRight"
							)
		addStyle(
					wb		=	wb,
					sheet	=	diff_ws,
					style	=	meanUpperStyle,
					cols	=	mean_col_range,
					rows	=	start_row,
					stack	=	T
					)
		
		# the style of the condition names associated with "Mean of reads"
		# in the second row
		meanLowerStyle <-
			createStyle(
							halign			=	"center",
							fontColour		=	mean_fg_color,
							textDecoration	=	c("bold", "italic")
							) 
		addStyle(
					wb		=	wb,
					sheet	=	diff_ws,
					style	=	meanLowerStyle,
					cols	=	mean_col_range,
					rows	=	start_row+1,
					stack	=	T
					)
		
		# the general style of the second row
		lowerStyle <-
			createStyle(
							halign			=	"center",
							textDecoration	=	c("bold", "italic")
							)
		addStyle(
					wb		=	wb,
					sheet	=	diff_ws,
					style	=	lowerStyle,
					cols	=	start_col:( ncol(df)+start_col ),
					rows	=	start_row+1,
					stack	=	T
					)
		
		# the style of the values of "Mean of reads"
		meanValueStyle <-
			createStyle(
							fontColour		=	mean_fg_color,
							textDecoration	=	"italic"
							)
		addStyle(
					wb				=	wb,
					sheet			=	diff_ws,
					style			=	meanValueStyle,
					cols			=	mean_col_range,
					rows			=	values_row_range,
					stack			=	T,
					gridExpand	=	T
					)
		
		
		# need to apply the process for each contrast (comparison)
		for (i in 0:(contrast_count-1)) {
		
			# the starting and the ending columns for each contrast
			start_column <- start_col + idcount + ncol(rgm) + (stat_count * i)
			end_column <- start_column + 2
			contrast_col_range <- start_column:end_column

			# columns of each statistic
			lfc_col <- start_column
			pvalue_col <- start_column + 1
			qvalue_col <- start_column + 2

			# fuses the cells that inform about the contrast in the first row
			mergeCells(
						  wb		=	wb,
						  sheet	=	diff_ws,
						  cols	=	contrast_col_range,
						  rows	=	start_row
						  )
		
			# the style of the fused-cell about the contrast in the first row
			upperStyle <-
				createStyle(
								halign			=	"center",
								fontColour		=	"#ffffff",
								fgFill			=	diff_bg_color,
								textDecoration	=	"bold",
								border			=	"TopBottomLeftRight"
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	upperStyle,
						cols	=	contrast_col_range,
						rows	=	start_row,
						stack	=	T
						)
		
			# the style of the statitic name : LFC
			lowerStyle_lfc <-
				createStyle(
								halign			=	"center",
								fontColour		=	diff_fg_color_lfc,
								textDecoration	=	c("bold", "italic")
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	lowerStyle_lfc,
						cols	=	lfc_col,
						rows	=	start_row+1,
						stack	=	T
						)
		
			# the style of the statitic name : pvalue
			lowerStyle_pvalue <-
				createStyle(
								halign			=	"center",
								fontColour		=	diff_fg_color_pvalue,
								textDecoration	=	c("bold", "italic")
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	lowerStyle_pvalue,
						cols	=	pvalue_col,
						rows	=	start_row+1,
						stack	=	T
						)
		
			# the style of the statitic name : pvalue
			lowerStyle_qvalue <-
				createStyle(
								halign			=	"center",
								fontColour		=	diff_fg_color_qvalue,
								textDecoration	=	c("bold", "italic")
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	lowerStyle_qvalue,
						cols	=	qvalue_col,
						rows	=	start_row+1,
						stack	=	T
						)
		
			# the style of the statitic values : LFC
			valueStyle_lfc <-
				createStyle(
								fontColour		=	diff_fg_color_lfc,
								textDecoration	=	"italic"
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	valueStyle_lfc,
						cols	=	lfc_col,
						rows	=	values_row_range,
						stack	=	T
						)
		
			# the style of the statitic values : pvalue
			valueStyle_pvalue <-
				createStyle(
								fontColour		=	diff_fg_color_pvalue,
								textDecoration	=	"italic"
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	valueStyle_pvalue,
						cols	=	pvalue_col,
						rows	=	values_row_range,
						stack	=	T
						)
		
			# the style of the statitic values : qvalue
			valueStyle_qvalue <-
				createStyle(
								fontColour		=	diff_fg_color_qvalue,
								textDecoration	=	"italic"
								)
			addStyle(
						wb		=	wb,
						sheet	=	diff_ws,
						style	=	valueStyle_qvalue,
						cols	=	qvalue_col,
						rows	=	values_row_range,
						stack	=	T
						)
		}
		
		# adjust columns width
		setColWidths(
						 wb		=	wb,
						 sheet	=	diff_ws,
						 cols		=	1:( ncol(df) + start_col ),
						 widths	=	rep("auto", ncol(df))
						 )
		

		# save the xlsx file
		dir <- file.path(xlsx_dir, analysis)
		if (!dir.exists(dir)) { dir.create(dir, recursive=T) }
		saveWorkbook(
						 wb			=	wb,
						 file			=	file.path(dir, "diff_genes.xlsx"),
						 overwrite	=	T
						 )

		return(wb)
	}


###
venn_plot <-
	# creates a venn diagram, saves it and returns it
	function(genes, name, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# the packages for the venn diagrams
		if ( ! "package::VennDiagram" %in% search() ) { library(VennDiagram) }

		# the packages for the venn diagrams
		if ( ! "package::grid" %in% search() ) { library(grid) }

		## xxxxxxxx #############################################################

		# the title of the plot
		title <- toupper(sub(COMPARISON_SEPARATOR, "\nVS\n", name))

		# the categories
		cat <- gsub("\\.", "\n", names(genes))

		# the colors
		color_choice <- c("blue", "red", "yellow", "purple", "orange")
		colors <- color_choice[1:length(genes)]

		# create the venn diagram
		v <- venn.diagram(
								x=genes,
								fill=colors,
								category.names=cat,
								filename=NULL,
								main=title,
								main.cex=2,
								sub.cex=1,
								cat.default.pos=c("outer"),
								fontfamily="sans",
								main.fontfamily="sans",
								sub.fontfamily="sans",
								cat.fontfamily="sans",
								margin=.1
								)

		# the output file name
		tag <- tolower(sub(COMPARISON_SEPARATOR, "_vs_", name))
		dir <- file.path(pdf_dir, "comparison")
		if (!dir.exists(dir)) { dir.create(dir, recursive=T) }
		path <- file.path(dir, paste0(tag, ".pdf"))

		# save the venn diagram
		pdf(path)
		grid.draw(v)
		dev.off()

		return(v)
	}


###
grid_venn <-
	# puts a set a venn diagram on the same grid
	function(venns, pdf_dir=PDF_DIRPATH) {

		## packages #############################################################

		# grid package for the gTree method
		if ( ! "package::grid" %in% search() ) { library(grid) }

		# gridExtra package for the grid.arrange method
		if ( ! "package::gridExtra" %in% search() ) { library(gridExtra) }

		# ggplot2 package for the ggsave method
		if ( ! "package::ggplot2" %in% search() ) { library(ggplot2) }

		## xxxxxxxx #############################################################

		# need to transform venn diagram into gTree objects
		gtrees <- lapply(venns,
							  function(x) {
								  g <- grid::gTree(
														 children=x,
														 vp=viewport(width=.95, height=.95)
														 )
								  return(g)
							  })

		# grid structure and dimensions
		n <- length(venns)
		scale <- 6
		ncol <- ceiling(sqrt(n))
		nrow <- ceiling(n/ncol)

		# put the gTree objects on a grid
		g <- gridExtra::grid.arrange(grobs=gtrees, ncol=ncol)

		# save the file
		p <- file.path(pdf_dir, "comparison.pdf")
		ggsave(p, plot=g, width=ncol*scale, height=nrow*scale)

		return(g)
	}


###
create_regex <-
	# returns a regex from a vector of words
	function(words) {
		names <- gsub("^", "\\^", gsub("$", "\\$", words) )
		regex <- paste(names, collapse="|")
		return(regex)
	}

###############################################################################
##                                                                           ##
##                         ARGUMENTS PARSING                                 ##
##                                                                           ##
###############################################################################

###################
## THE ARGUMENTS ##
###################

## the directory that contains the STAR results files
results.opt <- make_option(
								  c("-r", "--results-dir"),
								  type="character",
								  default=NULL,
								  help="results files directory absolute path"
								  )

## the design file that maps sample names to conditions
design.opt <- make_option(
								  c("-d", "--design-file"),
								  type="character",
								  default=NULL,
								  help="design file absolute path"
								  )

## the binomial which is required for gene ontology analysis
binom.opt <- make_option(
							  c("-b", "--binomial"),
							  type="character",
							  default=NULL,
							  help="binomial"
							  )

## the design file that maps sample names to conditions
annot.opt <- make_option(
							  c("-a", "--annot-file"),
							  type="character",
							  default=NULL,
							  help="annotation file absolute path"
							  )

## the q value threshold for the differentially expressed genes
thresh.opt <- make_option(
							  c("-t", "--q-value-threshold"),
							  type="double",
							  default=0.05,
							  help="q value threshold for significant genes"
							  )

## the output directory
output.opt <- make_option(
								  c("-o", "--output-dir"),
								  type="character",
								  default=".",
								  help="output directory absolute path"
								  )


#############
## PARSING ##
#############

## create the options and parse them
options = list(results.opt,
					design.opt,
					binom.opt,
					annot.opt,
					thresh.opt,
					output.opt)
parser = OptionParser(option_list=options)
args = parse_args(parser)

## count the minimal number of argument
default.values <- lapply(options, function(x)as.character(attr(x, "default")))
min.arg.count <- length(which(unlist(default.values)=="\001NULL\001"))

## it has to be at least min.arg.count arguments
if (length(args)-1<min.arg.count) { # "args" contains also $help, so remove 1
	print_help(parser)
	stop()
}

## check the paths provided in the arguments
check_path(args)

## after checking, everything shoud be ok
RESULTS_DIRPATH <- args[["results-dir"]]
DESIGN_FILEPATH <- args[["design-file"]]
BINOMIAL <- args[["binomial"]]
ANNOT_FILEPATH <- args[["annot-file"]]
Q_VALUE_THRESHOLD <- args[["q-value-threshold"]]
OUTPUT_DIRPATH <- args[["output-dir"]]



###############################################################################
##                                                                           ##
##                       OUTPUT DIRECTORIES CREATION                         ##
##                                                                           ##
###############################################################################


# directories by file types
PDF_DIRNAME <- "pdf"
XLSX_DIRNAME <- "xlsx"
CSV_DIRNAME <- "csv"

# directories path
PDF_DIRPATH <- file.path(OUTPUT_DIRPATH, PDF_DIRNAME)
XLSX_DIRPATH <- file.path(OUTPUT_DIRPATH, XLSX_DIRNAME)
CSV_DIRPATH <- file.path(OUTPUT_DIRPATH, CSV_DIRNAME)

# creation
if (!dir.exists(PDF_DIRPATH)) { dir.create(PDF_DIRPATH, recursive=T) }
if (!dir.exists(XLSX_DIRPATH)) { dir.create(XLSX_DIRPATH, recursive=T) }
if (!dir.exists(CSV_DIRPATH)) { dir.create(CSV_DIRPATH, recursive=T) }

## very annoying
dummy.path <- file.path(PDF_DIRPATH, "dummy.pdf")
pdf(dummy.path)



###############################################################################
##                                                                           ##
##                INFORMATION, DATA AND ANNOTATION LOADING                   ##
##                                                                           ##
###############################################################################

CONDITION_SEPARATOR <- "_"
COMPARISON_SEPARATOR <- ".vs."

#CONDITION_COLNAME <- "condition"
#SAMPLEID_COLNAME <- "sample_id"



###############################################################################
##                                                                           ##
##                INFORMATION, DATA AND ANNOTATION LOADING                   ##
##                                                                           ##
###############################################################################

################################
## DESIGN, MAP AND ANNOTATION ##
################################

### the map ##
#
## the map associate each sample to its original name after sequencing
#map <- read.csv(MAP_FILEPATH, header=T)
#map <- as.data.frame(apply(map, 2, as.character))
#
## the map should be bijective and with good names
#if ( ncol(map)!=2 || names(map) != c("sample", "sample_id") ) {
#	#stop("Error: the map file is not valid.")
#}


## the design ##

# the design gives the condition or treatments for each sample
design <- read.csv(DESIGN_FILEPATH, header=T)
design <- as.data.frame(apply(design, 2, as.character))

# a global condition vector need to be added
names(design) <- sub("^condition$", "condition.exp", names(design))
info_design_names <- c("sample", "file", "file1", "file2")
info_design_regex <- create_regex(info_design_names)
exp.conditions <- names(design)[-grep(info_design_regex, names(design))]

# check the number of experimental conditions
if (length(exp.conditions)==0) {

	quit(save='n')

} else if (length(exp.conditions)>1) {

	design$condition <- apply(design[,exp.conditions],
									  1,
									  paste,
									  collapse=CONDITION_SEPARATOR)
} else {
	design$condition <- design[,exp.conditions]
}

# the way we name the samples can be change if we want
design[,"sample_id"] <- design[,"sample"]

## merge map with design
#design <- merge(design, map)

## genome annotation ##

# the annotation is loading as data frame
message <- "Load annotation"
cat(paste0("\n", message, "\n"))
annot <- mcols(import(ANNOT_FILEPATH))


##################
## DATA RESULTS ##
##################

# get the sorted file names
results.files <- list.files(
									 path=RESULTS_DIRPATH,
									 recursive=T,
									 full.name=T,
									 pattern=".*genes.results$"
									 )
results.files <- sort(results.files)

# import the results as either complete data frames or count data frames
results <- lapply(results.files, get_counts)
names(results) <- unlist(lapply(results, function(x)x[["df"]][1,"sample_id"]))

# from the results list, create a complete long data frame
df <- do.call(rbind, lapply(results, function(x)return(x$df)))
df <- merge(df, design)

# from the results list, create a wide data frame with the just the counts
count <- do.call(cbind, lapply(results, function(x)return(x$count)))

# remove 0s
count <- count[rowMeans(count)>0,]



###############################################################################
##                                                                           ##
##                        GENERAL QUALITY CONTROL                            ##
##                                                                           ##
###############################################################################

## MA plots between replicates ##
ddply(df, "condition", maplot)

## normalisation methods visualisation ##
normalise_df <- normalisation(count)



###############################################################################
##                                                                           ##
##                               DESeq2                                      ##
##                                                                           ##
###############################################################################

####################
## GENERAL MATRIX ##
####################

# to be sure that each count matrix column is properly associated to its
# condition
count_cond <- design[ match(names(count), design[,"sample_id"]) , "condition" ]
count_cond <- as.factor(count_cond)

# progressing info
cat("\nCreate the general DESeq2 object\n")

# the general DESeq2 matrix
dds <- DESeqDataSetFromMatrix(countData=round(count),
										colData=DataFrame(condition=count_cond),
										design=as.formula("~ condition"))
dds <- DESeq(dds)

# transforms
rld <- rlog(dds, blind=T)
vst <- vst(dds, blind=T)

## replicate group means
rgmL <- lapply("condition", get_rgm, rld)
rgm <- do.call(cbind, rgmL)


# plot principal component analysis
p <-
	plot_pca(vst, design[,-grep(info_design_regex, names(design))], PDF_DIRPATH)

# plot the euclidian distance
hmap <- plot_heatmap(vst,
							design[,-grep(info_design_regex, names(design))],
							PDF_DIRPATH)

# check the best transform with the mean-sd plots
meansdgrob <- meansdplot(dds, PDF_DIRPATH)

##########################
## CONTRASTS GENERATION ##
##########################

# all the different conditions of the experiment
conditions <- sort(unique(design$condition))

# the number of comparisons
contrasts.count <- choose(length(conditions), 2)


# the list that contains all the comparison as a size 2 vector of conditions
contrasts <- combn(conditions, 2, simplify=F)
contrasts <- lapply(contrasts, sort)
names(contrasts) <- sapply(contrasts, paste, collapse=COMPARISON_SEPARATOR)
contrasts <- contrasts[sort(names(contrasts))]


##################################
## DESEQ DATA SET LIST CREATION ##
##################################

# for each comparison individually
ddsL.general <- lapply(contrasts, function(x)return(dds))
ddsL.individual <- lapply(contrasts, create_dds, count, design)

# get results for each comparison
resL.general <- mapply(get_results_from_dds,
								  ddsL.general,
								  contrasts[names(ddsL.general)])
resL.individual <- mapply(get_results_from_dds,
								  ddsL.individual,
								  contrasts[names(ddsL.individual)])

# get the significant results
decideL.general <- lapply(resL.general,
								  get_significant_results,
								  Q_VALUE_THRESHOLD)
decideL.individual <- lapply(resL.individual,
								  get_significant_results,
								  Q_VALUE_THRESHOLD)
decideL <- list(general=decideL.general, individual=decideL.individual)



###############################################################################
##                                                                           ##
##                            GENE ONTOLOGY                                  ##
##                                                                           ##
###############################################################################

# the binomial of the organism
binomial <- paste(substr(unlist(strsplit(BINOMIAL, " ")), 1, 1), collapse="")

# gene ontology and pathways
go <- mapply(function(analysis, analysis.name) { 

				 # for each case : general and individual
				 go_analysis <- mapply(function(results, contrast.name) {

											  # perform ontolgy and plot it
											  ontology(results=results,
														  analysis=analysis.name,
														  name=contrast.name,
														  binomial=binomial)

														 }, analysis, names(analysis))

								  }, decideL, names(decideL))



###############################################################################
##                                                                           ##
##                             RESULTS OUTPUT                                ##
##                                                                           ##
###############################################################################

###################
## BIOMART SETUP ##
###################

# binomial choice
if (tolower(BINOMIAL)=="homo sapiens") {

	mart.biomart <- "ensembl"
	mart.dataset <- "hsapiens_gene_ensembl"

} else if (tolower(BINOMIAL)=="mus musculus") {

	mart.biomart <- "ensembl"
	mart.dataset <- "mmusculus_gene_ensembl"

} else if (tolower(BINOMIAL)=="drosophila melanogaster") {

	mart.biomart <- "ensembl"
	mart.dataset <- "dmelanogaster_gene_ensembl"

} else if (tolower(BINOMIAL)=="saccharomyces cerevisiae") {

	mart.biomart <- "ensembl"
	mart.dataset <- "scerevisiae_gene_ensembl"

} else if (tolower(BINOMIAL)=="xenopus laevis") {

} else if (tolower(BINOMIAL)=="plasmodium falciparum") {

}

# the actual setup
mart <- useEnsembl(biomart=mart.biomart,
						 dataset=mart.dataset,
						 version=ENSEMBL_VERSION)


###############################
## FINAL DATA FRAME CREATION ##
###############################

dfL <- lapply(decideL, create_output_data_frame, annot, mart)


###################
## CSV EXPORTING ##
###################

mapply(
		 function(df, name) {

		 # progressing info
		 message <- paste("Save CSV file for", name)
		 cat(paste0("\n", message, "\n"))

		 # output directory
		 dir <- file.path(CSV_DIRPATH, name)
		 if (!dir.exists(dir)) { dir.create(dir, recursive=T) }

		 # save csv
		 write.csv(df, file.path(dir, "diff_genes.csv"), row.names=F)

						 }, dfL, names(dfL)
		 )


####################
## XLSX EXPORTING ##
####################

mapply(
		 function(df, name) {

		 # progressing info
		 message <- paste("Save XLSX file for", name)
		 cat(paste0("\n", message, "\n"))

		 # create summary
		 decide <- decideL[[name]]
		 summary <- as.data.frame( unlist(lapply(decide, nrow)) )
		 comparison <- sub(COMPARISON_SEPARATOR, " VS ", rownames(summary))
		 genes <- summary[,1]
		 summary <- data.frame(comparison=comparison, genes=genes)
		 summary <- summary[order(summary[,"genes"], decreasing=T),]

		 # save csv
		 wb <- export_xlsx(df, name, length(contrasts), summary, XLSX_DIRPATH) 

						 }, dfL, names(dfL)
		 )


######################################
## GENES COMMON TO THE TWO ANALYSIS ##
######################################

## check if the contrasts are the same in the two analysis ##

# general
general.contrasts <- sort( names( decideL[["general"]] ) )

# individual
individual.contrasts <- sort( names( decideL[["individual"]] ) )

# quit if not the same
if ( all( general.contrasts!=individual.contrasts ) ) {
	message <- "Error: general and individual don't have the same contrasts"
	stop(message)
}

## the list of the comon genes ##

# create a list of element that contains a list of the differentially genes
# for each analysis
common <- lapply(
					  general.contrasts,
					  function(contrast) {

						  # general
						  g <- rownames(decideL[["general"]][[contrast]])
						  g <- sort( unique( g ) )

						  # individual
						  i <- rownames(decideL[["individual"]][[contrast]])
						  i <- sort( unique( i ) )

						  return(list(general=g, individual=i))
					  })
names(common) <- general.contrasts

## the venn diagrams ##

# create venn diagrams for each contrast and save it
vennL <- mapply(
					 function(genes, name) {

						 total_genes_num <- sum(unlist( lapply(genes, length) ))

						 if (total_genes_num>0) {
							 v <- venn_plot(genes, name, PDF_DIRPATH)
							 return(v)
						 }

					 }, common, names(common)
					 , SIMPLIFY=F
					 )

# put all the venn diagram on the same grid and save the grid
venn.grid <- grid_venn(vennL, PDF_DIRPATH)

# venn diagrams creation generates log files to remove
file.remove(list.files(".", pattern="VennDiagram.*\\.log"))



## very annoying
dev.off()
file.remove(dummy.path)
