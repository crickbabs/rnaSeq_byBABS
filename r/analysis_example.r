
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

## after checking, everything shoud be ok
RESULTS_DIRPATH <- args[["results-dir"]]
DESIGN_FILEPATH <- args[["design-file"]]
BINOMIAL <- args[["binomial"]]
ANNOT_FILEPATH <- args[["annot-file"]]
Q_VALUE_THRESHOLD <- args[["q-value-threshold"]]
OUTPUT_DIRPATH <- args[["output-dir"]]


###############################################################################
##                                                                           ##
##                               DESeq2                                      ##
##                                                                           ##
###############################################################################

## very annoying
dummy.path <- file.path(OUTPUT_DIRPATH, "dummy.pdf")
pdf(dummy.path)

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

# from the results list, create a wide data frame with the just the counts
count <- do.call(cbind, lapply(results, function(x)return(x$count)))

# remove 0s
count <- count[rowMeans(count)>0,]

# design table
design <- read.csv(DESIGN_FILEPATH, header=T)
design <- as.data.frame(apply(design, 2, as.character))
design <- design[, c("sample", "treatment")]
condition <- design[ match(names(count), design[,"sample"]) , "treatment" ]
condition <- as.factor(condition)

# the DESeq2 matrix
dds <-
	DESeq2::DESeqDataSetFromMatrix(countData=round(count),
											 colData=DataFrame(condition=condition),
											 design=as.formula("~ condition"))

# estimating size factors and dispersion
dds <- DESeq2::DESeq(dds)

# the differentially expressed genes
results <- DESeq2::results(dds,
									alpha=.05,
									contrast=c("condition", "treated", "untreated"))
write.csv(results, file.path(OUTPUT_DIRPATH, "diff_genes.csv"))


# the pca plot
pca <- plotPCA(rlog(dds), intgroup="condition", ntop=1000)
pdf(file.path(OUTPUT_DIRPATH, "pca_plot.pdf"))
print(pca)
dev.off()
png(file.path(OUTPUT_DIRPATH, "pca_plot.png"), type="cairo")
print(pca)
dev.off()

# gene ontology for biological process
ontology <- clusterProfiler::enrichGO(rownames(results),
											  OrgDb="org.Hs.eg.db",
											  keytype="ENSEMBL",
											  ont="BP",
											  pAdjustMetho="BH",
											  pvalueCutoff=.01,
											  qvalueCutoff=.05,
											  readable=T)

# dot plot for gene ontology
dplot <- DOSE::dotplot(ontology)
dplot$data$Description <- substr(dplot$data$Description, 0, 50)
dplot$data$Description <- make.unique(dplot$data$Description, sep = "_")
idx <- order(dplot$data$GeneRatio, decreasing=F)
dplot$data$Description <-
	factor(dplot$data$Description, levels=dplot$data$Description[idx])
ggsave(file.path(OUTPUT_DIRPATH, "ontology_biological_process.pdf"), dplot)
ggsave(file.path(OUTPUT_DIRPATH, "ontology_biological_process.png"), dplot,
		 type="cairo")

## very annoying
dev.off()
file.remove(dummy.path)

