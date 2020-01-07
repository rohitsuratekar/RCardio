# Count matrix
#
# Usage : Rscript countmatrix.R BAM_FILE_PATH GTF_FILE_PATH IS_PAIRED_END
# Example: Rscript countmatrix.R /path/to/SRA.Aligned.sortedByCoord.out.bam
# /path/to/Danio_rerio.GRCz11.98.chr.gtf
#
# Arguments
# BAM_FILE_PATH: Path to sorted BAM file
# GTF_FILE_PATH: Path to GTF file

library("optparse") # Need for getting proper arguments
library("Rsubread") # Need for featureCount function
library("data.table") # To properly format output file

# Generate options
opt_bam <- optparse::make_option(c("-f", "--file"),
help = "Full path of input BAM file.")

opt_gtf <- optparse::make_option(c("-g", "--gtf", help = "Full
path of the GTF file."))

opt_sra <- optparse::make_option(c("-s", "--srr"), help = "SRR ID which will be
 used to make column name in the count matrix. If not provided, input
 filename will be used", default = NULL)

opt_out <- optparse::make_option(c("-o", "--out"), help = "Full path of the 
output file.", default = "countmatrix.csv")

opt_paired <- optparse::make_option(c("-p", "--paired"), type = "integer",
help = "Is current input file is paired-end? Use 1 for yes [default] and 0 for no
.", default = 1)

opt_parser <- optparse::OptionParser(option_list = list(opt_bam, opt_gtf,
opt_paired, opt_sra, opt_out))

# Parse the options
opt <- optparse::parse_args(opt_parser)

paried_end <- TRUE
if (opt$paired == 0) {
    paired_end <- FALSE
}

fc <- Rsubread::featureCounts(
opt$file,
annot.ext = opt$gtf,
isGTFAnnotationFile = TRUE,
isPairedEnd = paried_end)

df <- as.data.frame(fc$counts)

# If SRR ID is provided, change the column name.
# Else, input filename will be used as a column name
if (! is.null(opt$srr)) {
    colnames(df) <- c(opt$srr)
}

# Set index column name before saving to the file 
data.table::setDT(df, keep.rownames = "gene_id")

# Save to the file
write.csv(df, file = opt$out)

