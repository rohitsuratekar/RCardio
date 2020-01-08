# Title     : Metadata Generator
# Objective : Generate the meta-data and input for DESeq2 analysis
# Created by: Rohit Suratekar
# Created on: 1/7/20

library("optparse") # Need for getting proper arguments

# Generate Options
opt_sra_file <- optparse::make_option(c("-f", "--file"),
help = "Full path of the file containing all the SRA details.")

opt_sra_folder <- optparse::make_option(c("-s", "--sra"), help = "Full path of
 the folder containing all the SRA folders of CardioTrans pipeline")

opt_parser <- optparse::OptionParser(option_list = list(opt_sra_file))

# Parse the options
opt <- optparse::parse_args(opt_parser, positional_arguments = TRUE)

# Check if SRA ids have provided or not
if (length(opt$args) == 0) {
    stop("SRA ids have not provided for the analysis. Please pass all the
    SRA ids at the end of all arguments. e.g. RScript generate_meta.R -f
    /path/to/sra_file SRA_ID1 SRA_ID2 ...")
}
# Read the SRA file and keep only entries as provided for the analysis
sra_file <- read.csv(opt$options$file, header = TRUE)
df <- as.data.frame(sra_file)
df <- df[which(df$sra_id %in% opt$args),]


print(opt)
