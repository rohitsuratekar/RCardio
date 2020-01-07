library("tximeta")
library("SummarizedExperiment")
library("DESeq2")
library("data.table")

source("local.R") # All the local variables will be imported from this file
names <-  c("SRX4720625","SRX4720626",
            "SRX4720628", "SRX4720629",
            "SRX4720631", "SRX4720632")

files <- file.path(FOLDER_RNA_SEQ, names, "salmon", "quant.sf") 
coldata <- data.frame(names=names, files=files, stringsAsFactors=FALSE)
indexDir <- FOLDER_SALOM_INDEX
gtfPath <- file.path(FOLDER_GTF, "Danio_rerio.GRCz11.98.chr.gtf")
fastaPath <-file.path(FOLDER_TRANSCRIPT, "Danio_rerio.GRCz11.cdna.all.fa.gz")

tximeta::makeLinkedTxome(indexDir=indexDir,
                source="Ensembl",
                organism="Danio rerio",
                release="98",
                genome="GTCz11.98",
                fasta=fastaPath,
                gtf=gtfPath,
                write=FALSE)

se <- tximeta::tximeta(coldata)

gse <- tximeta::summarizeToGene(se)

gse$time <- rep(c("24hpf", "48hpf", "72hpf"), each=2)

dds <- DESeq2::DESeqDataSet(gse, design = ~ time)

dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~1)
res <- DESeq2::results(dds)

head(res, 2)

DESeq2::resultsNames(dds) # This will give Coefficient Name to fill in following command

coef_name <-  "time_72hpf_vs_24hpf"
resLFC <- DESeq2::lfcShrink(dds, coef=coef_name, type="apeglm")



# Plotting MA plot
plotMA(resLFC, ylim = c(-10, 10),  main=coef_name )

# Saving data to the file
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
data.table::setDT(normalized_counts, keep.rownames="genes")
write.csv(normalized_counts, file="normalized_counts.csv")

df <- as.data.frame(res)
data.table::setDT(df, keep.rownames="genes")
write.csv(df, file="r_results.csv")
rld = DESeq2::rlogTransformation(dds)

DESeq2::plotPCA(rld, intgroup=c("time"))

