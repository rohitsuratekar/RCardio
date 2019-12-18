library("tximeta")
library("SummarizedExperiment")
library("DESeq2")

source("local.R") # All the local variables will be imported from this file
names <-  c("SRX4720625","SRX4720626","SRX4720628", "SRX4720629")

files <- file.path(FOLDER_RNA_SEQ, names, "salmon", "quant.sf") 
coldata <- data.frame(names=names, files=files, stringsAsFactors=FALSE)
indexDir <- FOLDER_SALOM_INDEX
gtfPath <- file.path(FOLDER_GTF, "Danio_rerio.GRCz11.98.chr.gtf")
fastaPath <-file.path(FOLDER_TRANSCRIPT, "Danio_rerio.GRCz11.cdna.all.fa.gz")

makeLinkedTxome(indexDir=indexDir,
                source="Ensembl",
                organism="Danio rerio",
                release="98",
                genome="GTCz11.98",
                fasta=fastaPath,
                gtf=gtfPath,
                write=FALSE)

se <- tximeta(coldata)

gse <- summarizeToGene(se)

dds <- DESeqDataSet(gse)

gse$time <- rep(c("24hpf", "48hpf"), each=2)

dds <- DESeqDataSet(gse, design = ~ time)

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds) # This will give Coefficient Name to fill in following command
resLFC <- lfcShrink(dds, coef="time_48hpf_vs_24hpf", type="apeglm")
plotMA(resLFC, ylim = c(-10, 10))
write.csv(as.data.frame(res), file="results.csv")
