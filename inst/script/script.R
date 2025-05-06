install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/transgenomes')
library(transgenomes)

library(circlize)
library(metablastr)

setwd("path/to/transgenomes/dir")

bed_q <- read.csv("inst/extdata/mitogenome.gff3",sep="\t",header = F)
bed_s <- read.csv("inst/extdata/plastome.gff3",sep="\t",header = F)

bed_q <- bed_q[bed_q$V3 %in% c("gene","pseudogene"),]
bed_q$V9 <- gsub(".*Name=([^;]+).*", "\\1", bed_q$V9)
bed_q <- bed_q[,c("V1","V4","V5","V9")]
colnames(bed_q) <- c("V1","V2","V3","V4")

bed_s <- bed_s[bed_s$V3 %in% c("gene","pseudogene"),]
bed_s$V9 <- gsub(".*Name=([^;]+).*", "\\1", bed_s$V9)
bed_s <- bed_s[,c("V1","V4","V5","V9")]
colnames(bed_s) <- c("V1","V2","V3","V4")

table(bed_s$V4)

transfer_1 <- transgenomes::transfer_function(fasta_q = "inst/extdata/mitogenome.fasta",
                                              fasta_s = "inst/extdata/plastome.fasta",
                                        bed_q = bed_q ,bed_s = bed_s,evalue_cut_off = 0.00001)

circlize::circos.clear()
png("Ex_circos.png", width=5, height=5, units = "in", res = 300)
transgenomes::plot_transfers(transfer_function_out = transfer_1,color = "green4")
dev.off()