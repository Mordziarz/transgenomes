install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/transgenomes')
library(transgenomes)

library(circlize)
library(metablastr)
library(Biostrings)
library(dplyr)
library(rstatix)
library(ggplot2)
library(car)
library(ggpubr)

setwd("path/to/transgenomes/dir")

bed_q <- read.csv("/inst/extdata/mito_carrot.gff3",sep="\t",header = F)
bed_s <- read.csv("/inst/extdata/plas_carrot.gff3",sep="\t",header = F)

bed_q <- bed_q[bed_q$V3 %in% c("gene","pseudogene"),]
bed_q$V9 <- gsub(".*Name=([^;]+).*", "\\1", bed_q$V9)
bed_q <- bed_q[,c("V1","V4","V5","V9")]
colnames(bed_q) <- c("V1","V2","V3","V4")

bed_s <- bed_s[bed_s$V3 %in% c("gene","pseudogene"),]
bed_s$V9 <- gsub(".*Name=([^;]+).*", "\\1", bed_s$V9)
bed_s <- bed_s[,c("V1","V4","V5","V9")]
colnames(bed_s) <- c("V1","V2","V3","V4")


transfer_function_out <- transgenomes::transfer_function(fasta_mt = "/inst/extdata/mito_carrot.fasta",
                                              fasta_pt = "/inst/extdata/plastome_carrot.fasta",
                                              bed_mt = bed_q ,
                                              bed_pt = bed_s,
                                              evalue_cut_off = 0.000001,
                                              min_length = 40,
                                              min_identity = 70,
                                              gene_buffer = 20,
                                              trans_buffer = 20)


circlize::circos.clear()
png("Ex_plot_transfers.png", width=6, height=6, units = "in", res = 300)
transgenomes::plot_transfers(transfer_function_out = transfer_function_out)
dev.off()

transfer_function_out <- transfer_function_out[transfer_function_out$direction %in% c("PT -> MT","MT -> PT"),]

test_s <- transgenomes::extract_regions(transfer_function_out = transfer_function_out,
                                        fasta_path = "/inst/extdata/plastome_carrot.fasta")

test_q <- transgenomes::extract_regions(transfer_function_out = transfer_function_out,
                                        fasta_path =  "/inst/extdata/mito_carrot.fasta")



results_GC <- transgenomes::analyze_GC_content(extract_regions_1 = test_s, 
                                               extract_regions_2 = test_q,
                                               group1_name = "Plastome",
                                               group2_name = "Mitogenome",
                                               alpha = 0.05,
                                               group1_col = "darkgreen", 
                                               group2_col = "orange")

results_GC$plot
results_GC$test_result
results_GC$normality_check
results_GC$test_method
results_GC$gc_data

png("Ex_GC_content.png", width=5, height=5, units = "in", res = 300)
results_GC$plot
dev.off()
