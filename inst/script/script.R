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

transfer_1 <- transgenomes::transfer_function(fasta_mt = "inst/extdata/mitogenome.fasta",
                                              fasta_pt = "inst/extdata/plastome.fasta",
                                        bed_mt = bed_q ,bed_pt = bed_s,evalue_cut_off = 0.00001)

circlize::circos.clear()
png("Ex_circos.png", width=6, height=6, units = "in", res = 300)
transgenomes::plot_transfers(transfer_function_out = transfer_1)
dev.off()

test_s <- transgenomes::extract_regions(transfer_function_out = transfer_1,
                              fasta_path = "inst/extdata/plastome.fasta",
                              id_col = "subject_id",
                              start_col = "s_start",
                              end_col = "s_end")

test_q <- transgenomes::extract_regions(transfer_function_out = transfer_1,
                              fasta_path = "inst/extdata/mitogenome.fasta",
                              id_col = "query_id",
                              start_col = "q_start",
                              end_col = "q_end")



results_GC <- transgenomes::analyze_GC_content(extract_regions_1 = test_s, extract_regions_2 = test_q,group1_name = "Plastome",group2_name = "Mitogenome",alpha = 0.05)

results_GC$plot
results_GC$test_result
results_GC$normality_check
results_GC$test_method
results_GC$gc_data

png("Ex_GC.png", width=5, height=5, units = "in", res = 300)
results_GC$plot
dev.off()

#### Transfers with nuclear genome
tabela <- transfer_function_nuclear(fasta_nuc = "path/to/nuclear/fasta",
                                         fasta_mt ="/path/to/mitochondrion/fasta" ,
                                         fasta_pt = "path/to/plastid/fasta",
                                         bed_nuc = bed_nuc,
                                         bed_mt = bed_mt,
                                         bed_pt = bed_pt,
                                         min_length = 100,
                                         evalue_cut_off = 0.000001,
                                         min_identity = 70,
                                         gene_buffer = 20,
                                         trans_buffer = 20)


regions_nuc <- extract_regions(transfer_function_out = tabela,
                             fasta_path = "path/to/nuclear/fasta")

regions_mt <- extract_regions(transfer_function_out = tabela,
                               fasta_path = "/path/to/mitochondrion/fasta")

regions_pt <- extract_regions(transfer_function_out = tabela,
                               fasta_path = "path/to/plastid/fasta")


analyze_GC_three_genomes_our <- analyze_GC_three_genomes(extract_regions_mt = regions_mt,
                                                         extract_regions_pt = regions_pt,
                                                         extract_regions_nuc=regions_nuc)