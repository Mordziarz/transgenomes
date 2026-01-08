# transgenomes

A BLASTn-based tool for pairwise genomes comparison (https://github.com/drostlab/metablastr).

# Installation

```r
install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/transgenomes')
library(transgenomes)
```
To use all features of the program, you will need several libraries.

```r
library(metablastr)
library(circlize)
library(Biostrings)
library(dplyr)
library(rstatix)
library(ggplot2)
library(car)
library(ggpubr)
```

# Input data 

The transfer_function() function accepted two FASTA files (fasta_q and fasta_s) along with two BED files containing gene annotations (bed_q and bed_s). The evalue_cutoff parameter was used to specify the E-value threshold for BLAST results.

The bed should look like this: 
V1 - Genome name,
V2 - Start annotation,
V3 - End annotation,
V4 - Gene name

| V1  | V2 | V3 |   V4  |
| -------- | ----- |    -----   |   -----   |
| OR220799.1  | 0  | 74  | trnD1 |
| OR220799.1  | 270  |   342 |   trnN    |

```r
transfer_function_out <- transfer_function(fasta_mt = "fasta_q.fasta",
                                            fasta_pt = "fasta_s.fasta",
                                            bed_mt = bed_q,
                                            bed_pt = bed_s,
                                            evalue_cut_off = 0.0001 
                                            gene_buffer = 20, 
                                            trans_buffer = 20)
```

The transfer_function() identifies intergenomic transfers by comparing BLASTn alignments against Mitochondrial (MT) and Plastid (PT) feature sets. For each alignment, the function calculates gene completeness (g_perc) and transfer dominance (t_perc) for both genomes, recorded in the mt_genes and pt_genes columns (e.g., atp1 (g_len=1530, ov_len=790, g_perc=51.6%, t_perc=98.5%)). Crucially, the algorithm ignores tRNA genes (filtered via regex ^trn|tRNA) when determining direction, as these highly conserved sequences are often non-diagnostic for transfer events. The final direction—stored in the direction column as MT -> PT, PT -> MT, or unknown—is determined by comparing the maximum metrics from both sides. A transfer is assigned if the difference in completeness exceeds the gene_buffer or if the difference in dominance exceeds the trans_buffer. All decisions, including cases where regions are intergenic or metrics are too ambiguous to reach a conclusion, are fully documented in a dedicated reason column.

# Visualization

The program generated a basic visualization using the circlize package (https://github.com/jokergoo/circlize). This visualization was created based on the output from the transfer_function() function.

```r
plot_transfers(transfer_function_out = transfer_function_out)
```

Useful functions for image cleaning in R

```r
dev.off()
circlize::circos.clear()
```

![Circular](inst/graphs/Ex_circos.png)

# Extract transfer regions

The extract_regions() function allows you to extract FASTA sequences from transfer events. Simply provide the output of transfer_function() as the transfer_function_out argument, and specify the appropriate column names for IDs and coordinates.

```r
regions_s <- extract_regions(transfer_function_out = transfer_1,
                              fasta_path = "inst/extdata/plastome.fasta",
                              id_col = "subject_id",
                              start_col = "s_start",
                              end_col = "s_end")
```

# GC content

You can easily compare the GC content between two sets of transfer regions using the analyze_GC_content() function. Group names and significance level can be customized. The function automatically selects the appropriate statistical test based on data distribution and provides a publication-ready plot with p-value annotation, as well as detailed test results.

```r
results_GC <- analyze_GC_content(extract_regions_1 = regions_s, 
                                extract_regions_2 = regions_q,
                                group1_name = "Plastome",
                                group2_name = "Mitogenome",
                                alpha=0.05)

results_GC$plot
results_GC$test_result
results_GC$normality_check
results_GC$test_method
results_GC$gc_data
```

![GC](inst/graphs/Ex_GC.png)

# Citation

Paper in preparation

# Support
Any issues connected with the transgenomes should be addressed to Mateusz Mazdziarz (mazdziarzm@gmail.com).

# Usage in scientific papers
