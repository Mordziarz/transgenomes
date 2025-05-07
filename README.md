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
| Mitogenome  | 0  | 74  | trnD1 |
| Mitogenome  | 270  |   342 |   trnN    |

```r
transfer_function(fasta_q = "fasta_q.fasta",
                  fasta_s = "fasta_s.fasta",
                  bed_q = bed_q,
                  bed_s = bed_s,
                  evalue_cut_off = 0.0001,
                  cores = 10)
```

The transfer_function() function produced an output table enriched with annotations for genes that had undergone partial or complete transfer to the second genome. The q_genes and s_genes columns within this table provided details about these transfers. For instance, the entry "genes: atp1 (1530)/(790)" indicated that the gene atp1, with a total length of 1530 nucleotides, had been transferred, and the transferred portion was 790 nucleotides long.

# Visualization

The program generated a basic visualization using the circlize package (https://github.com/jokergoo/circlize). This visualization was created based on the output from the transfer_function() function.

```r
plot_transfers(transfer_function_out = transfer_function_out,
                gap=40,
                start_degree=90,
                transparency=0.7,
                color="green")
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
