#' Annotate BLASTn alignments with overlapping gene features
#'
#' Performs nucleotide BLAST alignment, adjusts coordinates for directionality,
#' and annotates results with overlapping gene features from BED files.
#'
#' @param fasta_1 Path to the first subject FASTA file (required)
#' @param bed_1 The first subject species BED file (required)
#' @param fasta_2 Path to the second subject FASTA file (required)
#' @param bed_2 The second subject species BED file (required)
#' @param fasta_3 Path to the third subject FASTA file (required)
#' @param bed_3 The third subject species BED file (required)
#' @param evalue_cut_off Maximum e-value threshold for BLAST hits (default: 0.001)
#'
#' @return A data frame containing:
#' - BLAST alignment details
#' - Annotated overlapping genes from query and subject
#' - Relative positions of overlaps in format: "gene (full_length/overlap_length)"
#' 
#' @details
#' ## Requirements
#' - All input files must exist and be non-empty
#' - BED files must contain at least 4 columns (chrom, start, end, name)
#' - Requires `metablastr` package
#'
#' @examples
#' \dontrun{
#' result <- transfer_function(
#'   fasta_1 = "fasta_1",
#'   fasta_2 = "fasta_2",
#'   fasta_2 = "fasta_2",
#'   bed_1 = "bed_1",
#'   bed_2 = "bed_2",
#'   bed_3 = "bed_3",
#'   evalue_cut_off = 1e-5
#' )
#' }
#' @export

transfer_function2 <- function(fasta_1="",fasta_2="",fasta_3="",bed_1=bed_1, bed_2=bed_2, bed_3=bed_3, evalue_cut_off=0.001) {
  
  if (base::missing(fasta_1)) {
    stop("The fasta_1 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(fasta_2)) {
    stop("The fasta_2 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(fasta_3)) {
    stop("The fasta_3 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(bed_1)) {
    stop("The bed_1 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(bed_2)) {
    stop("The bed_2 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(bed_3)) {
    stop("The bed_3 predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  blast_1 <- metablastr::blast_nucleotide_to_nucleotide(query = fasta_1,
                                                        subject = fasta_2,
                                                        db.import = F,
                                                        task="blastn",
                                                        evalue = evalue_cut_off)
  
  blast_1$q_start <- as.numeric(blast_1$q_start)
  blast_1$s_start <- as.numeric(blast_1$s_start)
  
  blast_1$q_end <- as.numeric(blast_1$q_end)
  blast_1$s_end <- as.numeric(blast_1$s_end)
  
  blast_1$q_start1 <- ifelse(blast_1$q_start < blast_1$q_end, blast_1$q_start,blast_1$q_end)
  blast_1$q_end1 <- ifelse(blast_1$q_start < blast_1$q_end, blast_1$q_end,blast_1$q_start)
  
  blast_1$s_start1 <- ifelse(blast_1$s_start < blast_1$s_end, blast_1$s_start, blast_1$s_end)
  blast_1$s_end1 <- ifelse(blast_1$s_start < blast_1$s_end, blast_1$s_end, blast_1$s_start)
  
  
  blast_1$q_start <- blast_1$q_start1
  blast_1$s_start <- blast_1$s_start1
  
  blast_1$q_end <- blast_1$q_end1
  blast_1$s_end <- blast_1$s_end1
  
  blast_1$q_start1 <- NULL
  blast_n$s_start1 <- NULL
  
  blast_1$q_end1 <- NULL
  blast_1$s_end1 <- NULL
  
  blast_1$q_genes <- "genes: "
  blast_1$s_genes <- "genes: "
  
  bed_1 <- bed_1
  bed_2 <- bed_2
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_1$q_start[j] < bed_1$V2[i] & blast_1$q_end[j] > bed_1$V3[i]) {
        blast_1$q_genes[j] <- paste0(blast_1$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_1$q_start[j] > bed_1$V2[i] & blast_1$q_start[j] < bed_1$V3[i] & blast_1$q_end[j] > bed_1$V3[i]) {
        blast_1$q_genes[j] <- base::paste0(blast_1$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",bed_1$V3[i]-blast_1$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_1$q_start[j] < bed_1$V2[i] & bed_1$V2[i] < blast_1$q_end[j] & blast_1$q_end[j] < bed_1$V3[i]) {
        blast_1$q_genes[j] <- base::paste0(blast_1$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",blast_1$q_end[j]- bed_1$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_1$q_start[j] > bed_1$V2[i] & blast_1$q_end[j] < bed_1$V3[i]) {
        blast_1$q_genes[j] <- base::paste0(blast_1$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",blast_1$q_end[j]- blast_1$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_1$s_start[j] < bed_2$V2[i] & blast_1$s_end[j] > bed_2$V3[i]) {
        blast_1$s_genes[j] <- base::paste0(blast_1$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_1$s_start[j] > bed_2$V2[i] & blast_1$s_start[j] < bed_2$V3[i] & blast_1$s_end[j] > bed_2$V3[i]) {
        blast_1$s_genes[j] <- base::paste0(blast_1$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",bed_2$V3[i]-blast_1$s_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_1$s_start[j] < bed_2$V2[i] & bed_2$V2[i] < blast_1$s_end[j] & blast_1$s_end[j] < bed_2$V3[i]) {
        blast_1$s_genes[j] <- base::paste0(blast_1$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_1$s_end[j]- bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_1)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_1$s_start[j] > bed_2$V2[i] & blast_1$s_end[j] < bed_2$V3[i]) {
        blast_1$s_genes[j] <- base::paste0(blast_1$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_1$s_end[j]- blast_1$s_start[j],")",",")
      }
    }
  }
  
  blast_2 <- metablastr::blast_nucleotide_to_nucleotide(query = fasta_1,
                                                        subject = fasta_3,
                                                        db.import = F,
                                                        task="blastn",
                                                        evalue = evalue_cut_off)
  
  blast_2$q_start <- as.numeric(blast_2$q_start)
  blast_2$s_start <- as.numeric(blast_2$s_start)
  
  blast_2$q_end <- as.numeric(blast_2$q_end)
  blast_2$s_end <- as.numeric(blast_2$s_end)
  
  blast_2$q_start1 <- ifelse(blast_2$q_start < blast_2$q_end, blast_2$q_start,blast_2$q_end)
  blast_2$q_end1 <- ifelse(blast_2$q_start < blast_2$q_end, blast_2$q_end,blast_2$q_start)
  
  blast_2$s_start1 <- ifelse(blast_2$s_start < blast_2$s_end, blast_2$s_start, blast_2$s_end)
  blast_2$s_end1 <- ifelse(blast_2$s_start < blast_2$s_end, blast_2$s_end, blast_2$s_start)
  
  
  blast_2$q_start <- blast_2$q_start1
  blast_2$s_start <- blast_2$s_start1
  
  blast_2$q_end <- blast_2$q_end1
  blast_2$s_end <- blast_2$s_end1
  
  blast_2$q_start1 <- NULL
  blast_n$s_start1 <- NULL
  
  blast_2$q_end1 <- NULL
  blast_2$s_end1 <- NULL
  
  blast_2$q_genes <- "genes: "
  blast_2$s_genes <- "genes: "
  
  bed_1 <- bed_1
  bed_3 <- bed_3
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_2$q_start[j] < bed_1$V2[i] & blast_2$q_end[j] > bed_1$V3[i]) {
        blast_2$q_genes[j] <- paste0(blast_2$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_2$q_start[j] > bed_1$V2[i] & blast_2$q_start[j] < bed_1$V3[i] & blast_2$q_end[j] > bed_1$V3[i]) {
        blast_2$q_genes[j] <- base::paste0(blast_2$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",bed_1$V3[i]-blast_2$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_2$q_start[j] < bed_1$V2[i] & bed_1$V2[i] < blast_2$q_end[j] & blast_2$q_end[j] < bed_1$V3[i]) {
        blast_2$q_genes[j] <- base::paste0(blast_2$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",blast_2$q_end[j]- bed_1$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_1)) {
      
      if (blast_2$q_start[j] > bed_1$V2[i] & blast_2$q_end[j] < bed_1$V3[i]) {
        blast_2$q_genes[j] <- base::paste0(blast_2$q_genes[j],bed_1$V4[i]," ","(",bed_1$V3[i]-bed_1$V2[i],")","/","(",blast_2$q_end[j]- blast_2$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_3)) {
      
      if (blast_2$s_start[j] < bed_3$V2[i] & blast_2$s_end[j] > bed_3$V3[i]) {
        blast_2$s_genes[j] <- base::paste0(blast_2$s_genes[j],bed_3$V4[i]," ","(",bed_3$V3[i]-bed_3$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_3)) {
      
      if (blast_2$s_start[j] > bed_3$V2[i] & blast_2$s_start[j] < bed_3$V3[i] & blast_2$s_end[j] > bed_3$V3[i]) {
        blast_2$s_genes[j] <- base::paste0(blast_2$s_genes[j],bed_3$V4[i]," ","(",bed_3$V3[i]-bed_3$V2[i],")","/","(",bed_3$V3[i]-blast_2$s_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_3)) {
      
      if (blast_2$s_start[j] < bed_3$V2[i] & bed_3$V2[i] < blast_2$s_end[j] & blast_2$s_end[j] < bed_3$V3[i]) {
        blast_2$s_genes[j] <- base::paste0(blast_2$s_genes[j],bed_3$V4[i]," ","(",bed_3$V3[i]-bed_3$V2[i],")","/","(",blast_2$s_end[j]- bed_3$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_2)) {
    
    for (i in 1:base::nrow(bed_3)) {
      
      if (blast_2$s_start[j] > bed_3$V2[i] & blast_2$s_end[j] < bed_3$V3[i]) {
        blast_2$s_genes[j] <- base::paste0(blast_2$s_genes[j],bed_3$V4[i]," ","(",bed_3$V3[i]-bed_3$V2[i],")","/","(",blast_2$s_end[j]- blast_2$s_start[j],")",",")
      }
    }
  }
  
  blast_3 <- metablastr::blast_nucleotide_to_nucleotide(query = fasta_1,
                                                        subject = fasta_2,
                                                        db.import = F,
                                                        task="blastn",
                                                        evalue = evalue_cut_off)
  
  blast_3$q_start <- as.numeric(blast_3$q_start)
  blast_3$s_start <- as.numeric(blast_3$s_start)
  
  blast_3$q_end <- as.numeric(blast_3$q_end)
  blast_3$s_end <- as.numeric(blast_3$s_end)
  
  blast_3$q_start1 <- ifelse(blast_3$q_start < blast_3$q_end, blast_3$q_start,blast_3$q_end)
  blast_3$q_end1 <- ifelse(blast_3$q_start < blast_3$q_end, blast_3$q_end,blast_3$q_start)
  
  blast_3$s_start1 <- ifelse(blast_3$s_start < blast_3$s_end, blast_3$s_start, blast_3$s_end)
  blast_3$s_end1 <- ifelse(blast_3$s_start < blast_3$s_end, blast_3$s_end, blast_3$s_start)
  
  
  blast_3$q_start <- blast_3$q_start1
  blast_3$s_start <- blast_3$s_start1
  
  blast_3$q_end <- blast_3$q_end1
  blast_3$s_end <- blast_3$s_end1
  
  blast_3$q_start1 <- NULL
  blast_n$s_start1 <- NULL
  
  blast_3$q_end1 <- NULL
  blast_3$s_end1 <- NULL
  
  blast_3$q_genes <- "genes: "
  blast_3$s_genes <- "genes: "
  
  bed_2 <- bed_2
  bed_2 <- bed_2
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$q_start[j] < bed_2$V2[i] & blast_3$q_end[j] > bed_2$V3[i]) {
        blast_3$q_genes[j] <- paste0(blast_3$q_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$q_start[j] > bed_2$V2[i] & blast_3$q_start[j] < bed_2$V3[i] & blast_3$q_end[j] > bed_2$V3[i]) {
        blast_3$q_genes[j] <- base::paste0(blast_3$q_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",bed_2$V3[i]-blast_3$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$q_start[j] < bed_2$V2[i] & bed_2$V2[i] < blast_3$q_end[j] & blast_3$q_end[j] < bed_2$V3[i]) {
        blast_3$q_genes[j] <- base::paste0(blast_3$q_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_3$q_end[j]- bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$q_start[j] > bed_2$V2[i] & blast_3$q_end[j] < bed_2$V3[i]) {
        blast_3$q_genes[j] <- base::paste0(blast_3$q_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_3$q_end[j]- blast_3$q_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$s_start[j] < bed_2$V2[i] & blast_3$s_end[j] > bed_2$V3[i]) {
        blast_3$s_genes[j] <- base::paste0(blast_3$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$s_start[j] > bed_2$V2[i] & blast_3$s_start[j] < bed_2$V3[i] & blast_3$s_end[j] > bed_2$V3[i]) {
        blast_3$s_genes[j] <- base::paste0(blast_3$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",bed_2$V3[i]-blast_3$s_start[j],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$s_start[j] < bed_2$V2[i] & bed_2$V2[i] < blast_3$s_end[j] & blast_3$s_end[j] < bed_2$V3[i]) {
        blast_3$s_genes[j] <- base::paste0(blast_3$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_3$s_end[j]- bed_2$V2[i],")",",")
      }
    }
  }
  
  for (j in 1:base::nrow(blast_3)) {
    
    for (i in 1:base::nrow(bed_2)) {
      
      if (blast_3$s_start[j] > bed_2$V2[i] & blast_3$s_end[j] < bed_2$V3[i]) {
        blast_3$s_genes[j] <- base::paste0(blast_3$s_genes[j],bed_2$V4[i]," ","(",bed_2$V3[i]-bed_2$V2[i],")","/","(",blast_3$s_end[j]- blast_3$s_start[j],")",",")
      }
    }
  }
  
  blast_end <- rbind(blast_1,blast_2,blast_3)
  base::message(base::paste0("Done!"))
  return(blast_end)
}
