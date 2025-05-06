#' Annotate BLASTn alignments with overlapping gene features
#'
#' Performs nucleotide BLAST alignment, adjusts coordinates for directionality,
#' and annotates results with overlapping gene features from BED files.
#'
#' @param fasta_q Path to query FASTA file (required)
#' @param fasta_s Path to subject FASTA file (required)
#' @param bed_q query species BED file (required)
#' @param bed_s subject species BED file (required)
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
#'   fasta_q = "query_genome.fna",
#'   fasta_s = "subject_genome.fna",
#'   bed_q = "query_genes.bed",
#'   bed_s = "subject_genes.bed",
#'   evalue_cut_off = 1e-5
#' )
#' }
#' @export

transfer_function <- function(fasta_q="",fasta_s="",bed_q=bed_q, bed_s=bed_s, evalue_cut_off=0.001) {

  if (base::missing(fasta_q)) {
    stop("The fasta_q predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(fasta_s)) {
    stop("The fasta_s predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(bed_q)) {
    stop("The bed_q predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(bed_s)) {
    stop("The bed_s predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  blast_n <- metablastr::blast_nucleotide_to_nucleotide(query = fasta_q,
                                                        subject =fasta_s,
                                                        db.import = F,
                                                        task="blastn",
                                                        evalue = evalue_cut_off)

  blast_n$q_start <- as.numeric(blast_n$q_start)
  blast_n$s_start <- as.numeric(blast_n$s_start)

  blast_n$q_end <- as.numeric(blast_n$q_end)
  blast_n$s_end <- as.numeric(blast_n$s_end)

  blast_n$q_start1 <- ifelse(blast_n$q_start < blast_n$q_end, blast_n$q_start,blast_n$q_end)
  blast_n$q_end1 <- ifelse(blast_n$q_start < blast_n$q_end, blast_n$q_end,blast_n$q_start)

  blast_n$s_start1 <- ifelse(blast_n$s_start < blast_n$s_end, blast_n$s_start, blast_n$s_end)
  blast_n$s_end1 <- ifelse(blast_n$s_start < blast_n$s_end, blast_n$s_end, blast_n$s_start)


  blast_n$q_start <- blast_n$q_start1
  blast_n$s_start <- blast_n$s_start1

  blast_n$q_end <- blast_n$q_end1
  blast_n$s_end <- blast_n$s_end1

  blast_n$q_start1 <- NULL
  blast_n$s_start1 <- NULL

  blast_n$q_end1 <- NULL
  blast_n$s_end1 <- NULL

  blast_n$q_genes <- "genes: "
  blast_n$s_genes <- "genes: "

  bed_q <- bed_q
  bed_s <- bed_s

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_q)) {

      if (blast_n$q_start[j] < bed_q$V2[i] & blast_n$q_end[j] > bed_q$V3[i]) {
        blast_n$q_genes[j] <- paste0(blast_n$q_genes[j],bed_q$V4[i]," ","(",bed_q$V3[i]-bed_q$V2[i],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_q)) {

      if (blast_n$q_start[j] > bed_q$V2[i] & blast_n$q_start[j] < bed_q$V3[i] & blast_n$q_end[j] > bed_q$V3[i]) {
        blast_n$q_genes[j] <- base::paste0(blast_n$q_genes[j],bed_q$V4[i]," ","(",bed_q$V3[i]-bed_q$V2[i],")","/","(",bed_q$V3[i]-blast_n$q_start[j],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_q)) {

      if (blast_n$q_start[j] < bed_q$V2[i] & bed_q$V2[i] < blast_n$q_end[j] & blast_n$q_end[j] < bed_q$V3[i]) {
        blast_n$q_genes[j] <- base::paste0(blast_n$q_genes[j],bed_q$V4[i]," ","(",bed_q$V3[i]-bed_q$V2[i],")","/","(",blast_n$q_end[j]- bed_q$V2[i],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_q)) {

      if (blast_n$q_start[j] > bed_q$V2[i] & blast_n$q_end[j] < bed_q$V3[i]) {
        blast_n$q_genes[j] <- base::paste0(blast_n$q_genes[j],bed_q$V4[i]," ","(",bed_q$V3[i]-bed_q$V2[i],")","/","(",blast_n$q_end[j]- blast_n$q_start[j],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_s)) {

      if (blast_n$s_start[j] < bed_s$V2[i] & blast_n$s_end[j] > bed_s$V3[i]) {
        blast_n$s_genes[j] <- base::paste0(blast_n$s_genes[j],bed_s$V4[i]," ","(",bed_s$V3[i]-bed_s$V2[i],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_s)) {

      if (blast_n$s_start[j] > bed_s$V2[i] & blast_n$s_start[j] < bed_s$V3[i] & blast_n$s_end[j] > bed_s$V3[i]) {
        blast_n$s_genes[j] <- base::paste0(blast_n$s_genes[j],bed_s$V4[i]," ","(",bed_s$V3[i]-bed_s$V2[i],")","/","(",bed_s$V3[i]-blast_n$s_start[j],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_s)) {

      if (blast_n$s_start[j] < bed_s$V2[i] & bed_s$V2[i] < blast_n$s_end[j] & blast_n$s_end[j] < bed_s$V3[i]) {
        blast_n$s_genes[j] <- base::paste0(blast_n$s_genes[j],bed_s$V4[i]," ","(",bed_s$V3[i]-bed_s$V2[i],")","/","(",blast_n$s_end[j]- bed_s$V2[i],")",",")
      }
    }
  }

  for (j in 1:base::nrow(blast_n)) {

    for (i in 1:base::nrow(bed_s)) {

      if (blast_n$s_start[j] > bed_s$V2[i] & blast_n$s_end[j] < bed_s$V3[i]) {
        blast_n$s_genes[j] <- base::paste0(blast_n$s_genes[j],bed_s$V4[i]," ","(",bed_s$V3[i]-bed_s$V2[i],")","/","(",blast_n$s_end[j]- blast_n$s_start[j],")",",")
      }
    }
  }

  base::message(base::paste0("Done!"))
  return(blast_n)
}
