#' Annotate BLASTn alignments and determine transfer direction (MT <-> PT)
#'
#' @param fasta_mt Path to mitochondrial FASTA file (Query).
#' @param fasta_pt Path to plastid FASTA file (Subject).
#' @param bed_mt Mitochondrial BED file (4 columns: chrom, start, end, name).
#' @param bed_pt Plastid BED file (4 columns: chrom, start, end, name).
#' @param evalue_cut_off Maximum e-value threshold (default: 1e-06).
#' @param min_length Minimum alignment length (default: 100).
#' @param min_identity Minimum percent identity (default: 70).
#' @param gene_buffer Min difference in sum of gene_perc (default: 20).
#' @param trans_buffer Min difference in sum of transfer_perc (default: 20).
#' @export

transfer_function <- function(fasta_mt, fasta_pt, bed_mt, bed_pt, 
                              evalue_cut_off = 1e-06, 
                              min_length = 100,
                              min_identity = 70,
                              gene_buffer = 20, 
                              trans_buffer = 20) {
  
  if (missing(fasta_mt) || missing(fasta_pt) || missing(bed_mt) || missing(bed_pt)) {
    stop("Wszystkie pliki wejściowe są wymagane.")
  }

  load_bed_local <- function(bed_input) {
    df <- if (is.character(bed_input)) read.table(bed_input, header = FALSE, stringsAsFactors = FALSE) else as.data.frame(bed_input)
    data.frame(
      chrom = tolower(as.character(df[[1]])),
      start = as.numeric(df[[2]]),
      end   = as.numeric(df[[3]]),
      name  = as.character(df[[4]]),
      stringsAsFactors = FALSE
    )
  }

  b_mt <- load_bed_local(bed_mt)
  b_pt <- load_bed_local(bed_pt)
  
  rna_regex <- "^trn|tRNA|^rrn|rRNA|[0-9]+S_rRNA"

  message("Uruchamianie BLASTn...")
  blast_n <- metablastr::blast_nucleotide_to_nucleotide(
    query = fasta_mt, subject = fasta_pt,
    db.import = FALSE, task = "blastn", evalue = evalue_cut_off
  )
  
  if (nrow(blast_n) == 0) {
    message("Brak dopasowań BLAST.")
    return(NULL)
  }

  blast_n$alig_length <- as.numeric(blast_n$alig_length)
  blast_n$perc_identity <- as.numeric(blast_n$perc_identity)
  
  blast_n <- blast_n[blast_n$alig_length >= min_length & 
                     blast_n$perc_identity >= min_identity, ]
  
  if (nrow(blast_n) == 0) return(NULL)

  blast_n$q_id_low <- tolower(blast_n$query_id)
  blast_n$s_id_low <- tolower(blast_n$subject_id)
  blast_n$q_start_fix <- pmin(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
  blast_n$q_end_fix   <- pmax(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
  blast_n$s_start_fix <- pmin(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))
  blast_n$s_end_fix   <- pmax(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))

  get_unique_overlap <- function(intervals) {
    if (nrow(intervals) == 0) return(0)
    intervals <- intervals[order(intervals$start), ]
    sum_len <- 0
    curr_start <- intervals$start[1]
    curr_end <- intervals$end[1]
    
    if(nrow(intervals) > 1) {
      for(i in 2:nrow(intervals)) {
        if(intervals$start[i] < curr_end) {
          curr_end <- max(curr_end, intervals$end[i])
        } else {
          sum_len <- sum_len + (curr_end - curr_start)
          curr_start <- intervals$start[i]
          curr_end <- intervals$end[i]
        }
      }
    }
    sum_len <- sum_len + (curr_end - curr_start)
    return(sum_len)
  }

  annotate_and_get_metrics <- function(df, bed, prefix) {
    chrom_col   <- if(prefix == "q") "q_id_low" else "s_id_low"
    b_start_col <- if(prefix == "q") "q_start_fix" else "s_start_fix"
    b_end_col   <- if(prefix == "q") "q_end_fix" else "s_end_fix"
    
    lapply(1:nrow(df), function(i) {
      hit <- df[i, ]
      overlaps <- bed[bed$chrom == hit[[chrom_col]] & 
                        bed$start < hit[[b_end_col]] & 
                        bed$end > hit[[b_start_col]], ]
      
      if (nrow(overlaps) == 0) {
        return(list(text = "none", sum_g_perc = 0, sum_t_perc = 0, has_cds = FALSE))
      }
      
      is_rna <- grepl(rna_regex, overlaps$name, ignore.case = TRUE)
      has_cds <- any(!is_rna)
      
      ov_intervals <- data.frame(
        start = pmax(overlaps$start, hit[[b_start_col]]),
        end = pmin(overlaps$end, hit[[b_end_col]])
      )
      unique_ov_len <- get_unique_overlap(ov_intervals)
      
      res_text <- sapply(1:nrow(overlaps), function(j) {
        g_len <- overlaps$end[j] - overlaps$start[j]
        ov_len <- max(0, min(overlaps$end[j], hit[[b_end_col]]) - max(overlaps$start[j], hit[[b_start_col]]))
        paste0(overlaps$name[j], " (ov=", round(ov_len, 0), ", g_perc=", round((ov_len/g_len)*100, 1), "%)")
      })
      
      total_g_perc <- sum(sapply(1:nrow(overlaps), function(j) {
        ov_len <- max(0, min(overlaps$end[j], hit[[b_end_col]]) - max(overlaps$start[j], hit[[b_start_col]]))
        (ov_len / (overlaps$end[j] - overlaps$start[j])) * 100
      }))

      return(list(text = paste(res_text, collapse = "; "), 
                  sum_g_perc = total_g_perc, 
                  sum_t_perc = (unique_ov_len / hit$alig_length) * 100, 
                  has_cds = has_cds))
    })
  }

  message("Analiza genów i kierunku...")
  mt_ann <- annotate_and_get_metrics(blast_n, b_mt, "q")
  pt_ann <- annotate_and_get_metrics(blast_n, b_pt, "s")
  
  blast_n$mt_genes <- sapply(mt_ann, function(x) x$text)
  blast_n$pt_genes <- sapply(pt_ann, function(x) x$text)
  blast_n$direction <- "Unidentified"

  for (i in 1:nrow(blast_n)) {
    m <- mt_ann[[i]]; p <- pt_ann[[i]]
    
    if (!m$has_cds || !p$has_cds) {
      blast_n$direction[i] <- "Unidentified"
      next
    }
    
    g_diff <- abs(m$sum_g_perc - p$sum_g_perc)
    if (g_diff >= gene_buffer) {
      blast_n$direction[i] <- if(m$sum_g_perc > p$sum_g_perc) "MT -> PT" else "PT -> MT"
    } else {
      t_diff <- abs(m$sum_t_perc - p$sum_t_perc)
      if (t_diff >= trans_buffer) {
        blast_n$direction[i] <- if(m$sum_t_perc > p$sum_t_perc) "MT -> PT" else "PT -> MT"
      }
    }
  }

  cols_to_remove <- c("q_start_fix", "q_end_fix", "s_start_fix", "s_end_fix", "q_id_low", "s_id_low")
  blast_n <- blast_n[, !(names(blast_n) %in% cols_to_remove)]
  
  message("Gotowe!")
  return(as.data.frame(blast_n))
}