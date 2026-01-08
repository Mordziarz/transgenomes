#' Intergenomic DNA Transfer Analysis (MT, PT, NUC)
#'
#' This function identifies and annotates DNA sequence transfers between mitochondrial (MT), 
#' plastid (PT), and nuclear (NUC) genomes. By aggregating metrics from all overlapping 
#' genomic features and applying strict quality filters (length and identity), it 
#' determines the most likely direction of transfer based on gene completeness 
#' and alignment dominance.
#'
#' @param fasta_mt Path to the mitochondrial genome FASTA file.
#' @param fasta_pt Path to the plastid genome FASTA file.
#' @param fasta_nuc Path to the nuclear genome FASTA file.
#' @param bed_mt Mitochondrial BED file (min. 4 columns: chrom, start, end, name).
#' @param bed_pt Plastid BED file (min. 4 columns).
#' @param bed_nuc Nuclear BED file (min. 4 columns).
#' @param evalue_cut_off E-value threshold for BLASTn (default: 1e-06).
#' @param min_length Minimum BLAST alignment length (alig_length) in bp (default: 100).
#' @param min_identity Minimum percent identity (perc_identity) (default: 70).
#' @param gene_buffer Threshold for the difference in total gene completeness (%) to assign direction (default: 20).
#' @param trans_buffer Threshold for the difference in total transfer dominance (%) to assign direction (default: 20).
#'
#' @return A data.frame containing BLASTn results merged with specific HGT metrics and predicted transfer directions.
#' @export

transfer_function_nuclear <- function(fasta_mt, fasta_pt, fasta_nuc, 
                                      bed_mt, bed_pt, bed_nuc,
                                      evalue_cut_off = 1e-06, 
                                      min_length = 100,
                                      min_identity = 70,
                                      gene_buffer = 20, 
                                      trans_buffer = 20) {
  
  load_bed <- function(bed_input) {
    if (is.character(bed_input)) {
      df <- read.table(bed_input, header = FALSE, stringsAsFactors = FALSE)
    } else {
      df <- as.data.frame(bed_input)
    }
    df[[1]] <- tolower(as.character(df[[1]]))
    res <- data.frame(
      chrom = df[[1]],
      start = as.numeric(df[[2]]),
      end   = as.numeric(df[[3]]),
      name  = as.character(df[[4]]),
      stringsAsFactors = FALSE
    )
    return(res)
  }
  
  b_mt  <- load_bed(bed_mt)
  b_pt  <- load_bed(bed_pt)
  b_nuc <- load_bed(bed_nuc)
  
  trna_regex <- "^trn|tRNA"
  
  run_single_transfer <- function(q_fasta, s_fasta, q_bed, s_bed, q_label, s_label) {
    message(sprintf("Running BLASTn: %s vs %s...", q_label, s_label))
    
    blast_n <- metablastr::blast_nucleotide_to_nucleotide(
      query = q_fasta, subject = s_fasta,
      db.import = FALSE, task = "blastn", evalue = evalue_cut_off
    )
    
    if (nrow(blast_n) == 0) return(NULL)
    
    blast_n$alig_length <- as.numeric(blast_n$alig_length)
    blast_n$perc_identity <- as.numeric(blast_n$perc_identity)
    
    blast_n <- blast_n[blast_n$alig_length >= min_length & 
                         blast_n$perc_identity >= min_identity, ]
    
    if (nrow(blast_n) == 0) {
      message(sprintf("No hits meeting quality criteria (len>=%d, id>=%d) for %s vs %s.", 
                      min_length, min_identity, q_label, s_label))
      return(NULL)
    }
    
    blast_n$query_id_lower <- tolower(blast_n$query_id)
    blast_n$subject_id_lower <- tolower(blast_n$subject_id)
    
    blast_n$q_start_fix <- pmin(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
    blast_n$q_end_fix   <- pmax(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
    blast_n$s_start_fix <- pmin(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))
    blast_n$s_end_fix   <- pmax(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))
    
    annotate_and_get_metrics <- function(df, bed, prefix) {
      chrom_col   <- if(prefix == "q") "query_id_lower" else "subject_id_lower"
      b_start_col <- if(prefix == "q") "q_start_fix" else "s_start_fix"
      b_end_col   <- if(prefix == "q") "q_end_fix" else "s_end_fix"
      
      lapply(1:nrow(df), function(i) {
        hit <- df[i, ]
        overlaps <- bed[bed$chrom == hit[[chrom_col]] & 
                          bed$start < hit[[b_end_col]] & 
                          bed$end > hit[[b_start_col]], ]
        
        if (nrow(overlaps) == 0) {
          return(list(text = "none", sum_g_perc = 0, sum_t_perc = 0, only_trna = FALSE))
        }
        
        is_trna <- grepl(trna_regex, overlaps$name, ignore.case = TRUE)
        only_trna <- all(is_trna)
        
        gene_results <- list()
        for(j in 1:nrow(overlaps)) {
          g_len  <- overlaps$end[j] - overlaps$start[j]
          ov_len <- max(0, min(overlaps$end[j], hit[[b_end_col]]) - max(overlaps$start[j], hit[[b_start_col]]))
          
          g_perc <- (ov_len / g_len) * 100
          t_perc <- (ov_len / hit$alig_length) * 100
          
          gene_results[[j]] <- list(g_perc = g_perc, t_perc = t_perc, g_len = g_len, 
                                    ov_len = ov_len, name = overlaps$name[j])
        }
        
        sum_g_perc <- sum(sapply(gene_results, function(x) x$g_perc))
        sum_t_perc <- sum(sapply(gene_results, function(x) x$t_perc))
        
        res_text <- sapply(gene_results, function(x) {
          paste0(x$name, " (g_len=", x$g_len, 
                 ", ov_len=", round(x$ov_len, 0), 
                 ", g_perc=", round(x$g_perc, 1), "%, ",
                 "t_perc=", round(x$t_perc, 1), "%)")
        })
        
        return(list(text = paste(res_text, collapse = "; "), 
                    sum_g_perc = sum_g_perc, 
                    sum_t_perc = sum_t_perc, 
                    only_trna = only_trna))
      })
    }
    
    q_ann <- annotate_and_get_metrics(blast_n, q_bed, "q")
    s_ann <- annotate_and_get_metrics(blast_n, s_bed, "s")
    
    blast_n[[paste0(tolower(q_label), "_genes")]] <- sapply(q_ann, function(x) x$text)
    blast_n[[paste0(tolower(s_label), "_genes")]] <- sapply(s_ann, function(x) x$text)
    
    blast_n$direction <- "unknown"
    
    for (i in 1:nrow(blast_n)) {
      m <- q_ann[[i]]; p <- s_ann[[i]]
      if (m$only_trna && p$only_trna) {
        blast_n$direction[i] <- "unknown"
        next
      }
      
      g_diff <- abs(m$sum_g_perc - p$sum_g_perc)
      if (g_diff >= gene_buffer) {
        blast_n$direction[i] <- if(m$sum_g_perc > p$sum_g_perc) paste(q_label, "->", s_label) else paste(s_label, "->", q_label)
      } else {
        t_diff <- abs(m$sum_t_perc - p$sum_t_perc)
        if (t_diff >= trans_buffer) {
          blast_n$direction[i] <- if(m$sum_t_perc > p$sum_t_perc) paste(q_label, "->", s_label) else paste(s_label, "->", q_label)
        }
      }
    }
    
    cols_to_remove <- c("q_start_fix", "q_end_fix", "s_start_fix", "s_end_fix", "query_id_lower", "subject_id_lower")
    return(blast_n[, !(names(blast_n) %in% cols_to_remove)])
  }
  
  res_mt_pt  <- run_single_transfer(fasta_mt, fasta_pt, b_mt, b_pt, "MT", "PT")
  res_mt_nuc <- run_single_transfer(fasta_mt, fasta_nuc, b_mt, b_nuc, "MT", "NUC")
  res_pt_nuc <- run_single_transfer(fasta_pt, fasta_nuc, b_pt, b_nuc, "PT", "NUC")
  
  final_res <- data.frame()
  combined <- list(res_mt_pt, res_mt_nuc, res_pt_nuc)
  for (res in combined) {
    if (!is.null(res)) {
      if (nrow(final_res) == 0) final_res <- res else final_res <- merge(final_res, res, all = TRUE)
    }
  }
  
  return(as.data.frame(final_res))
}