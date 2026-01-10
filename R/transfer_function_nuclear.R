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
    df <- if (is.character(bed_input)) read.table(bed_input, header = FALSE, stringsAsFactors = FALSE) else as.data.frame(bed_input)
    df[[1]] <- tolower(as.character(df[[1]]))
    return(data.frame(chrom = df[[1]], start = as.numeric(df[[2]]), 
                      end = as.numeric(df[[3]]), name = as.character(df[[4]]),
                      stringsAsFactors = FALSE))
  }
  
  b_mt  <- load_bed(bed_mt); b_pt  <- load_bed(bed_pt); b_nuc <- load_bed(bed_nuc)
  
  trna_regex <- "^trn|tRNA"
  rrna_regex <- "^rrn|rRNA|[0-9]+S_rRNA"
  
  run_single_transfer <- function(q_fasta, s_fasta, q_bed, s_bed, q_label, s_label, 
                                  validation_blast = NULL, validation_label = "3rd genome",
                                  q_label_short = NULL, s_label_short = NULL) {
    
    if (is.null(q_label_short)) q_label_short <- q_label
    if (is.null(s_label_short)) s_label_short <- s_label
    
    message(sprintf("Analysing: %s vs %s...", q_label, s_label))
    
    blast_n <- metablastr::blast_nucleotide_to_nucleotide(
      query = q_fasta, subject = s_fasta, db.import = FALSE, task = "blastn", evalue = evalue_cut_off
    )
    
    if (is.null(blast_n) || nrow(blast_n) == 0) return(NULL)
    
    blast_n$alig_length <- as.numeric(blast_n$alig_length)
    blast_n$perc_identity <- as.numeric(blast_n$perc_identity)
    blast_n$bit_score <- as.numeric(blast_n$bit_score)
    blast_n$q_start <- as.numeric(blast_n$q_start)
    blast_n$q_end <- as.numeric(blast_n$q_end)
    blast_n$s_start <- as.numeric(blast_n$s_start)
    blast_n$s_end <- as.numeric(blast_n$s_end)
    
    blast_n <- blast_n[blast_n$alig_length >= min_length & blast_n$perc_identity >= min_identity, ]
    if (nrow(blast_n) == 0) return(NULL)
    
    blast_n$query_id_lower <- tolower(blast_n$query_id)
    blast_n$subject_id_lower <- tolower(blast_n$subject_id)
    blast_n$q_start_fix <- pmin(blast_n$q_start, blast_n$q_end)
    blast_n$q_end_fix   <- pmax(blast_n$q_start, blast_n$q_end)
    blast_n$s_start_fix <- pmin(blast_n$s_start, blast_n$s_end)
    blast_n$s_end_fix   <- pmax(blast_n$s_start, blast_n$s_end)
    
    annotate_and_get_metrics <- function(df, bed, prefix) {
      chrom_col   <- if(prefix == "q") "query_id_lower" else "subject_id_lower"
      b_start_col <- if(prefix == "q") "q_start_fix" else "s_start_fix"
      b_end_col   <- if(prefix == "q") "q_end_fix" else "s_end_fix"
      
      lapply(1:nrow(df), function(i) {
        hit <- df[i, ]
        overlaps <- bed[bed$chrom == hit[[chrom_col]] & bed$start < hit[[b_end_col]] & bed$end > hit[[b_start_col]], ]
        
        if (nrow(overlaps) == 0) {
          return(list(text = NA, sum_g_perc = NA, sum_t_perc = NA, 
                      only_trna = FALSE, only_rrna = FALSE, only_rna = FALSE, 
                      has_genes = FALSE, has_coding = FALSE))
        }
        
        is_trna <- grepl(trna_regex, overlaps$name, ignore.case = TRUE)
        is_rrna <- grepl(rrna_regex, overlaps$name, ignore.case = TRUE)
        is_any_rna <- is_trna | is_rrna
        
        coding_genes <- overlaps[!is_any_rna, ]
        rna_genes <- overlaps[is_any_rna, ]
        
        target_ov <- if (nrow(coding_genes) > 0) coding_genes else rna_genes
        
        gene_results <- list()
        for(j in 1:nrow(target_ov)) {
          g_len  <- target_ov$end[j] - target_ov$start[j]
          ov_len <- max(0, min(target_ov$end[j], hit[[b_end_col]]) - max(target_ov$start[j], hit[[b_start_col]]))
          g_perc <- (ov_len / g_len) * 100
          t_perc <- (ov_len / hit$alig_length) * 100
          gene_results[[j]] <- list(g_perc = g_perc, t_perc = t_perc, g_len = g_len, ov_len = ov_len, name = target_ov$name[j])
        }
        
        sum_g_perc <- sum(sapply(gene_results, function(x) x$g_perc))
        sum_t_perc <- sum(sapply(gene_results, function(x) x$t_perc))
        res_text <- sapply(gene_results, function(x) {
          paste0(x$name, " (g_len=", x$g_len, ", ov_len=", round(x$ov_len, 0), 
                 ", g_perc=", round(x$g_perc, 1), "%, t_perc=", round(x$t_perc, 1), "%)")
        })
        
        return(list(text = paste(res_text, collapse = "; "), sum_g_perc = sum_g_perc, sum_t_perc = sum_t_perc, 
                    only_trna = all(is_trna), only_rrna = all(is_rrna), only_rna = all(is_any_rna), 
                    has_genes = TRUE, has_coding = nrow(coding_genes) > 0))
      })
    }
    
    q_ann <- annotate_and_get_metrics(blast_n, q_bed, "q")
    s_ann <- annotate_and_get_metrics(blast_n, s_bed, "s")
    
    blast_n[[paste0(tolower(q_label_short), "_genes")]] <- sapply(q_ann, function(x) x$text)
    blast_n[[paste0(tolower(s_label_short), "_genes")]] <- sapply(s_ann, function(x) x$text)
    
    blast_n$direction <- "Undefined"
    blast_n$excluded_coords <- NA_character_
    
    for (i in 1:nrow(blast_n)) {
      m <- q_ann[[i]]
      p <- s_ann[[i]]
      
      if (!is.null(validation_blast)) {
        cross <- validation_blast[
          validation_blast$query_id == blast_n$query_id[i] & 
          abs(validation_blast$q_start - blast_n$q_start[i]) < 50, 
        ]
        
        if (nrow(cross) > 0 && max(cross$bit_score) > blast_n$bit_score[i]) {
          best_cross <- cross[which.max(cross$bit_score), ]
          blast_n$direction[i] <- paste0("Excluded: Stronger match in ", validation_label)
          blast_n$excluded_coords[i] <- paste0(
            validation_label, " [", 
            as.numeric(best_cross$s_start), "-", 
            as.numeric(best_cross$s_end), 
            "] score=", 
            round(best_cross$bit_score, 1)
          )
          next
        }
      }
      
      q_has_genes <- !is.na(m$sum_g_perc) && m$has_genes
      s_has_genes <- !is.na(p$sum_g_perc) && p$has_genes
      q_has_coding <- isTRUE(m$has_coding)
      s_has_coding <- isTRUE(p$has_coding)
      
      if (!q_has_genes && s_has_genes) {
        blast_n$direction[i] <- paste(q_label, "->", s_label)
        next
      } else if (q_has_genes && !s_has_genes) {
        blast_n$direction[i] <- paste(s_label, "->", q_label)
        next
      } else if (!q_has_genes && !s_has_genes) {
        blast_n$direction[i] <- "Undefined (No genes)"
        next
      }
      
      if (m$only_rna && !p$only_rna) {
        blast_n$direction[i] <- paste(s_label, "->", q_label)
        next
      } else if (!m$only_rna && p$only_rna) {
        blast_n$direction[i] <- paste(q_label, "->", s_label)
        next
      }
      
      if (m$only_rna && p$only_rna) {
        if (m$only_trna && p$only_trna) {
          blast_n$direction[i] <- "Undefined (tRNA only)"
        } else if (m$only_rrna && p$only_rrna) {
          blast_n$direction[i] <- "Undefined (rRNA only)"
        } else {
          blast_n$direction[i] <- "Undefined (Mixed RNA)"
        }
        next
      }
      
      g_diff <- abs(m$sum_g_perc - p$sum_g_perc)
      t_diff <- abs(m$sum_t_perc - p$sum_t_perc)
      
      if (g_diff >= gene_buffer) {
        blast_n$direction[i] <- if(m$sum_g_perc > p$sum_g_perc) {
          paste(q_label, "->", s_label)
        } else {
          paste(s_label, "->", q_label)
        }
      } else if (t_diff >= trans_buffer) {
        blast_n$direction[i] <- if(m$sum_t_perc > p$sum_t_perc) {
          paste(q_label, "->", s_label)
        } else {
          paste(s_label, "->", q_label)
        }
      }
    }
    
    cols_to_remove <- c("q_start_fix", "q_end_fix", "s_start_fix", "s_end_fix", 
                        "query_id_lower", "subject_id_lower")
    return(blast_n[, !(names(blast_n) %in% cols_to_remove)])
  }
  
  raw_mt_pt  <- run_single_transfer(fasta_mt, fasta_pt, b_mt, b_pt, "MT", "PT",
                                    q_label_short = "mt", s_label_short = "pt")
  raw_mt_nuc <- run_single_transfer(fasta_mt, fasta_nuc, b_mt, b_nuc, "MT", "NUC",
                                    q_label_short = "mt", s_label_short = "nuc")
  raw_pt_nuc <- run_single_transfer(fasta_pt, fasta_nuc, b_pt, b_nuc, "PT", "NUC",
                                    q_label_short = "pt", s_label_short = "nuc")
  
  res_mt_pt  <- run_single_transfer(fasta_mt, fasta_pt, b_mt, b_pt, "MT", "PT", 
                                    raw_mt_nuc, "NUC",
                                    q_label_short = "mt", s_label_short = "pt")
  res_mt_nuc <- run_single_transfer(fasta_mt, fasta_nuc, b_mt, b_nuc, "MT", "NUC", 
                                    raw_mt_pt, "PT",
                                    q_label_short = "mt", s_label_short = "nuc")
  res_pt_nuc <- run_single_transfer(fasta_pt, fasta_nuc, b_pt, b_nuc, "PT", "NUC", 
                                    raw_mt_pt, "MT",
                                    q_label_short = "pt", s_label_short = "nuc")
  
  final_res <- data.frame()
  combined  <- list(res_mt_pt, res_mt_nuc, res_pt_nuc)
  
  for (res in combined) {
    if (!is.null(res)) {
      for(col in c("mt_genes", "pt_genes", "nuc_genes")) {
        if(!(col %in% names(res))) res[[col]] <- NA
      }
      if(!("excluded_coords" %in% names(res))) res$excluded_coords <- NA
      
      if (nrow(final_res) == 0) {
        final_res <- res
      } else {
        final_res <- rbind(final_res, res)
      }
    }
  }
  
  desired_order <- c(
    "query_id", "subject_id",
    "perc_identity", "num_ident_matches", "alig_length", "mismatches", 
    "gap_openings", "n_gaps", "pos_match", "ppos", 
    "q_start", "q_end", "q_len", "qcov", "qcovhsp",
    "s_start", "s_end", "s_len",
    "evalue", "bit_score", "score_raw",
    "mt_genes", "pt_genes", "nuc_genes",
    "direction", "excluded_coords"
  )
  
  for(col in desired_order) {
    if(!(col %in% names(final_res))) {
      final_res[[col]] <- NA
    }
  }
  
  final_res <- final_res[, desired_order]
  
  return(as.data.frame(final_res))
}
