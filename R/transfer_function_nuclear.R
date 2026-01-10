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
                                      trans_buffer = 20,
                                      nuclear_ratio_threshold = 1.1,
                                      nuclear_min_bit_score = 100) {
  
  if (!is.numeric(nuclear_ratio_threshold) || nuclear_ratio_threshold < 0 || nuclear_ratio_threshold > 5) {
    stop("nuclear_ratio_threshold must be a numeric value between 0 and 5")
  }
  if (!is.numeric(nuclear_min_bit_score) || nuclear_min_bit_score <= 0) {
    stop("nuclear_min_bit_score must be a positive numeric value")
  }
  
  load_bed <- function(bed_input) {
    df <- if (is.character(bed_input)) read.table(bed_input, header = FALSE, stringsAsFactors = FALSE) else as.data.frame(bed_input)
    df[[1]] <- tolower(as.character(df[[1]]))
    return(data.frame(chrom = df[[1]], start = as.numeric(df[[2]]), 
                      end = as.numeric(df[[3]]), name = as.character(df[[4]]),
                      stringsAsFactors = FALSE))
  }
  
  b_mt  <- load_bed(bed_mt); b_pt  <- load_bed(bed_pt); b_nuc <- load_bed(bed_nuc)
  trna_regex <- "^trn|tRNA"; rrna_regex <- "^rrn|rRNA|[0-9]+S_rRNA"
  
  run_single_transfer <- function(q_fasta, s_fasta, q_bed, s_bed, q_label, s_label, 
                                  validation_blast = NULL, validation_label = "3rd genome",
                                  q_label_short = NULL, s_label_short = NULL,
                                  nuc_validation_blast = NULL) {
    
    if (is.null(q_label_short)) q_label_short <- q_label
    if (is.null(s_label_short)) s_label_short <- s_label
    message(sprintf("Analysing: %s vs %s...", q_label, s_label))
    
    blast_n <- metablastr::blast_nucleotide_to_nucleotide(
      query = q_fasta, subject = s_fasta, db.import = FALSE, task = "blastn", evalue = evalue_cut_off)
    if (is.null(blast_n) || nrow(blast_n) == 0) return(NULL)
    
    cols_to_num <- c("alig_length", "perc_identity", "bit_score", "q_start", "q_end", "s_start", "s_end")
    blast_n[cols_to_num] <- lapply(blast_n[cols_to_num], as.numeric)
    blast_n <- blast_n[blast_n$alig_length >= min_length & blast_n$perc_identity >= min_identity, ]
    if (nrow(blast_n) == 0) return(NULL)
    
    blast_n$q_start_fix <- pmin(blast_n$q_start, blast_n$q_end)
    blast_n$q_end_fix   <- pmax(blast_n$q_start, blast_n$q_end)
    blast_n$s_start_fix <- pmin(blast_n$s_start, blast_n$s_end)
    blast_n$s_end_fix   <- pmax(blast_n$s_start, blast_n$s_end)
    
    annotate_metrics <- function(df, bed, prefix) {
      id_col <- if(prefix == "q") "query_id" else "subject_id"
      s_col  <- if(prefix == "q") "q_start_fix" else "s_start_fix"
      e_col  <- if(prefix == "q") "q_end_fix" else "s_end_fix"
      lapply(1:nrow(df), function(i) {
        hit <- df[i, ]
        ov  <- bed[bed$chrom == tolower(hit[[id_col]]) & bed$start < hit[[e_col]] & bed$end > hit[[s_col]], ]
        if (nrow(ov) == 0) return(list(text = NA, sum_g = 0, sum_t = 0, has_genes = FALSE, 
                                       only_rna = FALSE, only_trna = FALSE, only_rrna = FALSE))
        is_trna <- grepl(trna_regex, ov$name, ignore.case = TRUE)
        is_rrna <- grepl(rrna_regex, ov$name, ignore.case = TRUE)
        is_rna  <- is_trna | is_rrna
        coding_genes <- ov[!is_rna, ]
        target_ov <- if (nrow(coding_genes) > 0) coding_genes else ov
        gene_list <- lapply(1:nrow(target_ov), function(j) {
          g_len  <- target_ov$end[j] - target_ov$start[j]
          ov_len <- max(0, min(target_ov$end[j], hit[[e_col]]) - max(target_ov$start[j], hit[[s_col]]))
          c(g_perc = (ov_len / g_len) * 100, t_perc = (ov_len / hit$alig_length) * 100, 
            ov_len = ov_len, g_len = g_len)
        })
        gene_mat <- do.call(rbind, gene_list)
        res_text <- paste(sapply(1:nrow(target_ov), function(j) {
          paste0(target_ov$name[j], " (g_len=", gene_mat[j, "g_len"],
                 ", ov_len=", round(gene_mat[j, "ov_len"], 0),
                 ", g_perc=", round(gene_mat[j, "g_perc"], 1), 
                 "%, t_perc=", round(gene_mat[j, "t_perc"], 1), "%)")
        }), collapse = "; ")
        only_trna <- all(is_trna[is_rna]) && any(is_trna)
        only_rrna <- all(is_rrna[is_rna]) && any(is_rrna)
        only_rna_all <- all(is_rna) && any(is_rna)
        return(list(text = res_text, sum_g = sum(gene_mat[,"g_perc"]), sum_t = sum(gene_mat[,"t_perc"]),
                    has_genes = TRUE, only_rna = only_rna_all, only_trna = only_trna, only_rrna = only_rrna))
      })
    }
    
    q_ann <- annotate_metrics(blast_n, q_bed, "q")
    s_ann <- annotate_metrics(blast_n, s_bed, "s")
    blast_n[[paste0(tolower(q_label_short), "_genes")]] <- sapply(q_ann, function(x) x$text)
    blast_n[[paste0(tolower(s_label_short), "_genes")]] <- sapply(s_ann, function(x) x$text)
    blast_n$direction <- "Undefined"
    blast_n$excluded_coords <- NA_character_
    blast_n$nuclear_paralog_flag <- NA_character_
    
    for (i in 1:nrow(blast_n)) {
      m <- q_ann[[i]]; p <- s_ann[[i]]
      
      if (!is.null(nuc_validation_blast)) {
        tryCatch({
          relevant_nuc_query <- nuc_validation_blast[nuc_validation_blast$query_id == blast_n$query_id[i], ]
          relevant_nuc_subject <- nuc_validation_blast[nuc_validation_blast$query_id == blast_n$subject_id[i], ]
          
          if (nrow(relevant_nuc_query) > 0 && nrow(relevant_nuc_subject) > 0) {
            blast_n$nuclear_paralog_flag[i] <- "BOTH_ENDPOINTS_NUCLEAR"
            blast_n$direction[i] <- "Undefined (Both endpoints in nuclear)"
            next
          }
          
          if (nrow(relevant_nuc_query) > 0) {
            nuc_q_start <- pmin(as.numeric(relevant_nuc_query$q_start), as.numeric(relevant_nuc_query$q_end))
            nuc_q_end   <- pmax(as.numeric(relevant_nuc_query$q_start), as.numeric(relevant_nuc_query$q_end))
            nuc_cross <- relevant_nuc_query[nuc_q_start < blast_n$q_end_fix[i] & nuc_q_end > blast_n$q_start_fix[i], ]
            
            if (nrow(nuc_cross) > 0) {
              best_nuc_idx <- which.max(as.numeric(nuc_cross$bit_score))
              best_nuc <- nuc_cross[best_nuc_idx, ]
              current_bit_score <- as.numeric(blast_n$bit_score[i])
              current_alig_length <- as.numeric(blast_n$alig_length[i])
              nuc_bit_score <- as.numeric(best_nuc$bit_score)
              nuc_alig_length <- as.numeric(best_nuc$alig_length)
              
              if (!is.na(nuc_bit_score) && !is.na(current_bit_score) && 
                  !is.na(nuc_alig_length) && !is.na(current_alig_length) &&
                  current_alig_length > 0 && nuc_alig_length > 0) {
                nuc_density <- nuc_bit_score / nuc_alig_length
                cur_density <- current_bit_score / current_alig_length
                density_ratio <- nuc_density / cur_density
                if (nuc_density >= (nuclear_min_bit_score / 1000) && 
                    density_ratio >= nuclear_ratio_threshold) {
                  blast_n$nuclear_paralog_flag[i] <- paste0(
                    "QUERY_NUCLEAR_RISK [NUC_density=", round(nuc_density, 2),
                    " vs ", tolower(s_label), "_density=", round(cur_density, 2),
                    " ratio=", round(density_ratio, 2), "]")
                  blast_n$direction[i] <- "Undefined (Query nuclear paralog)"
                  next
                }
              }
            }
          }
          
          if (nrow(relevant_nuc_subject) > 0 && is.na(blast_n$nuclear_paralog_flag[i])) {
            nuc_q_start <- pmin(as.numeric(relevant_nuc_subject$q_start), as.numeric(relevant_nuc_subject$q_end))
            nuc_q_end   <- pmax(as.numeric(relevant_nuc_subject$q_start), as.numeric(relevant_nuc_subject$q_end))
            nuc_cross <- relevant_nuc_subject[nuc_q_start < blast_n$s_end_fix[i] & nuc_q_end > blast_n$s_start_fix[i], ]
            
            if (nrow(nuc_cross) > 0) {
              best_nuc_idx <- which.max(as.numeric(nuc_cross$bit_score))
              best_nuc <- nuc_cross[best_nuc_idx, ]
              current_bit_score <- as.numeric(blast_n$bit_score[i])
              current_alig_length <- as.numeric(blast_n$alig_length[i])
              nuc_bit_score <- as.numeric(best_nuc$bit_score)
              nuc_alig_length <- as.numeric(best_nuc$alig_length)
              
              if (!is.na(nuc_bit_score) && !is.na(current_bit_score) && 
                  !is.na(nuc_alig_length) && !is.na(current_alig_length) &&
                  current_alig_length > 0 && nuc_alig_length > 0) {
                nuc_density <- nuc_bit_score / nuc_alig_length
                cur_density <- current_bit_score / current_alig_length
                density_ratio <- nuc_density / cur_density
                if (nuc_density >= (nuclear_min_bit_score / 1000) && 
                    density_ratio >= nuclear_ratio_threshold) {
                  blast_n$nuclear_paralog_flag[i] <- paste0(
                    "SUBJECT_NUCLEAR_RISK [NUC_density=", round(nuc_density, 2),
                    " vs ", tolower(q_label), "_density=", round(cur_density, 2),
                    " ratio=", round(density_ratio, 2), "]")
                  blast_n$direction[i] <- "Undefined (Subject nuclear paralog)"
                  next
                }
              }
            }
          }
        }, error = function(e) {
          warning(sprintf("Nuclear validation error at row %d: %s", i, e$message))
        })
      }
      
      if (!is.null(validation_blast)) {
        tryCatch({
          v_q_start <- pmin(as.numeric(validation_blast$q_start), as.numeric(validation_blast$q_end))
          v_q_end   <- pmax(as.numeric(validation_blast$q_start), as.numeric(validation_blast$q_end))
          cross <- validation_blast[validation_blast$query_id == blast_n$query_id[i] & 
                                     v_q_start < blast_n$q_end_fix[i] & v_q_end > blast_n$q_start_fix[i], ]
          if (nrow(cross) > 0) {
            best_c_idx <- which.max(as.numeric(cross$bit_score))
            best_c <- cross[best_c_idx, ]
            best_c_score <- as.numeric(best_c$bit_score)
            if (!is.na(best_c_score) && best_c_score > blast_n$bit_score[i]) {
              blast_n$direction[i] <- paste0("Excluded: Stronger match in ", validation_label)
              blast_n$excluded_coords[i] <- paste0(validation_label, " [", best_c$s_start, "-", 
                                                    best_c$s_end, "] score=", round(best_c_score, 1))
              next
            }
          }
        }, error = function(e) {
          warning(sprintf("Cross-validation error at row %d: %s", i, e$message))
        })
      }
      
      m_has_genes <- m$has_genes; p_has_genes <- p$has_genes
      m_only_rna <- m$only_rna; p_only_rna <- p$only_rna
      m_only_trna <- m$only_trna; p_only_trna <- p$only_trna
      m_only_rrna <- m$only_rrna; p_only_rrna <- p$only_rrna
      
      if ((m_only_rna && !p_has_genes) || (!m_has_genes && p_only_rna)) {
        blast_n$direction[i] <- "Undefined (RNA vs NA)"; next
      }
      if (m_has_genes && !p_has_genes) {
        blast_n$direction[i] <- paste(q_label, "->", s_label); next
      } else if (!m_has_genes && p_has_genes) {
        blast_n$direction[i] <- paste(s_label, "->", q_label); next
      } else if (!m_has_genes && !p_has_genes) {
        blast_n$direction[i] <- "Undefined (No genes)"; next
      }
      
      if (m_only_rna && p_only_rna) {
        if (m_only_rrna && p_only_rrna) { blast_n$direction[i] <- "Undefined (rRNA vs rRNA)"; next }
        if (m_only_trna && p_only_trna) { blast_n$direction[i] <- "Undefined (tRNA vs tRNA)"; next }
        if ((m_only_trna && p_only_rrna) || (m_only_rrna && p_only_trna)) { 
          blast_n$direction[i] <- "Undefined (tRNA vs rRNA)"; next }
        blast_n$direction[i] <- "Undefined (Mixed RNA)"; next
      }
      
      if (!m_only_rna && p_only_rna) {
        blast_n$direction[i] <- paste(q_label, "->", s_label); next
      } else if (m_only_rna && !p_only_rna) {
        blast_n$direction[i] <- paste(s_label, "->", q_label); next
      }
      
      if (abs(m$sum_g - p$sum_g) >= gene_buffer) {
        blast_n$direction[i] <- if (m$sum_g > p$sum_g) paste(q_label, "->", s_label) else paste(s_label, "->", q_label)
      } else if (abs(m$sum_t - p$sum_t) >= trans_buffer) {
        blast_n$direction[i] <- if (m$sum_t > p$sum_t) paste(q_label, "->", s_label) else paste(s_label, "->", q_label)
      }
    }
    
    cols_to_del <- c("q_start_fix", "q_end_fix", "s_start_fix", "s_end_fix")
    return(blast_n[, !(names(blast_n) %in% cols_to_del)])
  }
  
  message("\n=== STAGE 1: Raw pairwise comparisons ===")
  r_mt_pt  <- run_single_transfer(fasta_mt, fasta_pt, b_mt, b_pt, "MT", "PT", q_label_short="mt", s_label_short="pt")
  r_mt_nuc <- run_single_transfer(fasta_mt, fasta_nuc, b_mt, b_nuc, "MT", "NUC", q_label_short="mt", s_label_short="nuc")
  r_pt_nuc <- run_single_transfer(fasta_pt, fasta_nuc, b_pt, b_nuc, "PT", "NUC", q_label_short="pt", s_label_short="nuc")
  
  message("\n=== STAGE 2: Validated transfers with DUAL nuclear paralog filtering ===")
  res_mt_pt  <- run_single_transfer(fasta_mt, fasta_pt, b_mt, b_pt, "MT", "PT", 
                                    r_mt_nuc, "NUC", "mt", "pt", r_mt_nuc)
  res_mt_nuc <- run_single_transfer(fasta_mt, fasta_nuc, b_mt, b_nuc, "MT", "NUC", 
                                    r_mt_pt, "PT", "mt", "nuc", r_pt_nuc)
  res_pt_nuc <- run_single_transfer(fasta_pt, fasta_nuc, b_pt, b_nuc, "PT", "NUC", 
                                    r_mt_pt, "MT", "pt", "nuc", r_mt_nuc)
  
  combined <- list(res_mt_pt, res_mt_nuc, res_pt_nuc)
  final_df <- data.frame()
  for (res in combined) {
    if (!is.null(res)) {
      for(col in c("mt_genes", "pt_genes", "nuc_genes", "excluded_coords", "nuclear_paralog_flag")) {
        if(!(col %in% names(res))) res[[col]] <- NA
      }
      final_df <- rbind(final_df, res)
    }
  }
  
  order_cols <- c("query_id", "subject_id", "perc_identity", "alig_length",
                  "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score",
                  "mt_genes", "pt_genes", "nuc_genes", "direction", "excluded_coords", "nuclear_paralog_flag")
  
  for(c in order_cols) if(!(c %in% names(final_df))) final_df[[c]] <- NA
  
  message("\n=== RESULTS SUMMARY ===")
  message(sprintf("Total hits processed: %d", nrow(final_df)))
  message(sprintf("High confidence transfers: %d", sum(is.na(final_df$nuclear_paralog_flag) & is.na(final_df$excluded_coords))))
  message(sprintf("Nuclear paralog risks detected: %d", sum(!is.na(final_df$nuclear_paralog_flag))))
  message(sprintf("Hits excluded by cross-validation: %d", sum(!is.na(final_df$excluded_coords))))
  
  return(final_df[, order_cols])
}


