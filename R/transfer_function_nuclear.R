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
    return(data.frame(chrom = tolower(as.character(df[[1]])), 
                      start = as.numeric(df[[2]]), 
                      end = as.numeric(df[[3]]), 
                      name = as.character(df[[4]]), 
                      stringsAsFactors = FALSE))
  }
  
  b_mt <- load_bed(bed_mt); b_pt <- load_bed(bed_pt); b_nuc <- load_bed(bed_nuc)
  

  rna_regex <- "^trn|^rrn|tRNA|rRNA|[0-9]+S_rRNA"

  run_blast <- function(q, s, lbl_q, lbl_s) {
    message(sprintf("Running BLASTn: %s vs %s...", lbl_q, lbl_s))
    res <- metablastr::blast_nucleotide_to_nucleotide(query = q, subject = s, 
                                                      evalue = evalue_cut_off)
    if (nrow(res) == 0) return(NULL)
    
    res$alig_length <- as.numeric(res$alig_length)
    res$perc_identity <- as.numeric(res$perc_identity)
    res$bit_score <- as.numeric(res$bit_score)
    
    res <- res[res$alig_length >= min_length & res$perc_identity >= min_identity, ]
    return(res)
  }

  raw_mt_pt  <- run_blast(fasta_mt, fasta_pt, "MT", "PT")
  raw_mt_nuc <- run_blast(fasta_mt, fasta_nuc, "MT", "NUC")
  raw_pt_nuc <- run_blast(fasta_pt, fasta_nuc, "PT", "NUC")

  process_results <- function(blast_df, q_bed, s_bed, q_lbl, s_lbl, validation_df = NULL) {
    if (is.null(blast_df) || nrow(blast_df) == 0) return(NULL)
    
    blast_df$q_start_fix <- pmin(as.numeric(blast_df$q_start), as.numeric(blast_df$q_end))
    blast_df$q_end_fix   <- pmax(as.numeric(blast_df$q_start), as.numeric(blast_df$q_end))
    blast_df$s_start_fix <- pmin(as.numeric(blast_df$s_start), as.numeric(blast_df$s_end))
    blast_df$s_end_fix   <- pmax(as.numeric(blast_df$s_start), as.numeric(blast_df$s_end))
    
    results <- lapply(1:nrow(blast_df), function(i) {
      hit <- blast_df[i, ]
      
      if (!is.null(validation_df)) {
        cross_match <- validation_df[validation_df$query_id == hit$query_id & 
                                     abs(as.numeric(validation_df$q_start) - as.numeric(hit$q_start)) < 50, ]
        if (nrow(cross_match) > 0 && max(cross_match$bit_score) > hit$bit_score) {
          hit$direction <- paste0("Excluded: Stronger match in 3rd genome")
          hit$q_genes <- "N/A"; hit$s_genes <- "N/A"
          return(hit)
        }
      }

      get_overlap_info <- function(bed, chrom, start, end, alig_len) {
        ov <- bed[bed$chrom == tolower(chrom) & bed$start < end & bed$end > start, ]
        if (nrow(ov) == 0) return(list(summary = "none", g_sum = 0, t_sum = 0, only_rna = FALSE))
        
        is_rna <- grepl(rna_regex, ov$name, ignore.case = TRUE)
        
        target_ov <- if (any(!is_rna)) ov[!is_rna, ] else ov
        
        gene_stats <- sapply(1:nrow(target_ov), function(j) {
          g_len <- target_ov$end[j] - target_ov$start[j]
          ov_len <- max(0, min(target_ov$end[j], end) - max(target_ov$start[j], start))
          c(g_perc = (ov_len/g_len)*100, t_perc = (ov_len/alig_len)*100)
        })
        
        return(list(summary = paste(target_ov$name, collapse = "; "), 
                    g_sum = sum(gene_stats["g_perc", ]), 
                    t_sum = sum(gene_stats["t_perc", ]), 
                    only_rna = all(is_rna)))
      }

      q_info <- get_overlap_info(q_bed, hit$query_id, hit$q_start_fix, hit$q_end_fix, hit$alig_length)
      s_info <- get_overlap_info(s_bed, hit$subject_id, hit$s_start_fix, hit$s_end_fix, hit$alig_length)
      
      hit$q_genes <- q_info$summary
      hit$s_genes <- s_info$summary
      
      if (q_info$only_rna && s_info$only_rna) {
        hit$direction <- "unknown (RNA-only)"
      } else {
        g_diff <- abs(q_info$g_sum - s_info$g_sum)
        if (g_diff >= gene_buffer) {
          hit$direction <- if(q_info$g_sum > s_info$g_sum) paste(q_lbl, "->", s_lbl) else paste(s_lbl, "->", q_lbl)
        } else {
          t_diff <- abs(q_info$t_sum - s_info$t_sum)
          if (t_diff >= trans_buffer) {
            hit$direction <- if(q_info$t_sum > s_info$t_sum) paste(q_lbl, "->", s_lbl) else paste(s_lbl, "->", q_lbl)
          } else {
            hit$direction <- "unknown"
          }
        }
      }
      return(hit)
    })
    
    res_df <- as.data.frame(do.call(rbind, results))
    cols_to_keep <- c("query_id", "subject_id", "perc_identity", "alig_length", "bit_score", "q_genes", "s_genes", "direction")
    return(res_df[, cols_to_keep])
  }
  
  final_mt_pt  <- process_results(raw_mt_pt, b_mt, b_pt, "MT", "PT", raw_mt_nuc)
  
  final_mt_nuc <- process_results(raw_mt_nuc, b_mt, b_nuc, "MT", "NUC", raw_mt_pt)
  
  final_pt_nuc <- process_results(raw_pt_nuc, b_pt, b_nuc, "PT", "NUC", raw_mt_pt)

  all_results <- do.call(rbind, list(final_mt_pt, final_mt_nuc, final_pt_nuc))
  
  message("Analysis complete. Check 'direction' column for validated events.")
  return(as.data.frame(all_results))
}