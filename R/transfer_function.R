#' Annotate BLASTn alignments and determine transfer direction (MT <-> PT)
#'
#' @param fasta_mt Path to mitochondrial FASTA file (Query).
#' @param fasta_pt Path to plastid FASTA file (Subject).
#' @param bed_mt Mitochondrial BED file (4 columns: chrom, start, end, name).
#' @param bed_pt Plastid BED file (4 columns: chrom, start, end, name).
#' @param evalue_cut_off Maximum e-value threshold (default: 1e-06).
#' @param gene_buffer Min difference in gene_perc (default: 20).
#' @param trans_buffer Min difference in transfer_perc (default: 20).
#' @export
transfer_function <- function(fasta_mt, fasta_pt, bed_mt, bed_pt, 
                              evalue_cut_off = 1e-06, 
                              gene_buffer = 20, 
                              trans_buffer = 20) {
  
  if (missing(fasta_mt) || missing(fasta_pt) || missing(bed_mt) || missing(bed_pt)) {
    stop("All input files are required.")
  }
  
  message("Running BLASTn (MT vs PT)...")
  blast_n <- metablastr::blast_nucleotide_to_nucleotide(
    query = fasta_mt, subject = fasta_pt,
    db.import = FALSE, task = "blastn", evalue = evalue_cut_off
  )
  
  if (nrow(blast_n) == 0) {
    message("No BLAST hits found.")
    return(blast_n)
  }
  
  blast_n$q_start_fix <- pmin(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
  blast_n$q_end_fix   <- pmax(as.numeric(blast_n$q_start), as.numeric(blast_n$q_end))
  blast_n$s_start_fix <- pmin(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))
  blast_n$s_end_fix   <- pmax(as.numeric(blast_n$s_start), as.numeric(blast_n$s_end))
  blast_n$alig_length <- as.numeric(blast_n$alig_length)
  
  trna_regex <- "^trn|tRNA"
  
  annotate_and_get_metrics <- function(blast_df, bed_df, prefix) {
    colnames(bed_df)[1:4] <- c("chrom", "start", "end", "name")
    
    chrom_col   <- if(prefix == "q") "query_id" else "subject_id"
    b_start_col <- if(prefix == "q") "q_start_fix" else "s_start_fix"
    b_end_col   <- if(prefix == "q") "q_end_fix" else "s_end_fix"
    
    lapply(1:nrow(blast_df), function(i) {
      hit <- blast_df[i, ]
      overlaps <- bed_df[bed_df$chrom == hit[[chrom_col]] & 
                           bed_df$start < hit[[b_end_col]] & 
                           bed_df$end > hit[[b_start_col]], ]
      
      if (nrow(overlaps) == 0) {
        return(list(text = "none", max_g_perc = 0, max_t_perc = 0, only_trna = FALSE))
      }
      
      is_trna <- grepl(trna_regex, overlaps$name, ignore.case = TRUE)
      only_trna <- all(is_trna)
      
      gene_results <- list()
      for(j in 1:nrow(overlaps)) {
        g_start <- as.numeric(overlaps$start[j])
        g_end   <- as.numeric(overlaps$end[j])
        g_len   <- g_end - g_start
        
        ov_start <- max(g_start, hit[[b_start_col]])
        ov_end   <- min(g_end, hit[[b_end_col]])
        ov_len   <- max(0, ov_end - ov_start)
        
        g_perc <- (ov_len / g_len) * 100
        t_perc <- (ov_len / hit$alig_length) * 100
        
        gene_results[[j]] <- c(g_perc = g_perc, t_perc = t_perc, g_len = g_len, ov_len = ov_len, name = overlaps$name[j])
      }
      
      all_g_perc <- sapply(gene_results, function(x) as.numeric(x["g_perc"]))
      all_t_perc <- sapply(gene_results, function(x) as.numeric(x["t_perc"]))
      
      res_text <- sapply(gene_results, function(x) {
        paste0(x["name"], " (g_len=", x["g_len"], 
               ", ov_len=", round(as.numeric(x["ov_len"]), 0), 
               ", g_perc=", round(as.numeric(x["g_perc"]), 1), "%, ",
               "t_perc=", round(as.numeric(x["t_perc"]), 1), "%)")
      })
      
      return(list(text = paste(res_text, collapse = "; "), 
                  max_g_perc = max(all_g_perc), 
                  max_t_perc = max(all_t_perc), 
                  only_trna = only_trna))
    })
  }
  
  message("Annotating Mitochondrial side...")
  mt_ann <- annotate_and_get_metrics(blast_n, bed_mt, "q")
  message("Annotating Plastid side...")
  pt_ann <- annotate_and_get_metrics(blast_n, bed_pt, "s")
  
  blast_n$mt_genes <- sapply(mt_ann, function(x) x$text)
  blast_n$pt_genes <- sapply(pt_ann, function(x) x$text)
  
  blast_n$direction <- "unknown"
  blast_n$reason <- ""
  
  for (i in 1:nrow(blast_n)) {
    m <- mt_ann[[i]]; p <- pt_ann[[i]]
    
    if (m$only_trna && p$only_trna) {
      blast_n$direction[i] <- "unknown"
      blast_n$reason[i] <- "tRNA only (non-diagnostic)"
      next
    }
    
    g_diff <- abs(m$max_g_perc - p$max_g_perc)
    if (g_diff >= gene_buffer) {
      blast_n$direction[i] <- if(m$max_g_perc > p$max_g_perc) "MT -> PT" else "PT -> MT"
      blast_n$reason[i] <- sprintf("Gene Completeness (diff %.1f%% > %d%%)", g_diff, gene_buffer)
    } else {
      t_diff <- abs(m$max_t_perc - p$max_t_perc)
      if (t_diff >= trans_buffer) {
        blast_n$direction[i] <- if(m$max_t_perc > p$max_t_perc) "MT -> PT" else "PT -> MT"
        blast_n$reason[i] <- sprintf("Transfer Dominance (diff %.1f%% > %d%%)", t_diff, trans_buffer)
      } else {
        blast_n$direction[i] <- "unknown"
        if (m$text == "none" && p$text == "none") {
          blast_n$reason[i] <- "Intergenic region (no features)"
        } else {
          blast_n$reason[i] <- "Ambiguous metrics (below buffers)"
        }
      }
    }
  }
  
  blast_n <- blast_n[, !(names(blast_n) %in% c("q_start_fix", "q_end_fix", "s_start_fix", "s_end_fix"))]
  message("Done!")
  blast_n <- as.data.frame(blast_n)
  return(blast_n)
}