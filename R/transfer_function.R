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
    stop("All input files (FASTA and BED) are required.")
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

  b_mt  <- load_bed_local(bed_mt)
  b_pt  <- load_bed_local(bed_pt)
  
  trna_regex <- "^trn|tRNA"
  rrna_regex <- "^rrn|rRNA|[0-9]+S_rRNA"
  non_coding_regex <- paste0(trna_regex, "|", rrna_regex)

  message("Running BLASTn...")
  blast_n <- metablastr::blast_nucleotide_to_nucleotide(
    query = fasta_mt, subject = fasta_pt,
    db.import = FALSE, task = "blastn", evalue = evalue_cut_off
  )
  
  if (nrow(blast_n) == 0) {
    message("No BLAST hits found.")
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
        return(list(text = "none", sum_g_perc = 0, sum_t_perc = 0, has_protein_coding = FALSE))
      }
      
      is_non_coding <- grepl(non_coding_regex, overlaps$name, ignore.case = TRUE)
      has_protein_coding <- any(!is_non_coding)
      
      hit_range <- hit[[b_start_col]]:hit[[b_end_col]]
      all_gene_results <- list()
      protein_ov_bases <- c()
      
      for(j in 1:nrow(overlaps)) {
        g_start <- overlaps$start[j]
        g_end <- overlaps$end[j]
        
        g_len  <- (g_end - g_start) + 1
        g_range <- g_start:g_end
        ov_bases <- intersect(hit_range, g_range)
        ov_len <- length(ov_bases)
        
        g_perc <- min(100, (ov_len / g_len) * 100)
        t_perc_single <- (ov_len / hit$alig_length) * 100
        
        all_gene_results[[j]] <- list(
          name = overlaps$name[j], 
          g_len = g_len, 
          ov_len = ov_len, 
          g_perc = g_perc, 
          t_perc = t_perc_single,
          is_nc = is_non_coding[j]
        )
        
        if(!is_non_coding[j]) {
          protein_ov_bases <- union(protein_ov_bases, ov_bases)
        }
      }

      protein_genes <- Filter(function(x) !x$is_nc, all_gene_results)
      sum_g_perc <- sum(sapply(protein_genes, function(x) x$g_perc))
      sum_t_perc <- (length(protein_ov_bases) / hit$alig_length) * 100
      
      res_text <- sapply(all_gene_results, function(x) {
        paste0(x$name, " (g_len=", x$g_len, 
               ", ov_len=", round(x$ov_len, 0), 
               ", g_perc=", round(x$g_perc, 1), "%, ",
               "t_perc=", round(x$t_perc, 1), "%)")
      })
      
      return(list(text = paste(res_text, collapse = "; "), 
                  sum_g_perc = sum_g_perc, 
                  sum_t_perc = sum_t_perc, 
                  has_protein_coding = has_protein_coding))
    })
  }

  message("Annotating results and determining direction...")
  mt_ann <- annotate_and_get_metrics(blast_n, b_mt, "q")
  pt_ann <- annotate_and_get_metrics(blast_n, b_pt, "s")
  
  blast_n$mt_genes <- sapply(mt_ann, function(x) x$text)
  blast_n$pt_genes <- sapply(pt_ann, function(x) x$text)
  blast_n$direction <- "Unidentified"

  for (i in 1:nrow(blast_n)) {
    m <- mt_ann[[i]]; p <- pt_ann[[i]]
    
    if (!m$has_protein_coding && !p$has_protein_coding) {
      next
    }
    
    if (m$has_protein_coding && !p$has_protein_coding) {
      blast_n$direction[i] <- "MT -> PT"
      next
    }
    if (!m$has_protein_coding && p$has_protein_coding) {
      blast_n$direction[i] <- "PT -> MT"
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
  
  message("Done!")
  return(as.data.frame(blast_n))
}