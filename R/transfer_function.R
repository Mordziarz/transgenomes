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
    stop("All input files are required.")
  }

  # 1. Ładowanie plików BED
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
  
  # Regexy dla tRNA i rRNA
  trna_regex <- "^trn|tRNA"
  rrna_regex <- "^rrn|rRNA|[0-9]+S_rRNA"
  nc_regex   <- paste0(trna_regex, "|", rrna_regex)

  # 2. Uruchomienie BLASTn
  message("Running BLASTn (MT vs PT)...")
  blast_res <- metablastr::blast_nucleotide_to_nucleotide(
    query = fasta_mt, subject = fasta_pt,
    db.import = FALSE, task = "blastn", evalue = evalue_cut_off
  )
  
  if (nrow(blast_res) == 0) {
    message("No BLAST hits found.")
    return(NULL)
  }

  # 3. Zmiana nazw kolumn na czytelne (MT i PT)
  colnames(blast_res)[colnames(blast_res) == "query_id"]   <- "mt_id"
  colnames(blast_res)[colnames(blast_res) == "subject_id"] <- "pt_id"
  colnames(blast_res)[colnames(blast_res) == "q_start"]    <- "mt_start"
  colnames(blast_res)[colnames(blast_res) == "q_end"]      <- "mt_end"
  colnames(blast_res)[colnames(blast_res) == "q_len"]      <- "mt_total_len"
  colnames(blast_res)[colnames(blast_res) == "s_start"]    <- "pt_start"
  colnames(blast_res)[colnames(blast_res) == "s_end"]      <- "pt_end"
  colnames(blast_res)[colnames(blast_res) == "s_len"]      <- "pt_total_len"

  # 4. Filtrowanie i czyszczenie
  blast_res$alig_length <- as.numeric(blast_res$alig_length)
  blast_res$perc_identity <- as.numeric(blast_res$perc_identity)
  blast_res <- blast_res[blast_res$alig_length >= min_length & blast_res$perc_identity >= min_identity, ]
  
  if (nrow(blast_res) == 0) return(NULL)

  # Naprawa współrzędnych (start < end)
  mt_s <- pmin(as.numeric(blast_res$mt_start), as.numeric(blast_res$mt_end))
  mt_e <- pmax(as.numeric(blast_res$mt_start), as.numeric(blast_res$mt_end))
  pt_s <- pmin(as.numeric(blast_res$pt_start), as.numeric(blast_res$pt_end))
  pt_e <- pmax(as.numeric(blast_res$pt_start), as.numeric(blast_res$pt_end))

  # 5. Funkcja do scalania nakładających się regionów
  merge_intervals <- function(starts, ends) {
    if (length(starts) == 0) return(0)
    ord <- order(starts)
    s <- starts[ord]; e <- ends[ord]
    m_start <- s[1]; m_end <- e[1]; total <- 0
    if(length(s) > 1) {
      for (i in 2:length(s)) {
        if (s[i] <= m_end) { m_end <- max(m_end, e[i]) } 
        else { total <- total + (m_end - m_start); m_start <- s[i]; m_end <- e[i] }
      }
    }
    return(total + (m_end - m_start))
  }

  # 6. Adnotacja
  annotate_organelle <- function(df_idx, bed, organelle_starts, organelle_ends, ids) {
    hit_start <- organelle_starts[df_idx]
    hit_end   <- organelle_ends[df_idx]
    hit_id    <- tolower(ids[df_idx])
    hit_len   <- blast_res$alig_length[df_idx]
    
    overlaps <- bed[bed$chrom == hit_id & bed$start < hit_end & bed$end > hit_start, ]
    
    if (nrow(overlaps) == 0) {
      return(list(text = "none", sum_g_perc = 0, sum_t_perc = 0, has_protein = FALSE, only_nc = FALSE))
    }
    
    is_nc <- grepl(nc_regex, overlaps$name, ignore.case = TRUE)
    has_protein <- any(!is_nc)
    
    # Statystyki: tylko białka jeśli są, inaczej nc
    stats_overlaps <- if(has_protein) overlaps[!is_nc, ] else overlaps
    rel_ov_starts  <- pmax(stats_overlaps$start, hit_start)
    rel_ov_ends    <- pmin(stats_overlaps$end, hit_end)
    
    u_ov_len <- merge_intervals(rel_ov_starts, rel_ov_ends)
    u_g_len  <- merge_intervals(stats_overlaps$start, stats_overlaps$end)
    
    text <- paste(sapply(1:nrow(overlaps), function(j) {
      g_l <- overlaps$end[j] - overlaps$start[j]
      ov  <- max(0, min(overlaps$end[j], hit_end) - max(overlaps$start[j], hit_start))
      paste0(overlaps$name[j], " (g_len=", g_l, ", ov_len=", round(ov, 0), 
             ", g_perc=", round((ov/g_l)*100, 1), "%, t_perc=", round((ov/hit_len)*100, 1), "%)")
    }), collapse = "; ")
    
    return(list(text = text, sum_g_perc = (u_ov_len / u_g_len) * 100, sum_t_perc = (u_ov_len / hit_len) * 100, 
                has_protein = has_protein, only_nc = all(is_nc)))
  }

  mt_ann <- lapply(1:nrow(blast_res), function(i) annotate_organelle(i, b_mt, mt_s, mt_e, blast_res$mt_id))
  pt_ann <- lapply(1:nrow(blast_res), function(i) annotate_organelle(i, b_pt, pt_s, pt_e, blast_res$pt_id))
  
  blast_res$mt_genes <- sapply(mt_ann, function(x) x$text)
  blast_res$pt_genes <- sapply(pt_ann, function(x) x$text)
  blast_res$direction <- "Unidentified"

  # 7. Logika kierunku
  for (i in 1:nrow(blast_res)) {
    m <- mt_ann[[i]]; p <- pt_ann[[i]]
    
    # Jeśli obie strony to tylko niekodujące -> Unidentified
    if (m$only_nc && p$only_nc) {
      blast_res$direction[i] <- "Unidentified"
      next
    }
    
    # Bezwzględne pierwszeństwo białka
    if (m$has_protein && !p$has_protein) {
      blast_res$direction[i] <- "MT -> PT"
      next
    }
    if (!m$has_protein && p$has_protein) {
      blast_res$direction[i] <- "PT -> MT"
      next
    }
    
    # Jeśli obie strony mają białka -> Statystyki
    if (m$has_protein && p$has_protein) {
      g_diff <- abs(m$sum_g_perc - p$sum_g_perc)
      if (g_diff >= gene_buffer) {
        blast_res$direction[i] <- if(m$sum_g_perc > p$sum_g_perc) "MT -> PT" else "PT -> MT"
      } else {
        t_diff <- abs(m$sum_t_perc - p$sum_t_perc)
        if (t_diff >= trans_buffer) {
          blast_res$direction[i] <- if(m$sum_t_perc > p$sum_t_perc) "MT -> PT" else "PT -> MT"
        }
      }
    }
  }

  message("Analysis complete.")
  return(as.data.frame(blast_res))
}