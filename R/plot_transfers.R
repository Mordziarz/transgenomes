#' Circular Visualization of Genomic Transfers (MT-PT)
#'
#' Generates a high-quality Circos plot visualizing genomic alignments and transfer directions
#' between mitochondrial (MT) and plastid (PT) genomes. The plot includes colored genome 
#' sectors, direction-specific links, and dual legends.
#'
#' @param transfer_function_out Data frame containing alignment data. Must include 
#' columns: `query_id`, `subject_id`, `q_start`, `q_end`, `s_start`, `s_end`, `q_len`, 
#' `s_len`, and `direction`.
#' 
#' @param gap Numeric. Distance between genome segments in degrees. Default is 10.
#' @param start_degree Numeric. Starting rotation of the plot (0-360). Default is 90.
#' @param label_cex Numeric. Character expansion (size) for chromosome/sector labels. Default is 0.8.
#' @param label_dist Numeric. Distance of the labels from the genome track. Default is 1.5.
#' @param transparency Numeric. Opacity of the ribbons (0 = transparent, 1 = opaque). Default is 0.5.
#' @param mt_sector_col Character. Color for the mitochondrial genome segments. Default is "orange".
#' @param pt_sector_col Character. Color for the plastid genome segments. Default is "darkgreen".
#' @param mt_to_pt_col Character. Color for ribbons representing MT to PT transfers. Default is "firebrick1".
#' @param pt_to_mt_col Character. Color for ribbons representing PT to MT transfers. Default is "dodgerblue1".
#' @param unknown_col Character. Color for ribbons where direction is unknown. Default is "grey80".
#'
#' @details 
#' The function uses the `circlize` package to initialize the circular layout. It automatically
#' identifies unique sequences for both genomes and assigns colors based on the `direction` 
#' column. Labels are curved to follow the track's arc using `bending.inside` facing.
#'
#' @return Invisible NULL. The function generates a plot in the active graphics device.
#' 
#' @importFrom circlize circos.clear circos.par circos.initialize circos.track circos.rect circos.text circos.axis circos.genomicLink CELL_META
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics legend
#' @export
#' 

plot_transfers <- function(transfer_function_out, 
                           gap = 10, 
                           start_degree = 90, 
                           label_cex = 0.8,
                           label_dist = 1.5, 
                           transparency = 0.5,
                           mt_sector_col = "orange",
                           pt_sector_col = "darkgreen",
                           mt_to_pt_col = "firebrick1",
                           pt_to_mt_col = "dodgerblue1",
                           unknown_col = "grey80") {
  
  if (missing(transfer_function_out)) {
    stop("Argument 'transfer_function_out' is missing.")
  }
  
  get_alpha_col <- function(col, alpha) {
    rgb_val <- col2rgb(col)
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255, alpha = alpha * 255)
  }
  
  col_mt_pt <- get_alpha_col(mt_to_pt_col, transparency)
  col_pt_mt <- get_alpha_col(pt_to_mt_col, transparency)
  col_unk   <- get_alpha_col(unknown_col, transparency)
  
  link_colors <- ifelse(transfer_function_out$direction == "MT -> PT", col_mt_pt,
                        ifelse(transfer_function_out$direction == "PT -> MT", col_pt_mt, col_unk))
  
  mt_ids <- unique(transfer_function_out$query_id)
  pt_ids <- unique(transfer_function_out$subject_id)
  
  mt_chroms <- unique(transfer_function_out[, c("query_id", "q_len")])
  colnames(mt_chroms) <- c("id", "len")
  pt_chroms <- unique(transfer_function_out[, c("subject_id", "s_len")])
  colnames(pt_chroms) <- c("id", "len")
  
  full_genome_info <- rbind(mt_chroms, pt_chroms)
  full_genome_info <- full_genome_info[!duplicated(full_genome_info$id), ]
  
  sector_colors <- ifelse(full_genome_info$id %in% mt_ids, mt_sector_col, pt_sector_col)
  names(sector_colors) <- full_genome_info$id
  
  circlize::circos.clear()
  
  circlize::circos.par(
    start.degree = start_degree, 
    gap.after = rep(gap, nrow(full_genome_info)),
    canvas.xlim = c(-1.2, 1.2), 
    canvas.ylim = c(-1.2, 1.2)
  )
  
  circlize::circos.initialize(
    factors = full_genome_info$id, 
    xlim = cbind(rep(0, nrow(full_genome_info)), full_genome_info$len)
  )
  
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = circlize::CELL_META$sector.index
    s_col = sector_colors[chr]
    
    circlize::circos.rect(circlize::CELL_META$xlim[1], 0, 
                          circlize::CELL_META$xlim[2], 1, 
                          col = s_col, border = "black")
    
    circlize::circos.text(
      circlize::CELL_META$xcenter, 
      circlize::CELL_META$ylim[2] + label_dist, 
      chr, 
      cex = label_cex, 
      facing = "bending.inside", 
      niceFacing = TRUE
    )
    
    circlize::circos.axis(labels.cex = 0.4, minor.ticks = 1)
    
  }, bg.border = NA, track.height = 0.1) # Grubszy pasek dla lepszej widoczności kolorów
  
  links_q_fix <- data.frame(
    chr = transfer_function_out$query_id,
    start = pmin(as.numeric(transfer_function_out$q_start), as.numeric(transfer_function_out$q_end)),
    end = pmax(as.numeric(transfer_function_out$q_start), as.numeric(transfer_function_out$q_end))
  )
  links_s_fix <- data.frame(
    chr = transfer_function_out$subject_id,
    start = pmin(as.numeric(transfer_function_out$s_start), as.numeric(transfer_function_out$s_end)),
    end = pmax(as.numeric(transfer_function_out$s_start), as.numeric(transfer_function_out$s_end))
  )
  
  circlize::circos.genomicLink(links_q_fix, links_s_fix, col = link_colors, border = NA)
  
  graphics::legend(x = -1.3, y = -0.7, 
                   legend = c("MT -> PT", "PT -> MT", "Unknown"), 
                   fill = c(mt_to_pt_col, pt_to_mt_col, unknown_col), 
                   title = "Transfer Direction", bty = "n", cex = 0.7)
  
  graphics::legend(x = -1.3, y = -1.0, 
                   legend = c("Mitochondrion", "Plastid"), 
                   fill = c(mt_sector_col, pt_sector_col), 
                   title = "Genomes", bty = "n", cex = 0.7)
  
  message("Plot created with maximized circle and side-by-side legends.")
}