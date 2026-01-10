#' Circular Visualization of Genomic Transfers (MT-PT)
#'
#' Generates a high-quality Circos plot visualizing genomic alignments and transfer directions
#' between mitochondrial (MT) and plastid (PT) genomes.
#'
#' @param transfer_function_out Data frame containing alignment data. Must include 
#' columns: `mt_id`, `pt_id`, `mt_start`, `mt_end`, `pt_start`, `pt_end`, `mt_total_len`, 
#' `pt_total_len`, and `direction`.
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
#' @param undefined_col Character. Color for ribbons where direction is Unidentified. Default is "grey80".

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
                           undefined_col = "grey80") {
  
  if (missing(transfer_function_out)) {
    stop("Argument 'transfer_function_out' is missing.")
  }
  
  get_alpha_col <- function(col, alpha) {
    rgb_val <- col2rgb(col)
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255, alpha = alpha * 255)
  }
  
  col_mt_pt <- get_alpha_col(mt_to_pt_col, transparency)
  col_pt_mt <- get_alpha_col(pt_to_mt_col, transparency)
  col_undef <- get_alpha_col(undefined_col, transparency)
  
  link_colors <- ifelse(transfer_function_out$direction == "MT -> PT", col_mt_pt,
                        ifelse(transfer_function_out$direction == "PT -> MT", col_pt_mt, col_undef))
  
  mt_ids <- unique(transfer_function_out$mt_id)
  pt_ids <- unique(transfer_function_out$pt_id)
  
  mt_chroms <- unique(transfer_function_out[, c("mt_id", "mt_total_len")])
  colnames(mt_chroms) <- c("id", "len")
  
  pt_chroms <- unique(transfer_function_out[, c("pt_id", "pt_total_len")])
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
    canvas.ylim = c(-1.5, 1.0)
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
    
  }, bg.border = NA, track.height = 0.1)
  
  links_mt <- data.frame(
    chr = transfer_function_out$mt_id,
    start = pmin(as.numeric(transfer_function_out$mt_start), as.numeric(transfer_function_out$mt_end)),
    end = pmax(as.numeric(transfer_function_out$mt_start), as.numeric(transfer_function_out$mt_end))
  )
  links_pt <- data.frame(
    chr = transfer_function_out$pt_id,
    start = pmin(as.numeric(transfer_function_out$pt_start), as.numeric(transfer_function_out$pt_end)),
    end = pmax(as.numeric(transfer_function_out$pt_start), as.numeric(transfer_function_out$pt_end))
  )
  
  circlize::circos.genomicLink(links_mt, links_pt, col = link_colors, border = NA)
  
  graphics::legend(x = -1.3, y = -0.7, 
                   legend = c("MT -> PT", "PT -> MT", "Unidentified"), 
                   fill = c(mt_to_pt_col, pt_to_mt_col, undefined_col), 
                   title = expression(bold("Transfer Direction")), bty = "n", cex = 0.7)
  
  graphics::legend(x = -1.3, y = -1.1, 
                   legend = c("Mitogenome", "Plastome"), 
                   fill = c(mt_sector_col, pt_sector_col), 
                   title = expression(bold("Genomes")), bty = "n", cex = 0.7)
  
  message("DONE !!!!! :)")
}