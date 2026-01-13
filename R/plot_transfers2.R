#' Linear Visualization of Genomic Transfers
#'
#' Generates a linear synteny plot (ribbon plot) visualizing genomic transfers 
#' between the mitogenome and the plastome. This function supports multiple 
#' chromosomes on both axes by introducing customizable gaps between them 
#' for enhanced clarity.
#'
#' @param transfer_function_out A data frame containing at least the following columns: 
#'   \code{mt_id}, \code{pt_id}, \code{mt_start}, \code{mt_end}, \code{pt_start}, 
#'   \code{pt_end}, \code{mt_total_len}, \code{pt_total_len}, and optionally \code{direction}.
#' @param normalization Logical. If \code{TRUE}, positions are scaled to percentages (0-100%). 
#'   If \code{FALSE}, raw base pair (bp) units are used. Defaults to \code{TRUE}.
#' @param chromosome_gap Numeric. The size of the gap between adjacent chromosomes. 
#'   When \code{normalization} is \code{TRUE}, this is a fraction of the total axis length 
#'   (e.g., 0.05 = 5%). Defaults to 0.05.
#' @param transparency Numeric. The opacity of the ribbon fill (0 to 1). Defaults to 0.5.
#' @param mt_sector_col Character. Color for the mitochondrial genome segments. Defaults to "orange".
#' @param pt_sector_col Character. Color for the plastid genome segments. Defaults to "darkgreen".
#' @param mt_to_pt_col Character. Color for ribbons representing MT to PT transfers. Defaults to "firebrick1".
#' @param pt_to_mt_col Character. Color for ribbons representing PT to MT transfers. Defaults to "dodgerblue1".
#' @param unidentified_col Character. Color for ribbons with an unknown transfer direction. Defaults to "grey80".
#'
#'
#' @examples
#' # Assuming 'my_data' is your data frame:
#' # plot_transfers2(transfer_function_out = my_data, normalization = TRUE)

plot_transfers2 <- function(transfer_function_out, 
                            normalization = TRUE,
                            chromosome_gap = 0.05, 
                            transparency = 0.5,
                            mt_sector_col = "orange", 
                            pt_sector_col = "darkgreen",
                            mt_to_pt_col = "firebrick1", 
                            pt_to_mt_col = "dodgerblue1", 
                            unidentified_col = "grey80") {
  
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required for this function.")
  }

  data_clean <- transfer_function_out %>%
    mutate(across(c(mt_start, mt_end, pt_start, pt_end, mt_total_len, pt_total_len), as.numeric))
  
  if (!"direction" %in% colnames(data_clean)) data_clean$direction <- "PT -> MT"
  
  mt_raw <- data_clean %>% group_by(mt_id) %>% summarise(len = first(mt_total_len), .groups = 'drop')
  pt_raw <- data_clean %>% group_by(pt_id) %>% summarise(len = first(pt_total_len), .groups = 'drop')
  
  mt_spacer <- if(normalization) chromosome_gap * sum(mt_raw$len) else chromosome_gap * (sum(mt_raw$len)/10)
  pt_spacer <- if(normalization) chromosome_gap * sum(pt_raw$len) else chromosome_gap * (sum(pt_raw$len)/10)
  
  mt_info <- mt_raw %>%
    mutate(offset = c(0, head(cumsum(len + mt_spacer), -1)))
  pt_info <- pt_raw %>%
    mutate(offset = c(0, head(cumsum(len + pt_spacer), -1)))
  
  total_mt_axis <- sum(mt_info$len) + (nrow(mt_info)-1) * mt_spacer
  total_pt_axis <- sum(pt_info$len) + (nrow(pt_info)-1) * pt_spacer
  
  scale_pos <- function(pos, offset, total_axis, norm) {
    global <- pos + offset
    if(norm) return((global / total_axis) * 100)
    return(global)
  }
  
  data_proc <- data_clean %>%
    left_join(mt_info %>% select(mt_id, mt_off = offset), by = "mt_id") %>%
    left_join(pt_info %>% select(pt_id, pt_off = offset), by = "pt_id") %>%
    mutate(m_s = scale_pos(pmin(mt_start, mt_end), mt_off, total_mt_axis, normalization),
           m_e = scale_pos(pmax(mt_start, mt_end), mt_off, total_mt_axis, normalization),
           p_s = scale_pos(pmin(pt_start, pt_end), pt_off, total_pt_axis, normalization),
           p_e = scale_pos(pmax(pt_start, pt_end), pt_off, total_pt_axis, normalization))
  
  ribbon_data <- data.frame(
    x = c(data_proc$m_s, data_proc$m_e, data_proc$p_e, data_proc$p_s),
    y = c(rep(5, nrow(data_proc)), rep(5, nrow(data_proc)), 
          rep(1, nrow(data_proc)), rep(1, nrow(data_proc))),
    group = rep(1:nrow(data_proc), 4),
    direction = rep(data_proc$direction, 4)
  )

  ggplot() +
    geom_diagonal_wide(data = ribbon_data, 
                       aes(x = x, y = y, group = group, fill = direction, color = direction), 
                       alpha = transparency, 
                       linewidth = 0.2) + 
    
    geom_segment(data = mt_info, 
                 aes(x = scale_pos(0, offset, total_mt_axis, normalization), 
                     xend = scale_pos(len, offset, total_mt_axis, normalization), y = 5, yend = 5), 
                 linewidth = 3, color = mt_sector_col, lineend = "butt") +
    
    geom_segment(data = pt_info, 
                 aes(x = scale_pos(0, offset, total_pt_axis, normalization), 
                     xend = scale_pos(len, offset, total_pt_axis, normalization), y = 1, yend = 1), 
                 linewidth = 3, color = pt_sector_col, lineend = "butt") +
    
    geom_text(data = mt_info, aes(x = scale_pos(len/2, offset, total_mt_axis, normalization), 
                                  y = 5.4, label = mt_id), size = 4.5, fontface = "italic") +
    geom_text(data = pt_info, aes(x = scale_pos(len/2, offset, total_pt_axis, normalization), 
                                  y = 0.6, label = pt_id), size = 4.5, fontface = "italic") +
    
    scale_fill_manual(values = c("MT -> PT" = mt_to_pt_col, 
                                 "PT -> MT" = pt_to_mt_col, 
                                 "Unidentified" = unidentified_col)) +
    scale_color_manual(values = c("MT -> PT" = mt_to_pt_col, 
                                  "PT -> MT" = pt_to_mt_col, 
                                  "Unidentified" = unidentified_col)) +
    
    scale_y_continuous(breaks = c(1, 5), labels = c("Plastome", "Mitogenome"), limits = c(0, 6)) + 
    scale_x_continuous(expand = c(0.02, 0.02)) +
    
    labs(x = ifelse(normalization, "Relative Position (%)", "Position (bp)"), 
         y = "", fill = "Transfer Direction", color = "Transfer Direction") +
    
    theme_minimal() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.major.y = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )
}