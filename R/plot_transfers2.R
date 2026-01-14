#' Linear Visualization of Genomic Transfers
#'
#' Generates a linear synteny plot visualizing genomic transfers 
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
                            transparency = 0.6,
                            mt_sector_col = "orange", 
                            pt_sector_col = "darkgreen",
                            mt_to_pt_col = "firebrick1", 
                            pt_to_mt_col = "dodgerblue1", 
                            unidentified_col = "grey80") {
  
  library(dplyr)
  library(ggplot2)
  
  data_clean <- transfer_function_out %>%
    mutate(across(c(mt_start, mt_end, pt_start, pt_end, mt_total_len, pt_total_len), as.numeric))
  
  if (!"direction" %in% colnames(data_clean)) data_clean$direction <- "PT -> MT"
  
  mt_raw <- data_clean %>% 
    group_by(mt_id) %>% 
    summarise(len = max(c(mt_total_len, mt_start, mt_end), na.rm = TRUE), .groups = 'drop')
  
  pt_raw <- data_clean %>% 
    group_by(pt_id) %>% 
    summarise(len = max(c(pt_total_len, pt_start, pt_end), na.rm = TRUE), .groups = 'drop')
  
  mt_spacer <- if(normalization) chromosome_gap * sum(mt_raw$len) else chromosome_gap * (sum(mt_raw$len)/10)
  pt_spacer <- if(normalization) chromosome_gap * sum(pt_raw$len) else chromosome_gap * (sum(pt_raw$len)/10)
  
  mt_info <- mt_raw %>% mutate(offset = c(0, head(cumsum(len + mt_spacer), -1)))
  pt_info <- pt_raw %>% mutate(offset = c(0, head(cumsum(len + pt_spacer), -1)))
  
  total_mt_axis <- sum(mt_info$len) + (nrow(mt_info)-1) * mt_spacer
  total_pt_axis <- sum(pt_info$len) + (nrow(pt_info)-1) * pt_spacer
  
  scale_pos <- function(pos, offset, total_axis, norm) {
    global <- pos + offset
    if(norm) return((global / total_axis) * 100)
    return(global)
  }
  
  data_proc <- data_clean %>%
    left_join(mt_info %>% select(mt_id, mt_off = offset, mt_actual_len = len), by = "mt_id") %>%
    left_join(pt_info %>% select(pt_id, pt_off = offset, pt_actual_len = len), by = "pt_id") %>%
    mutate(
      m_start_scaled = scale_pos(mt_start, mt_off, total_mt_axis, normalization),
      m_end_scaled   = scale_pos(mt_end, mt_off, total_mt_axis, normalization),
      p_start_scaled = scale_pos(pt_start, pt_off, total_pt_axis, normalization),
      p_end_scaled   = scale_pos(pt_end, pt_off, total_pt_axis, normalization),
      m_mid = (m_start_scaled + m_end_scaled) / 2,
      p_mid = (p_start_scaled + p_end_scaled) / 2,
      transfer_size = abs(mt_end - mt_start)
    )
  
  ggplot(data_proc) +
    geom_segment(aes(x = m_mid, y = 5, 
                     xend = p_mid, yend = 1, 
                     linewidth = transfer_size, 
                     color = direction), 
                 alpha = transparency) + 
    
    geom_segment(data = mt_info, 
                 aes(x = scale_pos(0, offset, total_mt_axis, normalization), 
                     xend = scale_pos(len, offset, total_mt_axis, normalization), 
                     y = 5, yend = 5), 
                 linewidth = 5, color = mt_sector_col, lineend = "butt") +
    
    geom_segment(data = pt_info, 
                 aes(x = scale_pos(0, offset, total_pt_axis, normalization), 
                     xend = scale_pos(len, offset, total_pt_axis, normalization), 
                     y = 1, yend = 1), 
                 linewidth = 5, color = pt_sector_col, lineend = "butt") +
    
    geom_text(data = mt_info, aes(x = scale_pos(len/2, offset, total_mt_axis, normalization), 
                                  y = 5.4, label = mt_id), fontface = "bold") +
    geom_text(data = pt_info, aes(x = scale_pos(len/2, offset, total_pt_axis, normalization), 
                                  y = 0.6, label = pt_id), fontface = "bold") +
    
    scale_color_manual(values = c("MT -> PT" = mt_to_pt_col, 
                                  "PT -> MT" = pt_to_mt_col, 
                                  "Unidentified" = unidentified_col)) +
    scale_y_continuous(breaks = c(1, 5), labels = c("Plastome", "Mitogenome"), limits = c(0, 6)) + 
    scale_x_continuous(expand = c(0.05, 0.05)) +
    labs(x = ifelse(normalization, "Relative Position (%)", "Position (bp)"), 
         y = "", color = "Transfer Direction") +
    guides(linewidth = "none", color = guide_legend(override.aes = list(linewidth = 3))) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
}