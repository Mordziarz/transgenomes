plot_transfers2 <- function(data, 
                                   normalization = TRUE,
                                   chromosome_gap = 0.05, 
                                   transparency = 0.5,
                                   mt_sector_col = "orange", 
                                   pt_sector_col = "darkgreen",
                                   mt_to_pt_col = "firebrick1", 
                                   pt_to_mt_col = "dodgerblue1", 
                                   unidentified_col = "grey80") {
  
  data <- data %>%
    mutate(across(c(mt_start, mt_end, pt_start, pt_end, mt_total_len, pt_total_len), as.numeric))
  
  if (!"direction" %in% colnames(data)) data$direction <- "PT -> MT"
  
  mt_raw <- data %>% group_by(mt_id) %>% summarise(len = first(mt_total_len), .groups = 'drop')
  pt_raw <- data %>% group_by(pt_id) %>% summarise(len = first(pt_total_len), .groups = 'drop')
  
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
  
  data_proc <- data %>%
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