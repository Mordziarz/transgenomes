#' Analyze and Compare GC Content Between Two Sets of DNA Sequences
#'
#' @param extract_regions_1 A DNAStringSet object (e.g., from extract_regions).
#' @param extract_regions_2 A DNAStringSet object (e.g., from extract_regions).
#' @param group1_name Character string, name for the first group (default: "Group 1").
#' @param group2_name Character string, name for the second group (default: "Group 2").
#' @param alpha Significance level for statistical tests (default: 0.05).
#'
#' @return A list containing the plot, test results, and raw GC data.
#' @import Biostrings ggplot2 dplyr rstatix car ggpubr
#' @export

analyze_GC_content <- function(
  extract_regions_1, 
  extract_regions_2, 
  group1_name = "Group 1", 
  group2_name = "Group 2",
  alpha = 0.05
) {
  
  prepare_seqs <- function(x) {
    if (inherits(x, "DNAStringSetList")) return(unlist(x))
    return(x)
  }

  seqs1 <- prepare_seqs(extract_regions_1)
  seqs2 <- prepare_seqs(extract_regions_2)
  
  if(length(seqs1) == 0 || length(seqs2) == 0) {
    stop("One of the input sequence sets is empty.")
  }

  calculate_gc <- function(x) {
    freqs <- Biostrings::letterFrequency(x, letters = "GC", as.prob = TRUE)
    return(as.numeric(freqs) * 100)
  }
  
  gc_data <- data.frame(
    GC = c(calculate_gc(seqs1), calculate_gc(seqs2)),
    group = factor(
      rep(c(group1_name, group2_name), 
          c(length(seqs1), length(seqs2))),
      levels = c(group1_name, group2_name)
    )
  )
  
  select_test <- function(data) {
    norm_test <- data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        n = dplyr::n(),
        p_val = if(dplyr::n() >= 3) shapiro.test(GC)$p.value else 0, # Force non-normal if n < 3
        is_normal = p_val > alpha,
        .groups = "drop"
      )
    
    if(all(norm_test$is_normal)) {
      var_test <- car::leveneTest(GC ~ group, data = data)
      p_var <- var_test$`Pr(>F)`[1]
      
      if(p_var > alpha) {
        test <- rstatix::t_test(GC ~ group, data = data, var.equal = TRUE)
        method <- "t-test"
      } else {
        test <- rstatix::t_test(GC ~ group, data = data, var.equal = FALSE)
        method <- "Welch t-test"
      }
    } else {
      test <- rstatix::wilcox_test(GC ~ group, data = data)
      method <- "Wilcoxon test"
    }
    list(test = test, method = method, norm_test = norm_test)
  }
  
  res_list <- select_test(gc_data)
  
  stat_test <- res_list$test %>%
    rstatix::add_y_position() %>%
    dplyr::mutate(label = paste0(res_list$method, "\np = ", signif(p, 3)))
  
  gc_plot <- ggplot(gc_data, aes(x = group, y = GC, fill = group)) +
    geom_boxplot(alpha = 0.6, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = group), width = 0.15, alpha = 0.5, size = 1.5) +
    ggpubr::stat_pvalue_manual(
      stat_test, 
      label = "label",
      tip.length = 0.02,
      label.size = 3.5
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "GC Content Comparison",
      y = "GC Content (%)",
      x = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black")
    )
  
  return(list(
    plot = gc_plot,
    test_result = res_list$test,
    test_method = res_list$method,
    normality_check = res_list$norm_test,
    gc_data = gc_data
  ))
}