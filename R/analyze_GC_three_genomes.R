#' Analyze and Compare GC Content Between Three Genomic Sources
#'
#' @param extract_regions_mt DNAStringSet(List) for Mitochondrial regions.
#' @param extract_regions_pt DNAStringSet(List) for Plastid regions.
#' @param extract_regions_nuc DNAStringSet(List) for Nuclear regions.
#' @param alpha Significance level (default: 0.05).
#'
#' @import Biostrings dplyr rstatix ggplot2 ggpubr
#' @export
analyze_GC_three_genomes <- function(
    extract_regions_mt, 
    extract_regions_pt, 
    extract_regions_nuc,
    alpha = 0.05
) {
  
  process_regions <- function(x) {
    if (is.null(x) || length(unlist(x)) == 0) return(numeric(0))
    Biostrings::letterFrequency(unlist(x), "GC", as.prob = TRUE) * 100
  }
  
  gc_mt <- process_regions(extract_regions_mt)
  gc_pt <- process_regions(extract_regions_pt)
  gc_nuc <- process_regions(extract_regions_nuc)
  
  gc_data <- data.frame(
    GC = c(gc_mt, gc_pt, gc_nuc),
    group = factor(
      rep(c("MT", "PT", "NUC"), c(length(gc_mt), length(gc_pt), length(gc_nuc))),
      levels = c("MT", "PT", "NUC")
    )
  )
  
  norm_test <- gc_data %>%
    group_by(group) %>%
    summarise(
      p.value = if(n() >= 3) shapiro.test(GC)$p.value else NA_real_,
      is_normal = p.value > alpha,
      .groups = "drop"
    )
  
  var_test <- gc_data %>% rstatix::levene_test(GC ~ group)
  is_homoscedastic <- var_test$p > alpha
  
  is_normal_all <- all(norm_test$is_normal, na.rm = TRUE)

  if (is_normal_all) {
    if (is_homoscedastic) {
      main_test <- gc_data %>% rstatix::anova_test(GC ~ group)
      post_hoc  <- gc_data %>% rstatix::tukey_hsd(GC ~ group)
      method    <- "ANOVA / Tukey HSD"
    } else {
      main_test <- gc_data %>% rstatix::welch_anova_test(GC ~ group)
      post_hoc  <- gc_data %>% rstatix::games_howell_test(GC ~ group)
      method    <- "Welch ANOVA / Games-Howell"
    }
  } else {
    main_test <- gc_data %>% rstatix::kruskal_test(GC ~ group)
    post_hoc  <- gc_data %>% rstatix::dunn_test(GC ~ group, p.adjust.method = "bonferroni")
    method    <- "Kruskal-Wallis / Dunn's Test"
  }
  
  stat_p_plot <- post_hoc %>% rstatix::add_y_position(fun = "max", step.increase = 0.12)
  
  gc_plot <- ggplot2::ggplot(gc_data, ggplot2::aes(x = group, y = GC, fill = group)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    ggplot2::geom_jitter(ggplot2::aes(color = group), width = 0.15, alpha = 0.4) +
    ggpubr::stat_pvalue_manual(
      stat_p_plot, 
      label = "p.adj.signif", 
      hide.ns = FALSE,
      tip.length = 0.01
    ) +
    ggplot2::scale_fill_manual(values = c("MT" = "#ff7f0e", "PT" = "#2ca02c", "NUC" = "#1f77b4")) +
    ggplot2::scale_color_manual(values = c("MT" = "#ff7f0e", "PT" = "#2ca02c", "NUC" = "#1f77b4")) +
    ggplot2::labs(
      title = "GC Content Comparison Across Genomes",
      subtitle = paste("Normality:", is_normal_all, "| Homoscedasticity:", is_homoscedastic, "\nMethod:", method),
      y = "GC Content (%)",
      x = "Genomic Source"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  return(list(
    plot = gc_plot,
    main_test = main_test,
    post_hoc = post_hoc,
    method = method,
    normality = norm_test,
    variance = var_test,
    gc_data = gc_data
  ))
}