#' Analyze and Compare GC Content Between Two Sets of DNA Sequences
#'
#' @param extract_regions_1 A DNAStringSet or DNAStringSetList (e.g., from extract_regions).
#' @param extract_regions_2 A DNAStringSet or DNAStringSetList.
#' @param group1_name Character string, name for the first group (default: "Group 1").
#' @param group2_name Character string, name for the second group (default: "Group 2").
#' @param alpha Significance level for statistical tests (default: 0.05).
#'
#' @return A list with plot, test results, method, normality check, and raw data.
#' @export

analyze_GC_content <- function(
  extract_regions_1, 
  extract_regions_2, 
  group1_name = "Group 1", 
  group2_name = "Group 2",
  alpha = 0.05
) {

  prepare_input <- function(x) {
    if (inherits(x, "DNAStringSetList")) return(unlist(x))
    return(x)
  }

  regions_1_unlisted <- prepare_input(extract_regions_1)
  regions_2_unlisted <- prepare_input(extract_regions_2)
  
  if(length(regions_1_unlisted) == 0 || length(regions_2_unlisted) == 0) {
    stop("One of the input sequence sets is empty.")
  }
  
  calculate_gc <- function(x) Biostrings::letterFrequency(x, "GC", as.prob = TRUE) * 100
  
  gc_data <- data.frame(
    GC = c(as.numeric(calculate_gc(regions_1_unlisted)), 
           as.numeric(calculate_gc(regions_2_unlisted))),
    group = factor(
      rep(c(group1_name, group2_name), 
          c(length(regions_1_unlisted), length(regions_2_unlisted))),
      levels = c(group1_name, group2_name)
    )
  )
  
  select_test <- function(data) {
    norm_test <- data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        p.value = if(dplyr::n() >= 3) shapiro.test(GC)$p.value else 0,
        .groups = "drop"
      )
    
    if(all(norm_test$p.value > alpha, na.rm = TRUE)) {
      var_test <- car::leveneTest(GC ~ group, data = data)
      if(var_test$`Pr(>F)`[1] > alpha) {
        test <- rstatix::t_test(GC ~ group, data = data, var.equal = TRUE)
        method <- "Student's t-test"
      } else {
        test <- rstatix::t_test(GC ~ group, data = data, var.equal = FALSE)
        method <- "Welch's t-test"
      }
    } else {
      test <- rstatix::wilcox_test(GC ~ group, data = data)
      method <- "Wilcoxon rank-sum test"
    }
    list(test = test, method = method, norm_test = norm_test)
  }
  
  test_result <- select_test(gc_data)
  
  test_data <- test_result$test %>%
    mutate(method = test_result$method)
  
  pval_position <- test_data %>% 
    rstatix::add_xy_position(x = "group", dodge = 0.8)
  
  gc_plot <- ggplot(gc_data, aes(x = group, y = GC, fill = group)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, aes(color = group), size = 2) +
    ggpubr::stat_pvalue_manual(
      pval_position,
      label = "{method}\np = {p}",
      tip.length = 0.01,
      bracket.nudge.y = 2,
      inherit.aes = FALSE
    ) +
    labs(
      title = "GC Content Comparison",
      subtitle = test_result$method,
      y = "GC Content (%)",
      x = "Genomic Source"
    ) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  list(
    plot = gc_plot,
    test_result = test_result$test,
    test_method = test_result$method,
    normality_check = test_result$norm_test,
    gc_data = gc_data
  )
}