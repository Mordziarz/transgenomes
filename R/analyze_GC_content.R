#' Analyze and Compare GC Content Between Two Sets of DNA Sequences
#'
#' This function calculates the GC content for two sets of DNA sequences
#' (provided as DNAStringSetList or DNAStringSet objects, e.g., from Biostrings),
#' assigns user-defined group names, performs an appropriate statistical test
#' (Shapiro-Wilk for normality, Levene's test for variance, then t-test, Welch's t-test, or Wilcoxon rank-sum test as appropriate),
#' and generates a publication-ready boxplot with jittered points and p-value annotation.
#'
#' @param extract_regions_1 A DNAStringSetList or DNAStringSet containing the first group of sequences.
#' @param extract_regions_2 A DNAStringSetList or DNAStringSet containing the second group of sequences.
#' @param group1_name Character string, name for the first group (default: "Group 1").
#' @param group2_name Character string, name for the second group (default: "Group 2").
#' @param alpha Significance level for statistical tests (default: 0.05).
#'
#' @return A list with the following elements:
#'   \item{plot}{A ggplot2 object: boxplot with jittered points and p-value annotation.}
#'   \item{test_result}{A data.frame with the result of the statistical test.}
#'   \item{test_method}{A character string: the name of the statistical test used.}
#'   \item{normality_check}{A data.frame with results of the Shapiro-Wilk normality test for each group.}
#'   \item{gc_data}{A data.frame with GC content and group labels for each sequence.}
#'
#' @examples
#' \dontrun{
#' results <- analyze_GC_content(
#'   extract_regions_1 = my_plastome_regions,
#'   extract_regions_2 = my_mitogenome_regions,
#'   group1_name = "Plastome",
#'   group2_name = "Mitogenome"
#' )
#' print(results$plot)
#' results$test_result
#' }
#'
#' @import Biostrings
#' @import dplyr
#' @import rstatix
#' @import ggplot2
#' @import car
#' @export
analyze_GC_content <- function(
  extract_regions_1, 
  extract_regions_2, 
  group1_name = "Group 1", 
  group2_name = "Group 2",
  alpha = 0.05
) {
  regions_1_unlisted <- unlist(extract_regions_1)
  regions_2_unlisted <- unlist(extract_regions_2)
  
  calculate_gc <- function(x) Biostrings::letterFrequency(x, "GC", as.prob = TRUE) * 100
  
  gc_data <- data.frame(
    GC = c(calculate_gc(regions_1_unlisted), calculate_gc(regions_2_unlisted)),
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
        p.value = if(dplyr::n() >= 3) shapiro.test(GC)$p.value else NA_real_,
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
      label = "method\np = {p}",
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
