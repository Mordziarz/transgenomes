#' Create a Circular Visualization for Genomic Alignment Data
#'
#' Generates a Circos plot to visualize genomic alignments between query and subject sequences using output from `transfer_function()`.
#'
#' @param transfer_function_out A data frame from [transfer_function()] containing alignment data. Must include columns: `query_id`, `subject_id`, `q_start`, `q_end`, `s_start`, `s_end`, `q_len`, `s_len`.
#' @param gap Numeric (30-100). Space between genome tracks in degrees (default: 40).
#' @param start_degree Numeric (0-360). Initial rotation angle for the circular plot (default: 90).
#' @param transparency Numeric (0-1). Opacity of alignment ribbons (0 = fully transparent, 1 = opaque, default: 0.7).
#' @param color Character. Color of the links (ribbons) connecting aligned regions. Can be any color name recognized by R (e.g., "green", "red") or a vector of colors (default: "green").

#' @return Circos plot visualizing genomic alignments
#' @export
#'
#' @examples
#' \dontrun{
#' # After running transfer_function():
#' plot_transfers(alignment_results, gap = 45, transparency = 0.5)
#' }
#' @importFrom circlize circos.par circos.initializeWithIdeogram circos.genomicLink

plot_transfers <- function(transfer_function_out = transfer_function_out, gap=40, start_degree=90,transparency=0.7,color="green") {

  if (missing(transfer_function_out)) {
    stop("Required argument 'transfer_function_out' is missing. Provide output from transfer_function().", call. = FALSE)
  }
  
  required_cols <- c("query_id", "subject_id", "q_start", "q_end", "s_start", "s_end", "q_len", "s_len")
  if (!all(required_cols %in% colnames(transfer_function_out))) {
    stop("Input data missing required columns. Verify transfer_function() output structure.", call. = FALSE)
  }

  blast_n <- base::as.data.frame(transfer_function_out)
  links_s <- transfer_function_out[,c("subject_id","s_start","s_end")]
  links_q <- transfer_function_out[,c("query_id","q_start","q_end")]

  genome_plot <- base::data.frame("V1" = c(transfer_function_out$query_id[1],
                                           transfer_function_out$subject_id[1]),
                                  "V2" = c(0,
                                           0),
                                  "V3" = c(transfer_function_out$q_len[1],
                                           transfer_function_out$s_len[1]),
                                  "V4" = c(transfer_function_out$query_id[1],
                                           transfer_function_out$subject_id[1]),
                                  "V5" = c(transfer_function_out$query_id[1],
                                           transfer_function_out$subject_id[1]))

  circlize::circos.par(start.degree = start_degree,"gap.degree" = c(gap,gap))
  circlize::circos.initializeWithIdeogram(genome_plot)
  circlize::circos.genomicLink(links_q,
                               links_s,
                               col = "green",
                               border = NA)
  base::message(base::paste0("Done!"))
}
