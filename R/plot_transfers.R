#' Creating a circular plot based on the transfer_function() result
#'
#' @param transfer_function_out Output from transfer_function()
#' @param gap Gap between two neighbour sectors.
#' @param start_degree The starting degree from which the circle begins to draw.
#' @param transparency Transparency of links
#'
#' @return A plot
#' @export
#'
#'
plot_transfers <- function(transfer_function_out = transfer_function_out, gap=40, start_degree=90,transparency=0.7) {

  if (base::missing(transfer_function_out)) {
    stop("The transfer_function_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
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
                               col = rand_color(nrow(links_q),
                                                transparency = transparency),
                               border = NA)
  base::message(base::paste0("Done!"))
}
