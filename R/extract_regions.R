#' Extract Genomic Regions from BLAST Results
#' 
#' Extracts DNA sequences from a FASTA file based on alignment coordinates 
#' from BLAST results, handling reverse orientations and ID mismatches.
#' 
#' @param transfer_function_out A data.frame/tibble containing BLAST results with 
#'                     sequence IDs and coordinates from transfer_function
#' @param fasta_path Path to the input FASTA file
#' @param id_col Column name in `transfer_function_out` containing sequence IDs 
#'               (default: "subject_id")
#' @param start_col Column name for alignment start positions 
#'                 (default: "s_start")
#' @param end_col Column name for alignment end positions 
#'               (default: "s_end")
#' 
#' @return A `DNAStringSetList` containing extracted sequences with metadata-rich names
#' 
#' @details
#' Key features:
#' - Automatic handling of reverse orientation regions
#' - ID format normalization for FASTA headers and BLAST results
#' - Efficient batch processing using Biostrings/IRanges
#' - Comprehensive error checking and warnings
#' 
#' @export
#' @importFrom Biostrings readDNAStringSet reverseComplement
#' @importFrom IRanges IRanges IRangesList
#' @examples
#' \dontrun{
#' # Example usage:
#' blast_results <-  transfer_function(
#'   fasta_q = "query_genome.fna",
#'   fasta_s = "subject_genome.fna",
#'   bed_q = "query_genes.bed",
#'   bed_s = "subject_genes.bed",
#'   evalue_cut_off = 1e-5,
#'   cores=4
#' )
#' fragments <- extract_regions(
#'   transfer_function_out = blast_results,
#'   fasta_path = "genome_assembly.fna",
#'   id_col = "subject_acc",
#'   start_col = "s_start",
#'   end_col = "s_end"
#' )
#' Biostrings::writeXStringSet(unlist(fragments), "output.fasta")
#' }
extract_regions <- function(
    transfer_function_out, 
    fasta_path, 
    id_col = "subject_id",
    start_col = "s_start",
    end_col = "s_end"
) {

  sequences <- Biostrings::readDNAStringSet(fasta_path)
  names(sequences) <- sapply(strsplit(names(sequences), "\\s+"), `[`, 1)

  blast_ids <- sapply(strsplit(transfer_function_out[[id_col]], "\\|"), `[`, 1)
  
  missing_seqs <- setdiff(blast_ids, names(sequences))
  if(length(missing_seqs) > 0) {
    warning("Missing sequences in FASTA: ", paste(missing_seqs, collapse = ", "))
    transfer_function_out <- transfer_function_out[blast_ids %in% names(sequences), ]
  }
  
  ranges_list <- IRanges::IRangesList(
    lapply(1:nrow(transfer_function_out), function(i) {
      start <- transfer_function_out[[start_col]][i]
      end <- transfer_function_out[[end_col]][i]
      IRanges::IRanges(
        start = min(start, end), 
        end = max(start, end)
      )
    })
  )
  
  fragments <- Biostrings::extractAt(
    sequences[blast_ids], 
    ranges_list
  )
  
  result <- Biostrings::DNAStringSetList(lapply(seq_along(fragments), function(i) {
    if(transfer_function_out[[start_col]][i] > transfer_function_out[[end_col]][i]) {
      Biostrings::reverseComplement(fragments[[i]])
    } else {
      fragments[[i]]
    }
  }))
  
  names(result) <- paste0(
    blast_ids, "|",
    pmin(transfer_function_out[[start_col]], transfer_function_out[[end_col]]), "-",
    pmax(transfer_function_out[[start_col]], transfer_function_out[[end_col]])
  )
  
  return(result)
}
