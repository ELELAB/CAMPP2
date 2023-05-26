#' @title Create a SummarizedExperiment object
#' @author Nikola Tom
#' @description This function creates a SummarizedExperiment object from given feature counts and sample metadata.
#' It is crucial to ensure that the row names and column names match between the feature counts and sample metadata.
#' If they do not match, the function will stop and throw an error.
#' @param feature_counts A data frame with feature counts. Each column represents a sample,
#' and each row represents a feature (eg. gene) Column names must match the row names in sample metadata.
#' @param sample_metadata A data frame with sample metadata. Each row represents a sample,
#' and each column represents a metadata category. Row names must match the column names in feature counts
#' @return A SummarizedExperiment object containing the given feature_counts and sample_metadata.
#' @examples
#' \dontrun{
#'   rownames(campp2_brca_1_meta)<-campp2_brca_1_meta$ID
#'   campp2_brca_se_1 <- create_SE(campp2_brca_1, campp2_brca_1_meta)
#' }
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
create_SE <- function(feature_counts, sample_metadata) {
    # Check if input data frames have the same number of samples
    if (ncol(feature_counts) != nrow(sample_metadata)) {
        stop("The number of columns in feature_counts must match the number of rows in sample_metadata.")
    }

    # Check if row names in sample_metadata match column names in feature_counts
    if (!all(rownames(sample_metadata) %in% colnames(feature_counts))) {
        stop("All row names in sample_metadata must be present in the column names of feature_counts.")
    }

    # Create the SummarizedExperiment object
    se <- SummarizedExperiment(assays = list(counts = feature_counts),
                               colData = sample_metadata)

    return(se)
}
