#' @title Fix zeros
#' @description This function detects zero/negative values and replaces them
#' in a SummarizedExperiment object. In the first step, it checks for the presence
#' of zero and negative values. Next, features with zero counts surpassing the
#' threshold set by the size of smallest sample group are removed. Subsequently,
#' the remaining zeros are substituted with the minimal positive values observed
#' for each feature across all samples.
#' @param SE a SummarizedExperiment object of gene/abundance counts.
#' @param group_colname a string specifying the column name in the colData(SE)
#' representing the sample group.
#' @param remove.sparse.features a logical argument (TRUE/FALSE) for removal of
#' features with a sum of zero counts exceeding the size of the smallest
#' sample group and for replacing the remaining zeros. Default is TRUE.
#' @export
#' @import SummarizedExperiment
#' @return a SummarizedExperiment object with adjusted zeros. Features with a sum of zero counts
#' exceeding the smallest sample group size are removed and the remaining zeros are replaced.
#' @examples {
#' # Assuming you have a SummarizedExperiment object and relevant metadata column
#' SE_fixed = FixZeros(SE = some_SE_object, group_colname = "diagnosis")
#' }


FixZeros <- function(SE, group_colname, remove.sparse.features=TRUE) {

    # Extracting the group data from SE
    if(!group_colname %in% colnames(colData(SE))){
        stop("The specified group column name is not found in the SE colData.")
    }
    group <- colData(SE)[[group_colname]]

    if(any(assay(SE) == 0)){
        message("Data includes 0-value(s)")
    } else {
        message("Data doesn't include 0-value(s)")
    }

    if(any(assay(SE) < 0)){
        message("Data includes negative value(s)")
    } else {
        message("Data doesn't include negative value(s)")
    }

    ###Removal of features with high 0-counts
    if(remove.sparse.features) {
        smallestGr <- min(table(group))
        greaterthanBG <- rowSums(assay(SE) > 0)
        lessthanBG <- which(greaterthanBG < smallestGr)

        if (length(lessthanBG) > 0) {
            SE <- SE[-lessthanBG,]
            message(paste0(length(lessthanBG), " features removed due to excessive zeros"))
        }
    }

    ###Replacement of zeros by the lowest value for each feature
    min_per_row <- apply(assay(SE), 1, function(x) min(x[x != 0]))
    zero_positions <- which(assay(SE) == 0, arr.ind=TRUE)
    assay(SE)[zero_positions] <- min_per_row[zero_positions[, 1]]

    return(SE)
}
