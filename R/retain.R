#' Decides if a file should be retiained or removed based on its status.
#'
#' @param meta_files dataframe from meta file that has been preprocessed by the \code{\link{goodfcs}} function.
#' @param make_decision decision to be made should more than one \eqn{cells/\mu L} be good.
#' @param Status column name in meta_files containing status obtained from the \code{\link{goodfcs}} function.
#' @param CellspML column name in meta_files containing \eqn{cells/\mu L} measurements.
#'
#' @return a character vector with entries "Retain" for a file to be retained or "No!" for a file to be discarded.
#'
#' @description Function to determine what files to retain and finally read from the flow cytometer FCS file.
#'
#' @details It is typically not known in advance which dilution level would result in the desired \eqn{cells/\mu L}, therefore
#'          the samples are ran through the flow cytometer at two or more dilution levels. Out of these, one has to decide which
#'          to retain and finally use for further analysis. This function and \code{\link{goodfcs}} are to help you decide that.
#'          If more than one of the dilution levels are judged good, the option \emph{make_decision = "maxi"} will give "Retain" to the
#'          row with the maximum \eqn{cells/\mu L} while the opposite occurs for \emph{make_decision = "mini"}.
#'
#' @seealso \code{\link{goodfcs}}
#'
#' @examples \dontrun{
#' retain(meta_files = dataframe, make_decision = "maxi")
#' }
#'
#' @export retain
#'
retain <- function(meta_files, make_decision = c("maxi", "mini"), Status = "Status", CellspML = "CellspML") {
    are_all_good <- sum(meta_files[Status] == "good")

    decision <- rep(NA, nrow(meta_files))
    decision[which(meta_files[Status] != "good")] <- "No!"
    good_pos <- which(meta_files[Status] == "good")

    if (make_decision == "mini") {

        if (are_all_good == 1) {

            decision <- ifelse(meta_files[Status] == "good", "Retain", "No!")

        } else if (are_all_good >= 2) {

            mini_pos <- which(meta_files[CellspML][good_pos] == min(meta_files[CellspML][good_pos]))
            nmini_pos <- which(meta_files[CellspML][good_pos] != min(meta_files[CellspML][good_pos]))

            decision[good_pos[mini_pos]] <- "Retain"

            decision[good_pos[nmini_pos]] <- "No!"


        } else {

            decision <- "All Files are bad!"

        }

    } else if (make_decision == "maxi") {

        if (are_all_good == 1) {

            decision <- ifelse(meta_files[Status] == "good", "Retain", "No!")

        } else if (are_all_good >= 2) {

            maxi_pos <- which(meta_files[CellspML][good_pos] == max(meta_files[CellspML][good_pos]))
            nmaxi_pos <- which(meta_files[CellspML][good_pos] != max(meta_files[CellspML][good_pos]))

            decision[good_pos[maxi_pos]] <- "Retain"

            decision[good_pos[nmaxi_pos]] <- "No!"


        } else {

            decision <- "All Files are bad!"

        }

    } else stop("Supply make_decision")


    return(decision)
}
