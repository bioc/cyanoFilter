# Decides if a file should be retiained or removed based on its status
# (from the good_measurement() function) as well as the number of particles
# it contains
#
# @param x dataframe result from good_measurement() function
# @param y dataframe result from good_measurement() function
# @param dil_x dilution level of x
# @param dil_y dilution level of y
#
# @return dataframe with status, dilution, filename of the retained file (either x or y)
#
# @examples
#
# \dontrun{
#   retain(goodFiles_200, goodFiles_500, "200", "500")
# }
#
# @export retain
#
retain <- function(meta_files, make_decision = c("maxi", "mini")) {
      are_both_good <- sum(meta_files$Status == "good")

      if(make_decision == "mini") {

        if(are_both_good == 1) {

          decision <- ifelse(meta_files$Status == "good", "Retain", "No!")

        } else if(are_both_good >= 2) {

          decision <- NULL
          decision[which(meta_files$CellspML == min(meta_files$CellspML))] <- "Retain"
          decision[which(meta_files$CellspML != min(meta_files$CellspML))] <- "No!"

        } else {

          decision <- "All Files are bad!"

        }

      } else {

        if(are_both_good == 1) {

          decision <- ifelse(meta_files$Status == "good", "Retain", "No!")

        } else if(are_both_good >= 2) {

          decision <- NULL
          decision[which(meta_files$CellspML == max(meta_files$CellspML))] <- "Retain"
          decision[which(meta_files$CellspML != max(meta_files$CellspML))] <- "No!"

        } else {

          decision <- "All Files are bad!"

        }

      }


  return(decision)
}
