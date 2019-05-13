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
      are_all_good <- sum(meta_files$Status == "good")

      decision <- rep(NA, nrow(meta_files))
      decision[which(meta_files$Status != "good")] <- "No!"
      good_pos <- which(meta_files$Status == "good")

      if(make_decision == "mini") {

        if(are_all_good == 1) {

          decision <- ifelse(meta_files$Status == "good", "Retain", "No!")

        } else if(are_all_good >= 2) {

          mini_pos <- which(meta_files$CellspML[good_pos] == min(meta_files$CellspML[good_pos]) )
          nmini_pos <- which(meta_files$CellspML[good_pos] != min(meta_files$CellspML[good_pos]) )

          decision[good_pos[mini_pos]] <- "Retain"

          decision[good_pos[nmini_pos]] <- "No!"


        } else {

          decision <- "All Files are bad!"

        }

      } else if(make_decision == "maxi") {

        if(are_all_good == 1) {

          decision <- ifelse(meta_files$Status == "good", "Retain", "No!")

        } else if(are_all_good >= 2) {

          maxi_pos <- which(meta_files$CellspML[good_pos] == max(meta_files$CellspML[good_pos]) )
          nmaxi_pos <- which(meta_files$CellspML[good_pos] != max(meta_files$CellspML[good_pos]) )

          decision[good_pos[maxi_pos]] <- "Retain"

          decision[good_pos[nmaxi_pos]] <- "No!"


        } else {

          decision <- "All Files are bad!"

        }

      } else stop("Supply make_decision")


   return(decision)
}
