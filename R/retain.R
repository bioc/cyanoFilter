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
retain <- function(meta_files, decision = c("maxi", "mini")) {


  return(needed)
}
