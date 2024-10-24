#' Generate a nice color pallette.
#'
#' The functions takes a list of colors and the number of times you want those colors to be re-used and
#' makes a larger pallete containing distinct colors.
#'  
#' @param color_ls A character list of the R colors you want to serve as the basis for your pallete.
#' These should be pretty distinct colors.
#' @param num_color_ls Each number in this numeric list should correspond to the number times you want 
#' to iterate each color in color_ls.
#' #' @return The output will be a vector containing unique colors.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarise_if
#' @export
makecolors <- function(color_ls, num_color_ls){
  palette <- c()
  for (x in 1:length(color_ls)){
    base <- color_ls[x]
    amount <- num_color_ls[x]
    colfunc <- colorRampPalette(c("black", base,"white"))
    palette <- append(palette, colfunc(amount+2)[2:(amount+1)])
  }
  return(palette)
}