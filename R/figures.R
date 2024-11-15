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

#' Make NMDS plots.
#'
#' The functions takes a list of colors and the number of times you want those colors to be re-used and
#' makes a larger pallete containing distinct colors.
#'  
#' @param taxa Dataframe containing the taxa of interest as rows and the abundances per sample as 
#' numeric columns. Additional character columns with taxonomic information can be left in as is 
#' the standard output of frac_abund().
#' @param num_color_ls Each number in this numeric list should correspond to the number times you want 
#' to iterate each color in color_ls.
#' @importFrom magrittr %>%
#' @importFrom vegan vegdist
#' @return The output will be a vector containing unique colors.
#' @export
plot_nmds <- function(taxa, groups, dis_ind = "bray", p_color, p_shape){
  mat <- as.matrix(t(select_if(taxa, is.numeric)))
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species", "bin")
  nm <- tail(intersect(ranks, colnames(taxa)[colnames(taxa) %in% ranks]), 1)
  colnames(mat) <- taxa[ ,nm]
  dmat <- vegan::vegdist(mat, method = dis_ind)
  nmds <- vegan::metaMDS(dmat, k = 2, maxit = 1000, trymax = 500)
  sites <- as.data.frame(vegan::scores(nmds, display = "sites"))
  
  if (missing(groups)){
    sites$sample <- rownames(sites)
    ggplot(sites, aes(x=NMDS1,y=NMDS2))+
      geom_point(size=4)+
      theme_bw()+
      ggtitle(paste0(stringr::str_to_title(dis_ind), 
                     ' (dis)-similarity, Euclidean Distance NMDS with stress of ',
                     nmds$stress))+
      ggrepel::geom_label_repel(aes(x=NMDS1,y=NMDS2, label = sample), 
                       box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', 
                       inherit.aes = F)
  } else 
  {
    nm <- colnames(groups)[1]
    sites[ ,nm] <- rownames(sites)
    sites <- merge(sites, groups, by = nm)
    sites %>% ggplot(aes(x=NMDS1,y=NMDS2, color =  .data[[p_color]],
                         shape = .data[[p_shape]]))+
      geom_point(size=4)+
      scale_color_brewer(palette='Set1')+
      theme_bw()+
      ggtitle(paste0(stringr::str_to_title(dis_ind), 
                     ' (dis)-similarity, Euclidean Distance NMDS with stress of ',
                     nmds$stress))+ 
      ggrepel::geom_label_repel(aes(x=NMDS1,y=NMDS2, label = .data[[nm]]), 
                                box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', 
                                inherit.aes = F)
  }

}