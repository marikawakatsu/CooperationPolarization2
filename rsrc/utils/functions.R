##################################################
#
# Various functions
#
################################################### 
require(reshape2)
require(igraph)
require(ggplot2)
require(msm)
require(tidyverse)
require(dplyr) # may need to revert
require(gtools)
require(R.matlab)
require(ggnetwork)
require(ggrepel)
require(viridis)
require(RColorBrewer)
require(varhandle)
require(data.table)
require(grid)
require(gtable)
require(gridExtra)
require(lemon)
require(ggpubr)
require(ggrepel)
require(stargazer)

#########################
# Multiple plot function
#########################
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# from https://github.com/christokita/mixing-model

# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

cool_warm <- function(n) {
  colormap <- Rgnuplot:::GpdivergingColormap(seq(0,1,length.out=n),
                                             rgb1 = colorspace::sRGB( 0.230, 0.299, 0.754),
                                             rgb2 = colorspace::sRGB( 0.706, 0.016, 0.150),
                                             outColorspace = "sRGB")
  colormap[colormap>1] <- 1 # sometimes values are slightly larger than 1
  colormap <- grDevices::rgb(colormap[,1], colormap[,2], colormap[,3])
  colormap
}

# for data cleaning
relabel_cols <- function(simdata){
  
  colnames(simdata)[which(names(simdata) == "β")] <- "beta"
  colnames(simdata)[which(names(simdata) == "ϵ")] <- "epsilon"
  colnames(simdata)[which(names(simdata) == "p1")] <- "p"
  colnames(simdata)[which(names(simdata) == "CC_mean")]  <- "CC"
  colnames(simdata)[which(names(simdata) == "CD_mean")]  <- "CD"
  colnames(simdata)[which(names(simdata) == "DC_mean")]  <- "DC"
  colnames(simdata)[which(names(simdata) == "DD_mean")]  <- "DD"
  colnames(simdata)[which(names(simdata) == "CC_mean_final")]  <- "CC_final"
  colnames(simdata)[which(names(simdata) == "CD_mean_final")]  <- "CD_final"
  colnames(simdata)[which(names(simdata) == "DC_mean_final")]  <- "DC_final"
  colnames(simdata)[which(names(simdata) == "DD_mean_final")]  <- "DD_final"
  colnames(simdata)[which(names(simdata) == "total_means_mean")] <- "topinion_mean"
  colnames(simdata)[which(names(simdata) == "total_vars_mean")]  <- "topinion_var"
  colnames(simdata)[which(names(simdata) == "total_cooperation_mean")] <- "cooperation_all"
  colnames(simdata)[which(names(simdata) == "party_cooperation_mean")] <- "cooperation_in"
  colnames(simdata)[which(names(simdata) == "enemy_cooperation_var")] <- "cooperation_out" # FIX LATER
  colnames(simdata)[which(names(simdata) == "cityblock_ind_mean")]   <- "cityblock_all"
  colnames(simdata)[which(names(simdata) == "cityblock_party_mean")] <- "cityblock_in"
  colnames(simdata)[which(names(simdata) == "cityblock_enemy_mean")] <- "cityblock_out"
  colnames(simdata)[which(names(simdata) == "hamming_ind_mean")]     <- "hamming_all"
  colnames(simdata)[which(names(simdata) == "hamming_party_mean")]   <- "hamming_in"
  colnames(simdata)[which(names(simdata) == "hamming_enemy_mean")]   <- "hamming_out"
  
  return(simdata)
}