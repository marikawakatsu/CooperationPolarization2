################################################################################
#
# Figure Z: checking calculations against simulations
# Last updated: 22 Feb 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

# data parameters
epsilon <- 1
u       <- 0.001 # fixed, for now
gens    <- 20000000
saveplots <- 1
threshold <- 1 # 0 = use all data, 1 = threshold data by min(COUNT)
# vs      <- c(0.001, 0.005, 0.025) 

# plotting parameters
ymax      <- 0.3    # max y for plotting
ymin      <- 0.2
yinc      <- 0.02   # y-axis increments
figW      <- 4     # figure width for printing
ratio     <- 2/3   # ratio of fig height to width

# axis ticks etc
qmax   <- 1
qmin   <- 0
qinc   <- 0.25

########################
# LOAD SIMULATION DATA 
########################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

# load data
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
pattern   <- sprintf( "neutral" ) # select files sweeping across M and K
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- data.frame() # initialize data frame

for (i in 1:length(file_list)){
  temp_data    <- read.csv( paste0(file_dir, file_list[i]), header = TRUE)
  temp_data$id <- paste0("run",i) # add id to identify data source
  
  # TEMPORARY WHILE RUNNING SIMULATIONS
  # if == 0:  ignore data sets where there is at least one line with 0's 
  #           (ie consider only complete data sets at a given time)
  # if != -1: include all simulation results (later thresholded by min(COUNT))
  if ( dim(temp_data[temp_data$N == 0,])[1] == 0 ){
    # select common columns
    if (i == 1){
      simdata    <- rbind(simdata, temp_data) #bind the new data to data
    }else if (i > 1){
      commoncols <- intersect(colnames(temp_data), colnames(simdata))
      simdata    <- rbind(simdata[,commoncols], temp_data[,commoncols]) #bind new data
    }
  }
}

#########################
# LOAD CALCULATION DATA 
#########################
calcdata <- read.csv( "data/calc_data/calc_data.csv", header = TRUE)
# calcdata <- read.csv( "data/calc_data/calc_data_finiteN_p_0.5.csv", header = TRUE)

#########################
# PREP DATA
#########################
# select rows with specific M, v, beta values
simdata <- simdata[simdata$v != 0,] # ignore sims that have not been run

# because the number of simulations is uneven at the moment,
# count the minimum number of simulations per case
# so that every case has the same number of simulations
casecount <- simdata %>% 
  group_by(M, K, v, p1, β) %>% 
  summarize(COUNT = n())

if(threshold == 0){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) 
}else if(threshold == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) %>% 
    slice_head( n = min(casecount$COUNT) ) 
}

# prep sim data
simdata        <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols   <- c("id","M","K","u","v","p","beta","epsilon",
                   "CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","sia_sid","sia_sjd")
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_vals  <- c("CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","sia_sid","sia_sjd")

# select columns
simdata_bymeasure_all <- threshdata_all[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)

# compute average value per strategy for all simulation data
simdata_plot <- 
  simdata_bymeasure_all[which(simdata_bymeasure_all$variable %in% measure_vals),] %>%
  rename( Fraction = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

# group calcdata
calcdata <- as.data.frame(calcdata)
calcdata_plot <- calcdata %>% 
  gather("variable","value",-u,-v) %>%
  rename( Value = value, Metric = variable ) %>%
  group_by(v)

###########################################
# Plotting functions
###########################################
plot_figZa <- function(simdata_plot, Metric, tag = "", wtitle = FALSE, Title = ""){
  
  simsubdata   <- simdata_plot[simdata_plot$Metric == Metric, ] # select metric to plot
  simsubdata$v <- factor(simsubdata$v)
  
  figZa <- ggplot(simsubdata,
                  aes(x = p, y = Mean, color = v, group = v)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc"),
           axis.text = element_text (size = 8),
           axis.title = element_text (size = 10),
    ) +
    labs(x = "Party bias (p)",
         y = "Value",
         tag = tag) +
    ggtitle( paste0("quantity: ", if(wtitle){Title}else{Metric} ) ) +
    scale_color_manual(values = rev(viridis(9)) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(-0.05, 1)) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.4, shape = 1)
  
  return(figZa)
}

# v on x-axis
plot_figZb <- function(simdata_plot, calcdata_plot, Metric, tag = "", wtitle = FALSE, Title = ""){
  
  simsubdata    <- simdata_plot[simdata_plot$Metric == Metric, ] # select metric to plot
  simsubdata$p  <- factor(simsubdata$p)
  
  calcsubdata   <- calcdata_plot[calcdata_plot$Metric == Metric, ]
  calcsubdata$u <- factor(calcsubdata$u) # dummy variable

  figZb <- ggplot(simsubdata,
                  aes(x = v, y = Mean, color = p, group = p)) +  
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc"),
           axis.text = element_text (size = 8),
           axis.title = element_text (size = 10),
    ) +
    labs(x = "Issue / opinion exploration (v)",
         y = "Value",
         tag = tag) +
    ggtitle( paste0("quantity: ", if(wtitle){Title}else{Metric} ) ) +
    scale_color_manual(values = rev(viridis(5)) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(-0.05, 1)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans = 'log10') +
    # plot calculation data
    geom_line(data = calcsubdata, 
              aes(x = v, y = Value, group = u, linetype = u),
              alpha = 0.9, size = 0.4, color = "orange") +
    scale_linetype(labels = c("theoretical\nprediction")) +
    guides(linetype = guide_legend("")) +
    # plot simulation data
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.4, shape = 1)
  
  return(figZb)
}

# without calcdata, for E and F
plot_figZc <- function(simdata_plot, Metric, tag = "", wtitle = FALSE, Title = ""){
  
  simsubdata    <- simdata_plot[simdata_plot$Metric == Metric, ] # select metric to plot
  simsubdata$p  <- factor(simsubdata$p)
  
  figZb <- ggplot(simsubdata,
                  aes(x = v, y = Mean, color = p, group = p)) +  
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc"),
           axis.text = element_text (size = 8),
           axis.title = element_text (size = 10),
    ) +
    labs(x = "Issue / opinion exploration (v)",
         y = "Value",
         tag = tag) +
    ggtitle( paste0("quantity: ", if(wtitle){Title}else{Metric} ) ) +
    scale_color_manual(values = rev(viridis(5)) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(-0.05, 1)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans = 'log10') +
    # plot simulation data
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.4, shape = 1)
  
  return(figZb)
}


###########################################
# Figure SZ
###########################################
# # save multiplot
# figZa <- plot_figZa(simdata_plot, "y", "A")
# figZb <- plot_figZa(simdata_plot, "z", "B")
# figZc <- plot_figZa(simdata_plot, "g", "C")
# figZd <- plot_figZa(simdata_plot, "h", "D")
# figZe <- plot_figZa(simdata_plot, "sia_sid", "E", TRUE, "<s_ia s_id>")
# figZf <- plot_figZa(simdata_plot, "sia_sjd", "F", TRUE, "<s_ia s_jd>")
# 
# if(saveplots == 1){
#   
#   threshcount <- 
#     if(threshold %in% c(0,1)){ 
#       paste0("thresh_", threshold, "_", min(casecount$COUNT))
#     }
#   
#   plottype <- paste0("figSZ_", threshcount)
#   
#   png(filename = paste0("plots/figs/", plottype, "_", 
#                         format(Sys.Date(), format="%y%m%d"), "_old.png"), 
#       width = figW*1.5, height = figW*ratio*2.7, units = "in", res = 600)
#   multiplot(figZa, figZb, figZc, figZd, figZe, figZf,
#             layout = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = TRUE))
#   dev.off()
#   
# }

###########################################
# Figure SZ, v2 with data
###########################################
# save multiplot
figZa <- plot_figZb(simdata_plot, calcdata_plot, "y", "A")
figZb <- plot_figZb(simdata_plot, calcdata_plot,"z", "B")
figZc <- plot_figZb(simdata_plot, calcdata_plot,"g", "C")
figZd <- plot_figZb(simdata_plot, calcdata_plot,"h", "D")
figZe <- plot_figZc(simdata_plot, "sia_sid", "E", TRUE, "<s_ia s_id>")
figZf <- plot_figZc(simdata_plot, "sia_sjd", "F", TRUE, "<s_ia s_jd>")

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), ".png"), 
      width = figW*1.5, height = figW*ratio*2.7, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, figZe, figZf,
            layout = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = TRUE))
  dev.off()
  
}

# version without E/F
if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_v2.png"), 
      width = figW*1.5, height = figW*ratio*1.8, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, # figZe, figZf,
            layout = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
  dev.off()
  
}



