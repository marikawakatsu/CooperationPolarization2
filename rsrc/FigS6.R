################################################################################
#
# Figure S6: checking calculations against simulations, with M = 1
# Last updated: 26 Jul 2021
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
# path to your directory
setwd("~/.julia/dev/CooperationPolarization2/")

# load data
file_dir <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
simdata  <- read.csv( paste0(file_dir, "della_neutral_merged_M1.csv"), header = TRUE)

# check that every case has the same number of simulations
casecount <- simdata %>% group_by(M, K, v, p1, β) %>% summarize(COUNT = n())

if( length(unique(casecount$COUNT)) == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) 
  print( sprintf("each parameter setting has %d runs", unique(casecount$COUNT)) )
}else{
  stop("!!! check casecount !!!")
}

#########################
# LOAD CALCULATION DATA 
#########################
calcdata   <- read.csv( "analytics/calc_data_yzgh_all.csv", header = TRUE)

# calcdatap1  <- read.csv( "analytics/calc_data_y_p1.csv", header = TRUE)
# calcdata$y1 <- calcdatap1$y # add p=1 data for y only

#########################
# PREP DATA
#########################
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
calcdata_plot <- calcdata[calcdata$M == 1,] %>% 
  gather("variable","value",-u,-v) %>%
  rename( Value = value, Metric = variable ) %>%
  group_by(v)

###########################################
# Plotting functions
###########################################
# plot simulated data with predictions, with v on x-axis
plot_figZb <- function(simdata_plot, calcdata_plot, Metric, tag = "", wtitle = FALSE, Title = ""){
  
  if(Metric == "yy"){
    # for y, plot two prediction lines + p <= 1 data
    calcsubdata   <- calcdata_plot[calcdata_plot$Metric %in% c("y","y1"), ]
    simsubdata    <- simdata_plot[simdata_plot$Metric == Metric, ]
  }else{
    # for g, h, z, plot one prediction line + p < 1 data only
    calcsubdata   <- calcdata_plot[calcdata_plot$Metric == Metric, ]
    simsubdata    <- simdata_plot[simdata_plot$Metric == Metric & simdata_plot$p < 1, ]
  }
  simsubdata$p  <- factor(simsubdata$p)
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
    scale_color_manual(values = rev(viridis(7)[2:6]) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(-0.05, 1)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans = 'log10') +
    # plot calculation data
    geom_line(data = calcsubdata, 
              aes(x = v, y = Value, group = Metric, linetype = Metric),
              alpha = 0.9, size = 0.4, color = "orange") +
    scale_linetype(labels = c("theoretical\nprediction\nfor p<1", 
                              "theoretical\nprediction\nfor p=1")) +
    guides(linetype = guide_legend("")) +
    # plot simulation data
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.4, shape = 1)
  
  return(figZb)
}

###########################################
# Figure SZ, v2 with data
###########################################
# save multiplot
figZa <- plot_figZb(simdata_plot, calcdata_plot, "y", "A")
figZb <- plot_figZb(simdata_plot, calcdata_plot,"z", "B")
figZc <- plot_figZb(simdata_plot, calcdata_plot,"g", "C")
figZd <- plot_figZb(simdata_plot, calcdata_plot,"h", "D")

# version without E/F
if(saveplots == 1){
  
  threshcount <- paste0("thresh_", threshold, "_", min(casecount$COUNT))
  plottype    <- paste0("figSZ_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_v2.png"), 
      width = figW*1.5, height = figW*ratio*1.8, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, # figZe, figZf,
            layout = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
  dev.off()
  
}



