################################################################################
#
# Figure ZZ: checking calculations against simulations, with M2 and M3
# Last updated: 19 Jun 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

# data parameters
epsilon <- 1
N       <- 40
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
pattern   <- sprintf( "neutral_M2|neutral_M3" ) # select files sweeping across M and K
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
# calcdata   <- read.csv( "analytics/calc_data_yzgh.csv", header = TRUE)
calcdata   <- read.csv( "analytics/calc_data_yzgh_all.csv", header = TRUE)
calcdatap1 <- read.csv( "analytics/calc_data_y_p1.csv", header = TRUE)

calcdata$y1 <- calcdatap1$y # add p=1 data for y only

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
    # slice_head( n = min(casecount$COUNT) ) 
    slice_head( n = 75 ) 
}

# prep sim data
simdata        <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols   <- c("id","M","K","u","v","p","beta","epsilon",
                   "CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","zp","gp","hp","sia_sid","sia_sjd")
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_vals  <- c("CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","zp","gp","hp","sia_sid","sia_sjd")

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
  gather("variable","value",-M,-K,-u,-v) %>%
  rename( Value = value, Metric = variable ) %>%
  group_by(v)

###########################################
# Plotting functions
###########################################
# plot simulated data with predictions, with v on x-axis
plot_figZb <- function(simdata_plot, calcdata_plot, Metric, M = 1, K = 1, tag = "", wtitle = FALSE, Title = ""){
  
  if(Metric == "yy"){ # old
    # for y, plot two prediction lines + p <= 1 data
    calcsubdata   <- calcdata_plot[calcdata_plot$Metric %in% c("y","y1") &
                                   calcdata_plot$M == M &
                                   calcdata_plot$K == K, ]
    simsubdata    <- simdata_plot[simdata_plot$Metric == Metric &
                                  simdata_plot$M == M &
                                  simdata_plot$K == K, ]
  }else{
    # for g, h, z, plot one prediction line + p < 1 data only
    calcsubdata   <- calcdata_plot[calcdata_plot$Metric == Metric &
                                   calcdata_plot$M == M &
                                   calcdata_plot$K == K, ]
    simsubdata    <- simdata_plot[simdata_plot$Metric == Metric & 
                                  simdata_plot$M == M &
                                  simdata_plot$K == K &
                                  simdata_plot$p < 1, ]
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

# without calcdata, for E and F
plot_figZc <- function(simdata_plot, Metric, M = 1, K = 1, tag = "", wtitle = FALSE, Title = ""){
  
  simsubdata    <- simdata_plot[simdata_plot$Metric == Metric & 
                                simdata_plot$M == M &
                                simdata_plot$K == K &
                                simdata_plot$p <= 1, ] # select metric to plot
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
    scale_y_continuous(limits = c(-0.10, 1)) +
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
# Figure SZZ without data
###########################################
# save multiplot
M <- 2
K <- 2
figZa <- plot_figZc(simdata_plot, "y", M, K, "A", FALSE)
figZb <- plot_figZc(simdata_plot, "z", M, K, "B", FALSE)
figZc <- plot_figZc(simdata_plot, "g", M, K, "C", FALSE)
figZd <- plot_figZc(simdata_plot, "h", M, K, "D", FALSE)
figZe <- plot_figZc(simdata_plot, "zp", M, K, "E", FALSE)
figZf <- plot_figZc(simdata_plot, "gp", M, K, "F", FALSE)
figZg <- plot_figZc(simdata_plot, "hp", M, K, "G", FALSE)
figZh <- plot_figZc(simdata_plot, "sia_sid", M, K, "H", FALSE)
figZi <- plot_figZc(simdata_plot, "sia_sjd", M, K, "I", FALSE)


if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount, "_M", M, "K", K)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD.png"), 
      width = figW*1.5*1.5, height = figW*ratio*2.7, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, figZe, figZf, figZg, figZh, figZi,
            layout = matrix(c(1,2,3,4,5,6,7,8,9), ncol = 3, byrow = TRUE))
  dev.off()
  
}

# with data
M <- 3
K <- 2
figZa <- plot_figZb(simdata_plot, calcdata_plot, "y", M, K, "A", FALSE)
figZb <- plot_figZb(simdata_plot, calcdata_plot, "z", M, K, "B", FALSE)
figZc <- plot_figZb(simdata_plot, calcdata_plot, "g", M, K, "C", FALSE)
figZd <- plot_figZb(simdata_plot, calcdata_plot, "h", M, K, "D", FALSE)
figZe <- plot_figZb(simdata_plot, calcdata_plot, "zp", M, K, "E", FALSE)
figZf <- plot_figZb(simdata_plot, calcdata_plot, "gp", M, K, "F", FALSE)
figZg <- plot_figZb(simdata_plot, calcdata_plot, "hp", M, K, "G", FALSE)
figZh <- plot_figZc(simdata_plot, "sia_sid", M, K, "H", FALSE)
figZi <- plot_figZc(simdata_plot, "sia_sjd", M, K, "I", FALSE)
emp <- ggplot() + theme_void()

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount, "_M", M, "K", K)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD+calc.png"), 
      width = figW*1.5, height = figW*ratio*3.2, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, figZe, figZf, figZg, emp,
            layout = matrix(c(1,2,3,4,5,6,7,8), ncol = 2, byrow = TRUE))
  dev.off()
  
}


# with data
M <- 1
K <- 1
figZa <- plot_figZb(simdata_plot, calcdata_plot, "y", M, K, "A", FALSE)
figZb <- plot_figZb(simdata_plot, calcdata_plot, "z", M, K, "B", FALSE)
figZc <- plot_figZb(simdata_plot, calcdata_plot, "g", M, K, "C", FALSE)
figZd <- plot_figZb(simdata_plot, calcdata_plot, "h", M, K, "D", FALSE)

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount, "_M", M, "K", K)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD+calc.png"), 
      width = figW*1.5, height = figW*ratio*2, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, 
            layout = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
  dev.off()
  
}

###########################################
# Figure SZ, v2 with data
###########################################
# # save multiplot
# figZa <- plot_figZb(simdata_plot, calcdata_plot, "y", "A")
# figZb <- plot_figZb(simdata_plot, calcdata_plot,"z", "B")
# figZc <- plot_figZb(simdata_plot, calcdata_plot,"g", "C")
# figZd <- plot_figZb(simdata_plot, calcdata_plot,"h", "D")
# figZe <- plot_figZc(simdata_plot, "sia_sid", "E", TRUE, "<s_ia s_id>")
# figZf <- plot_figZc(simdata_plot, "sia_sjd", "F", TRUE, "<s_ia s_jd>")
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
#                         format(Sys.Date(), format="%y%m%d"), ".png"), 
#       width = figW*1.5, height = figW*ratio*2.7, units = "in", res = 600)
#   multiplot(figZa, figZb, figZc, figZd, figZe, figZf,
#             layout = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = TRUE))
#   dev.off()
#   
# }
# 
# # version without E/F
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
#                         format(Sys.Date(), format="%y%m%d"), "_v2.png"), 
#       width = figW*1.5, height = figW*ratio*1.8, units = "in", res = 600)
#   multiplot(figZa, figZb, figZc, figZd, # figZe, figZf,
#             layout = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
#   dev.off()
#   
# }
# 

#####################################
# TEMPORARY - CHECKING CALCULATIONS #
#####################################
if(threshold == 0){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) 
}else if(threshold == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) %>% 
    slice_head( n = min(casecount$COUNT) ) 
}
threshdata_all <- threshdata_all %>% 
  mutate(
    CC_calc = 1/4 + N*0.001*(1-u)/u * (1/4)*M*( b*(gp-hp) + c*(hp-zp) ),
    CD_calc = 1/4 + N*0.001*(1-u)/u * (1/4)*M*( b*(g-h) + c*(h-z) ),
    DC_calc = 1/4 + N*0.001*(1-u)/u * (1/4)*M*( b*(-g+h) + c*(-h+z) ),
    DD_calc = 1/4 + N*0.001*(1-u)/u * (1/4)*M*( b*(-gp+hp) + c*(-hp+zp) )
  )
threshdata_all <- relabel_cols(threshdata_all)
select_cols2 <- c("id","M","K","u","v","p","beta","epsilon",
                  "CC_calc","CD_calc","DC_calc","DD_calc", # temp
                  "y","z","g","h","zp","gp","hp","sia_sid","sia_sjd")
measure_vals2 <- c("CC_calc","CD_calc","DC_calc","DD_calc")

simdata_bystrat_all <- threshdata_all[ select_cols2 ] %>% 
  gather("strategy","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bystrat_all$value = as.numeric(simdata_bystrat_all$value)

simdata_strat <- 
  simdata_bystrat_all[which(simdata_bystrat_all$strategy %in% measure_vals2),] %>%
  rename( Fraction = value, Metric = strategy ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

# check intermediate result
plot_figZd <- function(simdata_strat, M = 1, K = 1, p = 0, tag = "", Title = ""){
  
  simsubdata    <- simdata_strat[simdata_strat$M == M &
                                   simdata_strat$K == K &
                                   simdata_strat$p == p, ]
  
  figZd <- ggplot(simsubdata,
                  aes(x = v, y = Mean, color = Metric, group = Metric)) +  
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
    ggtitle( paste0("M = ", M, ", K = ", K, ", p = ", p) ) +
    scale_color_manual(values = rep(c("#0571b0","#92c5de","#f4a582","#ca0020")) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(0.13, 0.39),
                       breaks = seq(0.01, 0.53, 0.04)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans = 'log10') +
    # plot calculation data
    # geom_line(data = calcsubdata, 
    #           aes(x = v, y = Value, group = Metric, linetype = Metric),
    #           alpha = 0.9, size = 0.4, color = "orange") +
    scale_linetype(labels = c("theoretical\nprediction\nfor p<1", 
                              "theoretical\nprediction\nfor p=1")) +
    guides(linetype = guide_legend("")) +
    # plot simulation data
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.4, shape = 1)
  
  return(figZd)
}

# save multiplot
figZa <- plot_figZd(simdata_strat, 2, 1, 0, "A", FALSE)
figZb <- plot_figZd(simdata_strat, 2, 2, 0, "B", FALSE)
figZc <- plot_figZd(simdata_strat, 3, 1, 0,  "C", FALSE)
figZd <- plot_figZd(simdata_strat, 3, 2, 0,  "D", FALSE)
figZe <- plot_figZd(simdata_strat, 3, 3, 0,  "E", FALSE)
emp <- ggplot() + theme_void()

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSZ_", threshcount, "_intermediate")
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD.png"), 
      width = figW*1.5, height = figW*ratio*2.7, units = "in", res = 600)
  multiplot(figZa, figZb, figZc, figZd, figZe, emp,
            layout = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = TRUE))
  dev.off()
  
}