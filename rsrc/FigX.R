################################################################################
#
# Figure X: sweep through v
# Last updated: 11 Feb 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# LOAD DATA
#######################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

# parameters
beta    <- 0.001 # fixed, for now
gens    <- 20000000
saveplots <- 1
threshold <- 1 # 0 = use all data, 1 = threshold data by min(COUNT)
# 2 = use separate threshold for A-D and E/F
# p       <- 0.
# vs      <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.5)
Mmax    <- 5

# load data
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
pattern   <- sprintf( "vsweep" )
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- data.frame() # initialize data frame

for (i in 1:length(file_list)){
  temp_data    <- read.csv( paste0(file_dir, file_list[i]), header = TRUE)
  temp_data$id <- paste0("run",i) # add id to identify data source
  
  # TEMPORARY
  # if == 0:  ignore data sets where there is at least one line with 0's (out_of_memory)
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

# select rows with specific M, v, beta values
simdata <- simdata[(simdata$β == beta),]
simdata <- simdata[(simdata$M != 0) & (simdata$M <= Mmax),]

# because the number of simulations is uneven at the moment,
# count the minimum number of simulations per case
# so that every case has the same number of simulations
casecount <- simdata %>% 
  group_by(M, K, v, p1, β) %>% 
  summarize(COUNT = n())

if(threshold == 0){
  threshdata_all <- simdata %>% group_by(M, K, v, p1,β) 
}else if(threshold == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1, β) %>% 
    slice_head( n = min(casecount$COUNT) ) 
}

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

############################################################################################
# PREP DATA
############################################################################################
# prep sim data
simdata <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols <- c("id","M","K","u","v","p","beta","epsilon","CC","CD","DC","DD",
                 "CC_final","CD_final","DC_final","DD_final","topinion_mean","topinion_var",
                 "cooperation_all", "cooperation_in", "cooperation_out"#, # cooperation
                 # "opn_simpson_mean","sets_simpson_mean", # simpsons indices
                 # "cityblock_all", "cityblock_in", "cityblock_out", # city block distances
                 # "hamming_all",   "hamming_in",   "hamming_out"    # hamming distances
)
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_strat  <- c("CC","CD","DC","DD")
measure_strat2 <- c("CC_final","CD_final","DC_final","DD_final")
measure_coop   <- c("cooperation_all", "cooperation_in", "cooperation_out")
# measure_pol    <- c("set_1_mean", "topinion_mean","topinion_var",
#                     "opn_simpson_mean","sets_simpson_mean",
#                     "cooperation_all", "cooperation_in", "cooperation_out")
# measure_dist   <- c("cityblock_all", "cityblock_in", "cityblock_out",
#                     "hamming_all",   "hamming_in",   "hamming_out"   )

# select columns
simdata_bymeasure_all <- threshdata_all[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)

# compute average value per strategy for all simulation data
simdata_strat_all <-
  simdata_bymeasure_all[which(simdata_bymeasure_all$variable %in% measure_strat),] %>%
  rename( Fraction = value, Strategy = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Strategy ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

simdata_coop_all <- 
  simdata_bymeasure_all[which(simdata_bymeasure_all$variable %in% measure_coop),] %>%
  rename( Fraction = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

###########################################
# Add extra labels #
###########################################
simdata_coop_all  <- simdata_coop_all  %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )
simdata_strat_all <- simdata_strat_all %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )

###########################################
# Plotting functions
###########################################
plot_figSXcoop <- function(simdata_coop, p, tag = "A", legend = TRUE){
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$p == p & 
                          simdata_coop$M <= Mmax, ]
  
  subdata <- subdata %>% mutate( MK2 = paste0(M2, ", ", K2) )
  
  figSXcoop<- ggplot(data = subdata,
                  aes(x = v, y = Mean, color = MK2, fill = MK2)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text(size = 8),
           legend.key.width = unit(0.02, "npc"),
           legend.key.height = unit(0.4, "cm"),
           panel.spacing = unit(0, "lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    ggtitle( paste0("p = ", p) ) +
    scale_color_manual( values = viridis(7)[1:6], name = "") +
    scale_fill_manual( values = viridis(7)[1:6], name = "") +
    labs(x = "Set mutation rate (v)",
         y = "Effective cooperation",
         tag = tag) + 
    # geom_hline(yintercept = 0.5, color = "gray80") + 
    scale_y_continuous(limits = c(0.2, 0.6), # c(0.47, 0.51),
                       breaks = seq(0.2, 0.6, 0.1)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans  = 'log10') +
    # geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.07, color = NA) +
    geom_line(size = 0.4, alpha = 0.7, lty = 1) +
    geom_point(size = 1.6, alpha = 0.9, stroke = 0.0) #, shape = 1)
  
  return(figSXcoop)
}

plot_figSXstrat <- function(simdata_strat, p, M = 1, K = 1, tag = "B", labeled = TRUE, wlegend = TRUE){
  
  subdata <- simdata_strat[simdata_strat$p == p & 
                           simdata_strat$M == M &
                           simdata_strat$K == K, ]
  
  figSXstrat <- ggplot(subdata,
                  aes(x = v, y = Mean, color = Strategy, group = Strategy, fill = Strategy)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           # legend.key.size = unit(0.025, "npc"),
           legend.key.width = unit(0.02, "npc"),
           legend.key.height = unit(0.4, "cm"),           
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    labs(x = "Opinion mutation rate (v)",
         y = "Relative abundance",
         tag = tag) +
    ggtitle( paste0("M = ", M, ", K = ", K, ", p = ", p) ) +
    scale_color_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_fill_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_y_continuous(limits = c(0.17, 0.33),
                       breaks = seq(0.01, 0.53, 0.04)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans  = 'log10') +
    # geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.4) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.07, color = NA) +
    geom_line(size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.6, alpha = 1, stroke = 0.0) #, shape = 1)
  
  return(figSXstrat)
}


###########################################
# Figure X
###########################################
# save multiplot
figSXa <- plot_figSXcoop( simdata_coop_all, 0.0, "A")
figSXb <- plot_figSXcoop( simdata_coop_all, 0.25, "B")
figSXc <- plot_figSXcoop( simdata_coop_all, 0.5, "C")
figSXd <- plot_figSXcoop( simdata_coop_all, 0.75, "D")
figSXe <- plot_figSXcoop( simdata_coop_all, 1.0, "E")
figSXf <- plot_figSXstrat( simdata_strat_all, 0.0, 1, 1, "F")
figSXg <- plot_figSXstrat( simdata_strat_all, 0.25, 1, 1,  "G")
figSXh <- plot_figSXstrat( simdata_strat_all, 0.5, 1, 1, "H")
figSXi <- plot_figSXstrat( simdata_strat_all, 0.75, 1, 1, "I ")
figSXj <- plot_figSXstrat( simdata_strat_all, 1.0, 1, 1, "J")
emp <- ggplot() + theme_void()

if(saveplots == 1){
  
  # SI Figure with all combos
  threshcount <- if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figSX_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD.png"), 
      width = figW*1.5, height = figW*ratio*1.5/2*5, units = "in", res = 600)
  # multiplot(figSXa, figSXb, figSXc, figSXd, figSXe, emp,
  #           layout = matrix(c(1,2,3,4,5,6), ncol = 3, byrow = TRUE))
  multiplot(figSXa, figSXb, figSXc, figSXd, figSXe, 
            figSXf, figSXg, figSXh, figSXi, figSXj, 
            layout = matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 2, byrow = FALSE))
  dev.off()

  
}

###########################################
# Figure 4
###########################################
# MS Figure with p = 0, 0.5, 1.0
fig4a <- plot_figSXcoop( simdata_coop_all, 0.0, "A")
fig4b <- plot_figSXcoop( simdata_coop_all, 0.5, "B")
fig4c <- plot_figSXcoop( simdata_coop_all, 1.0, "C")
fig4d <- plot_figSXstrat( simdata_strat_all, 0.0, 1, 1, "D")
fig4e <- plot_figSXstrat( simdata_strat_all, 0.5, 1, 1, "E")
fig4f <- plot_figSXstrat( simdata_strat_all, 1.0, 1, 1, "F")
emp <- ggplot() + theme_void()

if(saveplots == 1){
  
  # MS Figure with p = 0, 0.5, 1.0
  threshcount <- if(threshold %in% c(0,1)){ 
    paste0("thresh_", threshold, "_", min(casecount$COUNT))
  }
  plottype <- paste0("fig4_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_SD.png"), 
      width = figW*1.5, height = figW*ratio*1.5/2*3, units = "in", res = 600)
  # multiplot(fig4a, fig4b, fig4c, emp,
  #           layout = matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
  multiplot(fig4a, fig4b, fig4c, fig4d, fig4e, fig4f,
            layout = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = FALSE))
  dev.off()
}


