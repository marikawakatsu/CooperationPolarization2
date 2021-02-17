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
threshold <- 0 # 0 = use all data, 1 = threshold data by min(COUNT)
# 2 = use separate threshold for A-D and E/F
# p       <- 0.
# vs      <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.5)
Mmax    <- 2

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
  if ( dim(temp_data[temp_data$N == 0,])[1] != -1 ){
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
  threshdata_p0  <- threshdata_all
}else if(threshold == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1, β) %>% 
    # slice_head( n = round( min(casecount$COUNT), digits = -2) )
    slice_head( n = min(casecount$COUNT) ) 
  threshdata_p0  <- threshdata_all
}else if(threshold == 2){
  threshdata_all <- simdata %>% group_by(M, K, v, p1, β) %>% slice_head( n = min(casecount$COUNT) )
  casecount2    <- simdata[simdata$p == p,] %>% group_by(M, K, v, p1, β) %>% summarize(COUNT = n())
  threshdata_p0 <- simdata %>% group_by(M, K, v, p1, β) %>% slice_head( n = min(casecount2$COUNT) )
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
threshdata_p0  <- relabel_cols(threshdata_p0)

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

simdata_bymeasure_p0 <- threshdata_p0[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_p0$value = as.numeric(simdata_bymeasure_p0$value)

# compute average value per strategy for all simulation data
# simdata_strat <- 
#   simdata_bymeasure[which(simdata_bymeasure$variable %in% measure_strat),] %>%
#   rename( Fraction = value, Strategy = variable ) %>%
#   group_by( M, K, u, v, p, beta, epsilon, Strategy ) %>%
#   summarise(
#     Mean = mean(Fraction),
#     SD = sd(Fraction),
#     SE = sd(Fraction) / sqrt(length(Fraction)),
#     numCases = length(Fraction)
#   )

simdata_strat_p0 <-
  simdata_bymeasure_p0[which(simdata_bymeasure_p0$variable %in% measure_strat),] %>%
  rename( Fraction = value, Strategy = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Strategy ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

simdata_coop_p0 <- 
  simdata_bymeasure_p0[which(simdata_bymeasure_p0$variable %in% measure_coop),] %>%
  rename( Fraction = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
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

# simdata_strat2 <- 
#   simdata_bymeasure[which(simdata_bymeasure$variable %in% measure_strat2),] %>%
#   rename( Fraction = value, Strategy = variable ) %>%
#   group_by( M, K, u, v, p, beta, epsilon Strategy ) %>%
#   summarise(
#     Mean = mean(Fraction),
#     SD = sd(Fraction),
#     SE = sd(Fraction) / sqrt(length(Fraction)),
#     numCases = length(Fraction)
#   )

###########################################
# Add extra labels #
###########################################
simdata_coop_p0  <- simdata_coop_p0  %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )
simdata_coop_all <- simdata_coop_all %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )
simdata_strat_p0 <- simdata_strat_p0 %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )

###########################################
# Plotting functions
###########################################
plot_figSXe <- function(simdata_coop, p, tag = "E", legend = TRUE){
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$p == p & 
                          simdata_coop$M <= Mmax, ]
  
  subdata <- subdata %>% mutate( MK2 = paste0(M2, ", ", K2) )
  
  figSXe <- ggplot(data = subdata,
                  aes(x = v, y = Mean, color = MK2, fill = MK2)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text(size = 8),
           legend.key.width = unit(0.02, "npc"),
           legend.key.height = unit(0.04, "npc"),
           panel.spacing = unit(0, "lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    ggtitle( paste0("p = ", p) ) +
    scale_color_manual( values = rev(viridis(6)), name = "") +
    scale_fill_manual( values = rev(viridis(6)), name = "") +
    labs(x = "Set mutation rate (v)",
         y = "Effective cooperation",
         tag = tag) + 
    # geom_hline(yintercept = 0.5, color = "gray80") + 
    scale_y_continuous(limits = c(0.2, 0.6), # c(0.47, 0.51),
                       breaks = seq(0.2, 0.6, 0.1)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans  = 'log10') +
    # geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.1, color = NA) +
    geom_line(size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1., alpha = 1, stroke = 0.5)
  
  return(figSXe)
}

###########################################
# Figure X
###########################################
# save multiplot
figSXa <- plot_figSXe( simdata_coop_all, 0.0, "A")
figSXb <- plot_figSXe( simdata_coop_all, 0.25, "B")
figSXc <- plot_figSXe( simdata_coop_all, 0.5, "C")
figSXd <- plot_figSXe( simdata_coop_all, 0.75, "D")
figSXe <- plot_figSXe( simdata_coop_all, 1.0, "E")
emp <- ggplot() + theme_void()

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      # paste0("thresh_", threshold, "_", round( min(casecount$COUNT), digits = -2 ))
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }else if(threshold == 2){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT), "_", min(casecount2$COUNT)) 
    }
  
  plottype <- paste0("figSX_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), ".png"), 
      width = figW*2.25, height = figW*ratio*1.5, units = "in", res = 300)
  multiplot(figSXa, figSXb, figSXc, figSXd, figSXe, emp,
            layout = matrix(c(1,2,3,4,5, 6), ncol = 3, byrow = TRUE))
  dev.off()
  
}


