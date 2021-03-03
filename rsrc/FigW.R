################################################################################
#
# Figure W: checking calculations against simulations
# Last updated: 1 Mar 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# LOAD DATA
#######################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

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
pattern   <- sprintf( "vsweep_M1" )
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

#########################
# PREP DATA
#########################
# prep sim data
simdata <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols <- c("id","M","K","u","v","p","beta","epsilon","CC","CD","DC","DD",
                 "CC_final","CD_final","DC_final","DD_final",
                 "cooperation_all", "cooperation_in", "cooperation_out"#, # cooperation
)
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_strat  <- c("CC","CD","DC","DD")
measure_strat2 <- c("CC_final","CD_final","DC_final","DD_final")
measure_coop   <- c("cooperation_all", "cooperation_in", "cooperation_out")

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

simdata_coop_all  <- simdata_coop_all  %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )
simdata_plot <- simdata_strat_all %>% mutate( M2 = paste0( "M=", M ), K2 = paste0( "K=", K ) )

###########################################
# LOAD CALC DATA
###########################################
calcdata <- read.csv( "analytics/calc_data_strat.csv", header = TRUE)

# group calcdata
calcdata <- as.data.frame(calcdata)
calcdata_plot <- calcdata %>% 
  gather("variable","value",-u,-v) %>%
  rename( Value = value, Metric = variable ) %>%
  group_by(v)

############################################################################################
# PLOTTING FUNCTIONS
############################################################################################
plot_figSWstrat <- function(simdata_plot, calcldata_plot, p, tag = "B", labeled = TRUE, wlegend = TRUE){
  
  subdata <- simdata_plot[simdata_plot$p == p, ]
  
  calcsubdata   <- calcdata_plot
  calcsubdata$u <- factor(calcdata_plot$u) # dummy variable
  
  figSWstrat <- ggplot(subdata,
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
    labs(x = "Issue / opinion exploration (v)",
         y = "Relative abundance",
         tag = tag) +
    ggtitle( paste0("p = ", p) ) +
    scale_color_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_fill_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_y_continuous(limits = c(0.163, 0.337),
                       breaks = seq(0.01, 0.53, 0.04)) +
    scale_x_continuous(limits = c(0.001, 0.625),
                       breaks = c(0.001, 0.005, 0.025, 0.125, 0.625),
                       trans  = 'log10') +
    # plot calculation data
    geom_line(data = calcsubdata, aes(x = v, y = Value, linetype = u),
              alpha = 0.9, size = 0.4) +
    scale_linetype(labels = c("theoretical\nprediction")) +
    guides(linetype = guide_legend("")) +
    # plot simulation data
    # geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0., size = 0.3) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.1, color = NA) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.3, alpha = 0.9, stroke = 0.4, shape = 21)
  
  return(figSWstrat)
}

