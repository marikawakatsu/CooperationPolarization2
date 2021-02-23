################################################################################
#
# Figure Z: checking calculations against simulations
# Last updated: 22 Feb 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# LOAD DATA
#######################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

# parameters
epsilon <- 1
u       <- 0.001 # fixed, for now
gens    <- 20000000
saveplots <- 1
threshold <- 1 # 0 = use all data, 1 = threshold data by min(COUNT)
vs      <- c(0.001, 0.005, 0.025) 

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

# select rows with specific M, v, beta values
simdata <- simdata[simdata$v %in% vs,]

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
simdata        <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols   <- c("id","M","K","u","v","p","beta","epsilon",
                   "CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","sia_sid","sia_sjd"
                   )
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_vals  <- c("CC_final","CD_final","DC_final","DD_final",
                   "y","z","g","h","sia_sid","sia_sjd")

# select columns
simdata_bymeasure_all <- threshdata_all[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)

# compute average value per strategy for all simulation data
simdata_plot_all <- 
  simdata_bymeasure_all[which(simdata_bymeasure_all$variable %in% measure_vals),] %>%
  rename( Fraction = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    SE = sd(Fraction) / sqrt(length(Fraction)),
    numCases = length(Fraction)
  )

###########################################
# Plotting functions
###########################################
plot_figZa <- function(simdata_plot_all, Metric, tag = "", wtitle = FALSE, Title = ""){
  
  subdata   <- simdata_plot_all[simdata_plot_all$Metric == Metric, ] # select metric to plot
  subdata$v <- factor(subdata$v)
  
  figZb <- ggplot(subdata,
                  aes(x = p, y = Mean, color = v, group = v)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           # legend.key.size = unit(0.025, "npc"),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    labs(x = "Party bias (p)",
         y = "Value",
         tag = tag) +
    ggtitle( paste0("Metric: ", if(wtitle){Title}else{Metric} ) ) +
    scale_color_manual(values = rev(viridis(5)[1:4]) ) + # magma(6)[2:5]) +
    scale_y_continuous(limits = c(-0.1, 1)) +
    # scale_x_continuous(limits = c(qmin, qmax), 
    #                    breaks = seq(qmin, qmax, qinc)) +
    geom_errorbar(aes(ymin = Mean - 1.96*SE, ymax = Mean + 1.96*SE), width = 0) +
    # geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1.3, alpha = 1, stroke = 0.5, shape = 1)
  
  return(figZb)
}

###########################################
# Figure Z
###########################################
# save multiplot
figZa <- plot_figZa(simdata_plot_all, "y", "A")
figZb <- plot_figZa(simdata_plot_all, "z", "B")
figZc <- plot_figZa(simdata_plot_all, "g", "C")
figZd <- plot_figZa(simdata_plot_all, "h", "D")
figZe <- plot_figZa(simdata_plot_all, "sia_sid", "E", TRUE, "<s_ia s_id>")
figZf <- plot_figZa(simdata_plot_all, "sia_sjd", "F", TRUE, "<s_ia s_jd>")

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  
  plottype <- paste0("figZ_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), ".png"), 
      width = figW*1.75*1.65, height = figW*ratio*1.8, units = "in", res = 300)
  multiplot(figZa, figZb, figZc, figZd, figZe, figZf,
            layout = matrix(c(1,2,3,4,5,6), ncol = 3, byrow = TRUE))
  dev.off()
  
}

