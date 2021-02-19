################################################################################
#
# Figure Z
# Last updated: 18 Feb 2021
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
beta    <- 0.001 # fixed, for now
gens    <- 20000000
saveplots <- 1
threshold <- 0 # 0 = use all data, 1 = threshold data by min(COUNT)
# p       <- 0.
vs      <- c(0.001, 0.005, 0.025) #, 0.1)
Mmax    <- 3

# load data
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
pattern   <- sprintf( "run_multi" )
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
simdata <- simdata[(simdata$v %in% vs) & (simdata$β == beta),]
# simdata <- simdata[(simdata$M != 0) & (simdata$M <= Mmax),]

# because the number of simulations is uneven at the moment,
# count the minimum number of simulations per case
# so that every case has the same number of simulations
casecount <- simdata %>% 
  group_by(M, K, v, p1, β) %>% 
  summarize(COUNT = n())

threshdata_all <- simdata %>% group_by(M, K, v, p1) 

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

# select columns
simdata_bymeasure_all <- threshdata_all[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)

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
simdata_coop_all <- simdata_coop_all %>% 
  mutate( M2 = paste0( "M=", M ), 
          K2 = paste0( "K=", K ),
          M3 = paste0( M ),
          K3 = paste0( K )
  )
simdata_strat_all <- simdata_strat_all %>% 
  mutate( M2 = paste0( "M=", M ), 
          K2 = paste0( "K=", K ),
          M3 = paste0( M ),
          K3 = paste0( K )
  )

###########################################
# Plotting functions
###########################################
plot_fig2b <- function(simdata_strat, v, tag = "B", labeled = TRUE, wlegend = TRUE){
  
  subdata <- simdata_strat[simdata_strat$v == v &
                           simdata_strat$M <= Mmax, ]
  
  fig2b <- ggplot(subdata,
                  aes(x = p, y = Mean - 0.25, label = K2, color = Strategy)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
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
         y = "Relative abundance",
         tag = tag) +
    ggtitle( paste0("v = ", v) ) +
    # scale_color_manual(values = rev(viridis(5)[1:4]), # magma(6)[2:5]) +
    scale_color_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020"),
    ) +
    # scale_y_continuous(limits = c(0.16, 0.34),
    #                    breaks = seq(0.16, 0.34, 0.02)) +
    # scale_x_continuous(limits = c(qmin, qmax), 
    #                    breaks = seq(qmin, qmax, qinc)) +
    geom_errorbar(aes(ymin = Mean - 0.25 - 1.96*SE, ymax = Mean - 0.25 + 1.96*SE), width = 0) +
    geom_line(aes(group = Strategy), size = 0.4, alpha = 1, lty = 1) +
    geom_point(stat="identity", size = 1.3, alpha = 1, stroke = 0.5, shape = 1)
  
  return(fig2b)
}

###########################################
# Figure 2
###########################################
# save multiplot
p  <- 0.
v1 <- 0.001
v2 <- 0.005 # 0.025
fig2c <- plot_fig2b(simdata_strat_all, v2, "C", TRUE, TRUE)

if(saveplots == 1){
  
  threshcount <- 
    if(threshold %in% c(0,1)){ 
      # paste0("thresh_", threshold, "_", round( min(casecount$COUNT), digits = -2 ))
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }else if(threshold == 2){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT), "_", min(casecount2$COUNT)) 
    }
  
  plottype <- paste0("fig2_p_", p, "_", threshcount)
  
  png(filename = paste0("plots/figs/", plottype, "_", 
                        format(Sys.Date(), format="%y%m%d"), "_alt.png"), 
      width = figW*1.75*1.65, height = figW*ratio*1.8, units = "in", res = 300)
  multiplot(fig2a, fig2b, fig2c, fig2d, fig2e, fig2f,
            layout = matrix(c(1,1,1,3,3,3,5,5,2,2,2,4,4,4,6,6), ncol = 8, byrow = TRUE))
  # layout = matrix(c(1,2,3,4,5,6), ncol = 3, byrow = FALSE))
  dev.off()
  
}

