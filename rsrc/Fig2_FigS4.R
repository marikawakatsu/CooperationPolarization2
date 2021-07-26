################################################################################
#
# Figure 2 and Figure S4
# Last updated: 26 Jul 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# SET PARAMETERS
#######################
# parameters
epsilon <- 1
u       <- 0.001 # fixed, for now
beta    <- 0.001 # fixed, for now
gens    <- 20000000
saveplots <- 1
threshold <- 1 # 0 = use all data, 1 = threshold data by min(COUNT)
vs        <- c(0.001, 0.005, 0.025) 
Mmax      <- 5

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
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
simdata   <- read.csv( paste0(file_dir, "della_all_merged.csv"), header = TRUE)

# select rows with specific M, v, beta values
simdata   <- simdata[(simdata$v %in% vs) & (simdata$β == beta),]
simdata   <- simdata[(simdata$M != 0) & (simdata$M <= Mmax),]

# check that every case has the same number of simulations
casecount <- simdata %>% group_by(M, K, v, p1, β) %>% summarize(COUNT = n())

if( length(unique(casecount$COUNT)) == 1){
  threshdata_all <- simdata %>% group_by(M, K, v, p1) 
  print( sprintf("each parameter setting has %d runs", unique(casecount$COUNT)) )
}else{
  stop("!!! check casecount !!!")
}

############################################################################################
# PREP DATA
############################################################################################
# prep sim data
simdata        <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols <- c("id","M","K","u","v","p","beta","epsilon","CC","CD","DC","DD",
                 "CC_final","CD_final","DC_final","DD_final","topinion_mean","topinion_var",
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
set.seed(42) # for gg_repel

plot_fig2a <- function(simdata_coop, p, tag){
  
  simdata_coop$v <- factor(simdata_coop$v)
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$p == p &
                          simdata_coop$M <= Mmax, ]
  
  fig2a <- ggplot(subdata,
                  aes(x = K2, y = Mean, label = K2, color = v)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           # legend.key.size = unit(0.025, "npc"),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc"),
           axis.text.x = element_text(angle = 60, vjust = 0.5, color = "gray50")
    ) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    labs(x = "",
         y = "Effective cooperation",
         tag = tag) +
    ggtitle( paste0("p = ", p, ", grouped by M") ) +
    scale_color_manual(values = rev(viridis(5)[1:4]), # magma(6)[2:5],
                       name = "Issue/opinion\nexploration (v)") +
    scale_y_continuous(limits = c(0.1, 0.6), 
                       breaks = seq(0.1, 0.6, 0.1)) +
    # geom_hline(yintercept = 0.5, color = "gray80") +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.3) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 1, color = NA) +
    geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(stat="identity", size = 1.1, alpha = 1, stroke = 0.5) + #, shape = 21, fill = "white") +
    facet_grid( ~ M2, space="free_x",
                switch = "x", scales="free_x") +
    theme(strip.placement = "outside") +
    theme(axis.title.x = element_blank(), 
          strip.background = element_blank(),
          strip.text.x.bottom = element_text(angle = 60))
  
  return(fig2a)
}

plot_fig2a_v2 <- function(simdata_coop, p, tag){
  
  simdata_coop$v <- factor(simdata_coop$v)
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                            simdata_coop$p == p &
                            simdata_coop$M <= Mmax, ]
  
  fig2a <- ggplot(subdata,
                  aes(x = M2, y = Mean, label = M2, color = v)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text (size = 8),
           # legend.key.size = unit(0.025, "npc"),
           legend.key.width = unit(0.015, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.2,  "lines"),
           legend.margin = margin(t = 0, unit="npc"),
           axis.text.x = element_text(angle = 60, vjust = 0.5, color = "gray50")
    ) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    labs(x = "",
         y = "Effective cooperation",
         tag = tag) +
    ggtitle( paste0("p = ", p, ", grouped by K") ) +
    scale_color_manual(values = rev(viridis(5)[1:4]), # magma(6)[2:5],
                       name = "Issue/opinion\nexploration (v)") +
    scale_y_continuous(limits = c(0.1, 0.6),  
                       breaks = seq(0.1, 0.6, 0.1)) +
    # geom_hline(yintercept = 0.5, color = "gray80") +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.3) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 1) +
    geom_line(aes(group = v), size = 0.4, alpha = 1, lty = 1) +
    geom_point(stat="identity", size = 1.1, alpha = 1, stroke = 0.5) + #, shape = 21, fill = "white") +
    facet_grid( ~ K2, space="free_x",
                switch = "x", scales="free_x") +
    theme(strip.placement = "outside") +
    theme(axis.title.x = element_blank(), 
          strip.background = element_blank(),
          strip.text.x.bottom = element_text(angle = 60))
  
  return(fig2a)
}

plot_fig2b <- function(simdata_strat, p, v, tag = "B", labeled = TRUE, wlegend = TRUE){
  
  subdata <- simdata_strat[simdata_strat$p == p & 
                           simdata_strat$v == v &
                           simdata_strat$M <= Mmax, ]
  
  ylabel <- if(labeled){"Frequency"}else{""}
  
  fig2b <- ggplot(subdata,
                  aes(x = K2, y = Mean, label = K2, color = Strategy)) +
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
           legend.margin = margin(t = 0, unit="npc"),
           axis.text.x = element_text(angle = 60, vjust = 0.5, color = "gray50")
    ) +
    labs(x = "",
         y = ylabel,
         tag = tag) +
    ggtitle( paste0("p = ", p, ", v = ", v) ) +
    scale_color_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_y_continuous(limits = c(0.0, 0.55),
                       breaks = seq(0.0, 0.5, 0.1)) +
    # scale_x_continuous(limits = c(qmin, qmax), 
    #                    breaks = seq(qmin, qmax, qinc)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.3) +
    geom_line(aes(group = Strategy), size = 0.4, alpha = 1, lty = 1) +
    geom_point(stat="identity", size = 1.1, alpha = 1, stroke = 0.5) + #, shape = 21, fill = "white") + 
    facet_grid( ~ M2, space="free_x",
                switch = "x", scales="free_x") +
    theme(strip.placement = "outside") +
    theme(axis.title.x = element_blank(), 
          strip.background = element_blank(),
          strip.text.x.bottom = element_text(angle = 60))
  
  return(fig2b)
}

# heatmap
plot_fig2e <- function(simdata_coop, v, tag = "E", legend = TRUE){
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$v == v & 
                          simdata_coop$M <= Mmax, ]
  subdata$p <- factor(subdata$p)
  
  fig2e <- ggplot(data = subdata,
                  aes(x = K3, y = p, fill = Mean)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text(size = 8),
           legend.key.size = unit(0.02, "npc"),
           panel.spacing = unit(0.25,"lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    ggtitle( paste0("v = ", v) ) +
    scale_fill_gradient2(low = "#BD1513", mid = "#CFD5D9", high = "#00428B",
                         midpoint = 0.5,
                         limit = c(0.2,0.55),
                         space = "Lab",
                         name = "Effective\ncooperation") +
    labs(x = "Partisan bias (p)",
         y = "",
         tag = tag) +
    geom_tile( show.legend = legend ) + 
    facet_grid(. ~ M2, space="free", 
               switch = "x", scales="free") +
    theme(strip.placement = "outside") +
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank())
  
  return(fig2e)
}

# heatmap, v2
plot_fig2e_v3 <- function(simdata_coop, v, tag = "E", legend = TRUE){
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$v == v & 
                          simdata_coop$M <= Mmax, ] %>% 
    arrange(p) # sort by p
  subdata$pmin <- rep( c(-0.025, 0.125, 0.375, 0.625, 0.775, 0.825, 0.875, 0.925, 0.975, 0.975), each = 15 )
  subdata$pmax <- rep( c( 0.125, 0.375, 0.625, 0.775, 0.825, 0.875, 0.925, 0.975, 0.975, 1.025), each = 15 )
  # subdata$p <- factor(subdata$p)
  labels <- subdata %>% 
    group_by(M, K) %>% 
    summarize( M2 = unique(M2), 
               K2 = unique(K2))
  
  fig2e <- ggplot(data = subdata,
                  aes(x = p, y = M, fill = Mean)) +
    theme_classic() +
    theme(plot.title = element_text(hjust  = 0.5,
                                    size   = 10,  
                                    margin = margin(-10,0,3,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text(size = 8),
           legend.key.width = unit(0.01, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.01, "npc"),
           legend.margin = margin(t=0, r=0, b=0.5, l=0, unit="cm"),
           axis.text.y = element_text(size = 7),
           axis.text.x = element_text(size = 7),
    ) +
    ggtitle( paste0("v = ", v) ) +
    scale_fill_viridis(limit = c(0.2, 0.55),
                       direction = -1,
                       name = "Effective\ncooperation") + # option 3
    labs(x = "Partisan bias (p)",
         y = "",
         tag = tag) +
    # geom_tile( show.legend = legend, 
               # width = rep(c(0.05, 0.25, 0.25, 0.25, 0.05, 0.05, 0.05, 0.05, 0.01, 0.05), 15)) + 
    geom_rect( aes(xmin = pmin, xmax = pmax, ymin = M-0.5, ymax = M+0.5) ) + 
    scale_y_continuous( breaks=1:5, 
                        labels=unique(labels$M2), 
                        expand = c(0, 0)
                        ) +
    scale_x_continuous( limits = c(-0.05, 1.05) ) +
    facet_grid(K2 ~ ., space="free", 
               switch = "y", scales="free_y") +
    theme(strip.placement = "outside") +
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0))
  
  return(fig2e)
}

# lineplot
plot_fig2e_v2 <- function(simdata_coop, v, tag = "E", legend = TRUE){
  
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$v == v & 
                          simdata_coop$M < 4, ]
  
  subdata <- subdata %>% mutate( MK2 = paste0(M2, ", ", K2) )
  # subdata$q <- factor(subdata$q)
  
  fig2e <- ggplot(data = subdata,
                  aes(x = p, y = Mean, color = MK2, fill = MK2)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 10, 
                                    margin = margin(0,0,0,0) ) ) +
    theme (legend.text = element_text (size = 7),
           legend.title = element_text(size = 8),
           legend.key.width = unit(0.01, "npc"),
           legend.key.height = unit(0.03, "npc"),
           panel.spacing = unit(0.25,"lines"),
           legend.margin = margin(t = 0, unit="npc")
    ) +
    ggtitle( paste0("v = ", v) ) +
    scale_color_manual( values = rev(viridis(6)), name = "") +
    scale_fill_manual( values = rev(viridis(6)), name = "") +
    labs(x = "Partisan bias (p)",
         y = "Effective cooperation",
         tag = tag) +
    # geom_hline(yintercept = 0.5, color = "gray80") + 
    scale_y_continuous(limits = c(0.1, 0.6), 
                       breaks = seq(0.1, 0.6, 0.1)) +
    # geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.3) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.1, color = NA) +
    geom_line(size = 0.4, alpha = 1, lty = 1) +
    geom_point(size = 1., alpha = 1, stroke = 0.5)
  
  return(fig2e)
}

# effect size, with M fixed
tableS2_M <- function(simdata_coop, p){
  
  # extract data with specified p
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                          simdata_coop$p == p &
                          simdata_coop$M <= Mmax, ]
  
  effectdata <- subdata %>%
    group_by(M, v) %>%
    summarise(
      minK        = min(K),
      Mean_minK   = Mean[K == min(K)],
      maxK        = max(K),
      Mean_maxK   = Mean[K == max(K)],
      effect_size = (Mean[K == max(K)] - Mean[K == min(K)]) / Mean[K == min(K)]*100
    ) %>%
    arrange(v)
  
  stargazer(as.data.frame(effectdata), summary=FALSE, rownames=FALSE)
  return(effectdata)
}

# effect size, with K fixed
tableS2_K <- function(simdata_coop, p){
  
  # extract data with specified p
  subdata <- simdata_coop[simdata_coop$Metric == "cooperation_all" & 
                            simdata_coop$p == p &
                            simdata_coop$M <= Mmax, ]
  
  effectdata <- subdata %>%
    group_by(K, v) %>%
    summarise(
      minM        = min(M),
      Mean_minM   = Mean[M == min(M)],
      maxM        = max(M),
      Mean_maxM   = Mean[M == max(M)],
      effect_size = (Mean[M == max(M)] - Mean[M == min(M)]) / Mean[M == min(M)]*100,
    ) %>%
    arrange(v)

  return(stargazer(as.data.frame(effectdata), summary=FALSE, rownames=FALSE))
}


###########################################
# Figure 2
###########################################
# save multiplot
p  <- 0.
v1 <- 0.001
v2 <- 0.025
fig2b <- plot_fig2a(   simdata_coop_all,  p, "B")
fig2a <- plot_fig2a_v2(simdata_coop_all,  p, "A")
fig2c <- plot_fig2b(simdata_strat_all, p, v1, "C", TRUE, TRUE)
fig2d <- plot_fig2b(simdata_strat_all, p, v2, "D", TRUE, TRUE)
fig2e <- plot_fig2e_v3( simdata_coop_all, v1, "E")
fig2f <- plot_fig2e_v3( simdata_coop_all, v2, "F")

if(saveplots == 1){
  
  threshcount <- if(threshold %in% c(0,1)){ 
      paste0("thresh_", threshold, "_", min(casecount$COUNT))
    }
  plottype <- paste0("fig2_p_", p, "_", threshcount)
  
  # png(filename = paste0("plots/figs/", plottype, "_", 
  #                       format(Sys.Date(), format="%y%m%d"), "_SD.png"), # !!! change !!!
  pdf(file = paste0("plots/figs/", plottype, "_", 
                    format(Sys.Date(), format="%y%m%d"), "_SD.pdf"), #  PDF
      width = figW*1.75*1.5, height = figW*ratio*1.8)
  multiplot(fig2a, fig2b, fig2c, fig2d, fig2e, fig2f,
            layout = matrix(c(1,1,1,3,3,3,5,5,2,2,2,4,4,4,6,6), ncol = 8, byrow = TRUE))
            # layout = matrix(c(1,1,1,1,3,3,3,5,5,2,2,2,2,4,4,4,6,6), ncol = 9, byrow = TRUE))
  dev.off()
  
} 

###########################################
# Figure S4 with p = 1
###########################################
# save multiplot
p  <- 1.
v1 <- 0.001
v2 <- 0.025
fig2b <- plot_fig2a(   simdata_coop_all,  p, "B")
fig2a <- plot_fig2a_v2(simdata_coop_all,  p, "A")
fig2c <- plot_fig2b(simdata_strat_all, p, v1, "C", TRUE, TRUE)
fig2d <- plot_fig2b(simdata_strat_all, p, v2, "D", TRUE, TRUE)

if(saveplots == 1){

  threshcount <- if(threshold %in% c(0,1)){
    paste0("thresh_", threshold, "_", min(casecount$COUNT))
  }
  plottype <- paste0("figSA_p_", p, "_", threshcount)

  png(filename = paste0("plots/figs/", plottype, "_",
                        format(Sys.Date(), format="%y%m%d"), "_SD.png"), # !!! change !!!
      width = figW*1.75*1.25, height = figW*ratio*1.8, units = "in", res = 600)
  multiplot(fig2a, fig2b, fig2c, fig2d,
            layout = matrix(c(1,3,2,4), ncol = 2, byrow = TRUE))
  dev.off()

}



