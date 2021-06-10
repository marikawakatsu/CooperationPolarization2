################################################################################
#
# Figure 4
# Last updated: 1 May 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# LOAD DATA
#######################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

######################################################
# MS version, showing behavior across values of v
######################################################

# parameters # not changed
beta  <- 0.001
u     <- 0.001 # fixed, for now
gens  <- 20000000
saveplots <- 1
Mmax      <- 5
threshold <- 1

vs     <- c(0.001, 0.005, 0.025, 0.125, 0.625) # for plotting vsweep data

# load data
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
pattern   <- sprintf( "vsweep" ) # specify data type
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- data.frame() # initialize data frame

for (i in 1:length(file_list)){
  temp_data    <- read.csv( paste0(file_dir, file_list[i]), header = TRUE)
  temp_data$id <- paste0("run",i) # add id to identify data source
  
  # TEMPORARY
  # ignore data sets where there is at least one line with 0's (out_of_memory)
  if ( dim(temp_data[temp_data$N == 0,])[1] == 0){
    # select common columns
    if (i == 1){
      simdata    <- rbind(simdata, temp_data) #bind the new data to data
    }else if (i > 1){
      commoncols <- intersect(colnames(temp_data), colnames(simdata))
      simdata    <- rbind(simdata[,commoncols], temp_data[,commoncols]) #bind new data
    }
  }
}

# select vs
simdata   <- simdata[(simdata$v %in% vs) & (simdata$β == beta),]
casecount <- simdata %>% group_by(M, K, v, p1, β) %>% summarize(COUNT = n())
if(threshold == 1){simdata <- simdata %>% group_by(M, K, v, p1, β) %>% slice_head( n = min(casecount$COUNT) ) }

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

# color palettes
cityblock_colors <- c("#7bccc4","#4eb3d3","#2b8cbe","#0868ac","#084081") 
hamming_colors   <- c("#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b") 

############################################################################################
# PREP DATA
############################################################################################
# prep sim data
simdata        <- relabel_cols(simdata)
threshdata_all <- relabel_cols(threshdata_all)

# melt data
select_cols <- c("id","M","K","u","v","p","beta","epsilon","CC","CD","DC","DD",
                 "CC_final","CD_final","DC_final","DD_final","topinion_mean","topinion_var",
                 "cooperation_all", "cooperation_in", "cooperation_out", # cooperation
                 # "opn_simpson_mean","sets_simpson_mean", # simpsons indices
                 "cityblock_all", "cityblock_in", "cityblock_out", # city block distances
                 "hamming_all",   "hamming_in",   "hamming_out"    # hamming distances
)
id_vars       <- c("M","K","u","v","p","beta","epsilon")

# different sets of measurements to plot
measure_strat  <- c("CC","CD","DC","DD")
measure_coop   <- c("cooperation_all", "cooperation_in", "cooperation_out")
measure_pol    <- c("topinion_mean","topinion_var",
                    "opn_simpson_mean","sets_simpson_mean",
                    "topinion_sd", "topinion_sd_normed"
)
measure_strat2 <- c("CC_final","CD_final","DC_final","DD_final")
measure_dist   <- c("cityblock_all", "cityblock_in", "cityblock_out",
                    "hamming_all",   "hamming_in",   "hamming_out",
                    "cityblock_all_normed", "cityblock_in_normed", "cityblock_out_normed",
                    "hamming_all_normed",   "hamming_in_normed",   "hamming_out_normed" 
)

simdata_bymeasure <- simdata[ select_cols ] %>% 
  # normalize avg opinion distance, opinion
  mutate(
    topinion_sd          = sqrt(topinion_var),
    topinion_sd_normed   = sqrt(topinion_var) / (2*M+1),
    cityblock_all_normed = cityblock_all / (2*K),
    cityblock_in_normed  = cityblock_in / (2*K),
    cityblock_out_normed = cityblock_out / (2*K),
    hamming_all_normed   = hamming_all / max(1, ( 2 * abs(M %/% 2 - floor(abs(M/2 - K))) ) ),
    hamming_in_normed    = hamming_in / max(1, ( 2 * abs(M %/% 2 - floor(abs(M/2 - K))) ) ),
    hamming_out_normed   = hamming_out / max(1, ( 2 * abs(M %/% 2 - floor(abs(M/2 - K))) ) )
  ) %>%
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)

simdata_bymeasure$value = as.numeric(simdata_bymeasure$value)

# compute average value per strategy for all simulation data
simdata_pol <- 
  simdata_bymeasure[which(simdata_bymeasure$variable %in% measure_pol),] %>%
  rename( Value = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    SE = sd(Value) / sqrt(length(Value)),
    numCases = length(Value)
  )

simdata_dist <- 
  simdata_bymeasure[which(simdata_bymeasure$variable %in% measure_dist),] %>%
  rename( Value = value, Metric = variable ) %>%
  group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    SE = sd(Value) / sqrt(length(Value)),
    numCases = length(Value)
  )


#####################################
# Add extra labels 
#####################################
simdata_pol <- simdata_pol %>%
  mutate( M2 = paste0( "M=", M),
          K2 = paste0( "K=",K ) ) %>%
  mutate( MK2 = paste0( M2,", ", K2) )


simdata_dist <- simdata_dist %>%
  mutate( M2 = paste0( "M=", M),
          K2 = paste0( "K=", K ) ) %>%
  mutate( MK2 = paste0( M2,", ", K2) )


#####################################
# Plotting functions
#####################################
plotA_normed <- function(simdata_dist, tag = "A", option = "all", colors = cityblock_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  if(option == "one"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_all_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level")
    linetypes <- c("solid")
  }else if(option == "two"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_in_normed", "cityblock_out_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("within-party","between-party")
    linetypes <- c("solid", "dotted","dashed")
  }else if(option == "all"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_all_normed", "cityblock_in_normed", "cityblock_out_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level","within-party","between-party")
    linetypes <- c("solid","dotted","dashed")
  }
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  figSa <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text  = element_text (size = 7),
          legend.title = element_text (size = 8),
          panel.spacing = unit(0.2/3*Mmax,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Partisan bias (p)",
         y = "Average opinion distance",
         tag = tag) +
    scale_linetype_manual(values = linetypes, 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Issue/opinion\nexploration (v)") +
    scale_fill_manual(values = colors, 
                      name   = "Issue/opinion\nexploration (v)")  +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1.0) +
    geom_point(aes(y = Mean), size = 0.8, alpha = 1, stroke = 0.4) + #, shape = 1) + ) +
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 11),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.key.size = unit(0.03, "npc"),
          legend.position = c(0.70, 0.97), 
          legend.justification = c(0,1))
  
  figSa_grob <- ggplotGrob(figSa)
  figSa_grob <- gtable_filter(
    figSa_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  figSa <- as_ggplot(figSa_grob)
  
  return(figSa)
}

plotB_normed <- function(simdata_dist, tag = "B", option = "all", colors = hamming_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  if(option == "one"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_all_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level")
    linetypes <- c("solid")
  }else if(option == "two"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_in_normed", "hamming_out_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("within-party","between-party")
    linetypes <- c("solid", "dotted","dashed")
  }else if(option == "all"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_all_normed", "hamming_in_normed", "hamming_out_normed") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level","within-party","between-party")
    linetypes <- c("solid","dotted","dashed")
  }
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  figSb <- ggplot(data = subdata, 
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 7),
          legend.title = element_text (size = 8),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2/3*Mmax,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Partisan bias (p)",
         y = "Average issue distance",
         tag = tag) + 
    scale_linetype_manual(values = linetypes, 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Issue/opinion\nexploration (v)" ) +
    scale_fill_manual(values = colors,
                      name   = "Issue/opinion\nexploration (v)" )  +
    scale_y_continuous(limits = c(0, 0.8),
                       breaks = seq(0, 0.8, 0.2)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1.0) +
    geom_point(aes(y = Mean), size = 0.8, alpha = 1, stroke = 0.4) + #, shape = 1) + ) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 11),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.key.size = unit(0.03, "npc"),
          legend.position = c(0.70, 0.97), 
          legend.justification = c(0,1))
  
  figSb_grob <- ggplotGrob(figSb)
  figSb_grob <- gtable_filter(
    figSb_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  figSb <- as_ggplot(figSb_grob)
  
  return(figSb)
}

# without normalization
plotA_unnormed <- function(simdata_dist, tag = "A", option = "all", colors = cityblock_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  if(option == "one"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_all") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level")
    linetypes <- c("solid")
  }else if(option == "two"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_in", "cityblock_out") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("within-party","between-party")
    linetypes <- c("solid", "dotted","dashed")
  }else if(option == "all"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("cityblock_all", "cityblock_in", "cityblock_out") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level","within-party","between-party")
    linetypes <- c("solid","dotted","dashed")
  }
  
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  fig3b <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 7),
          legend.title = element_text (size = 8),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2/3*Mmax,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Partisan bias (p)",
         y = "Average opinion distance",
         tag = tag) +
    scale_linetype_manual(values = linetypes, 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Issue/opinion\nexploration (v)") +
    scale_fill_manual(values = colors, 
                      name   = "Issue/opinion\nexploration (v)")  +
    scale_y_continuous(limits = c(0, 6),
                       breaks = seq(0, 10, 1)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1.0) +
    geom_point(aes(y = Mean), size = 0.8, alpha = 1, stroke = 0.4) + #, shape = 1) + ) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 11),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.key.size = unit(0.03, "npc"),
          legend.position = c(0.70, 0.97), 
          legend.justification = c(0,1))
  
  fig3b_grob <- ggplotGrob(fig3b)
  fig3b_grob <- gtable_filter(
    fig3b_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  fig3b <- as_ggplot(fig3b_grob)
  
  return(fig3b)
}

plotB_unnormed <- function(simdata_dist, tag = "B", option = "all", colors = hamming_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  if(option == "one"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_all") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level")
    linetypes <- c("solid")
  }else if(option == "two"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_in", "hamming_out") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("within-party","between-party")
    linetypes <- c("solid", "dotted","dashed")
  }else if(option == "all"){
    subdata   <- simdata_dist[simdata_dist$Metric %in% c("hamming_all", "hamming_in", "hamming_out") & 
                                simdata_dist$M <= Mmax, ]
    labels    <- c("population-level","within-party","between-party")
    linetypes <- c("solid","dotted","dashed")
  }
  
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  fig4a <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 7),
          legend.title = element_text (size = 8),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2/3*Mmax,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Partisan bias (p)",
         y = "Average issue distance",
         tag = tag) + 
    scale_linetype_manual(values = linetypes, 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Issue/opinion\nexploration (v)" ) +
    scale_fill_manual(values = colors,
                      name   = "Issue/opinion\nexploration (v)" )  +
    scale_y_continuous(limits = c(0, 2.4),
                       breaks = seq(0, 2.4, 0.4)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    # geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1.0) +
    geom_point(aes(y = Mean), size = 0.8, alpha = 1, stroke = 0.4) + #, shape = 1) + ) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 11),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.key.size = unit(0.03, "npc"),
          legend.position = c(0.70, 0.97), 
          legend.justification = c(0,1))
  
  fig4a_grob <- ggplotGrob(fig4a)
  fig4a_grob <- gtable_filter(
    fig4a_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  fig4a <- as_ggplot(fig4a_grob)
  
  return(fig4a)
}


###########################################
# MS Figure with vsweep
###########################################
figSa <- plotA_normed(simdata_dist, "A", "two")
figSb <- plotB_normed(simdata_dist, "B", "two")

if(saveplots == 1){
  
  plottype <- paste0("fig4_normed")
  png(filename = paste0("plots/figs/", plottype, "_vsweep_",
                        format(Sys.Date(), format="%y%m%d"), "_", min(casecount$COUNT), "_SD.png"), # !!! change !!! 
      width = figW*1.9/1.5, height = figW*ratio*1.75*2, units = "in", res = 600)
  multiplot(figSa, figSb, cols = 1)
  dev.off()
  
}

###########################################
# SI Figure with vsweep, unnormalized
###########################################
figSa_unnormed <- plotA_unnormed(simdata_dist, "A", "two")
figSb_unnormed <- plotB_unnormed(simdata_dist, "B", "two")

if(saveplots == 1){
  
  plottype <- paste0("figSC_unnormed")
  png(filename = paste0("plots/figs/", plottype, "_vsweep_",
                        format(Sys.Date(), format="%y%m%d"), "_", min(casecount$COUNT), "_SD.png"), # !!! change !!! 
      width = figW*1.9/1.5, height = figW*ratio*1.75*2, units = "in", res = 600)
  multiplot(figSa_unnormed, figSb_unnormed, cols = 1)
  dev.off()
  
}

###########################################
# SI Figure with vsweep, population-level only
###########################################
figSa_normed <- plotA_normed(simdata_dist, "A", "one")
figSb_normed <- plotB_normed(simdata_dist, "B", "one")
# figSa_unnormed <- plotA_unnormed(simdata_dist, "C", "one")
# figSb_unnormed <- plotB_unnormed(simdata_dist, "D", "one")

if(saveplots == 1){
  
  plottype <- paste0("figSC_poplevel")
  png(filename = paste0("plots/figs/", plottype, "_vsweep_",
                        format(Sys.Date(), format="%y%m%d"), "_", min(casecount$COUNT), "_SD.png"), # !!! change !!! 
      width = figW*1.9/0.75/2, height = figW*ratio*1.75*2, units = "in", res = 600)
  multiplot(figSa_normed, figSb_normed, 
            # figSa_unnormed, figSb_unnormed, cols = 2)
            cols = 1)
  dev.off()
  
}

