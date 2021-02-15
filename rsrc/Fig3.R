################################################################################
#
# Figure 3
# Last updated: 14 Feb 2020
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

#######################
# LOAD DATA
#######################

setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later

# parameters
beta  <- 0.001
u     <- 0.001 # fixed, for now
gens  <- 20000000
saveplots <- 1
Mmax      <- 5

vs     <- c(0.001, 0.025)
vlabel <- if(length(vs) < 3){ paste0("v_",vs[1],"+",vs[2]) }else{ "all" }

# load data
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )
pattern   <- sprintf( "M" )
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- data.frame() # initialize data frame

for (i in 1:length(file_list)){
  temp_data    <- read.csv( paste0(file_dir, file_list[i]), header = TRUE)
  temp_data$id <- paste0("run",i) # add id to identify data source
  
  # TEMPORARY
  # ignore data sets where there is at least one line with 0's (out_of_memory)
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

# for fixed gamma, consider only the rows with gamma == gamma
simdata <- simdata[(simdata$v %in% vs) & (simdata$β == beta),]

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
if(length(vs) == 3){
  cityblock_colors <- rev(heat.colors(7)[c(1,3,5)])
}else if(0.025 %in% vs){
  cityblock_colors <- rev(heat.colors(7)[c(3,5)])
}else if(0.1 %in% vs){
  cityblock_colors <- rev(heat.colors(7)[c(1,5)])
}

if(length(vs) == 3){
  hamming_colors <- rev(magma(7)[3:5])
}else if(0.025 %in% vs){
  hamming_colors <- rev(magma(7)[4:5])
}else if(0.1 %in% vs){
  hamming_colors <- rev(magma(7)[c(3,5)])
}

if(length(vs) == 3){
  topinion_colors <- rev( viridis(6)[2:5] )
}else if(0.025 %in% vs){
  topinion_colors <- rev( viridis(6)[3:4] )
}else if(0.1 %in% vs){
  topinion_colors <- rev( viridis(6)[c(2,4)] )
}

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
# Figure 3A: Total opinion dispersion  #
#####################################
simdata_pol <- simdata_pol %>%
  mutate( M2 = paste0( "M=", M),
          K2 = paste0( "K=",K ) ) %>%
  mutate( MK2 = paste0( M2,", ", K2) )

###########################################
# Figure 3B: 
###########################################
simdata_dist <- simdata_dist %>%
  mutate( M2 = paste0( "M=", M),
          K2 = paste0( "K=", K ) ) %>%
  mutate( MK2 = paste0( M2,", ", K2) )

###########################################
# Plotting functions
###########################################
plot_fig3b <- function(simdata_dist, tag = "B", colors = cityblock_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  subdata <- simdata_dist[simdata_dist$Metric %in% c("cityblock_all", "cityblock_in", "cityblock_out") & 
                            simdata_dist$M <= Mmax, ]
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  labels <- c("population-level","within-party","between-party")
  
  fig3b <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 8),
          legend.title = element_text (size = 9),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Party bias (p)",
         y = "Average opinion distance",
         tag = tag) +
    scale_linetype_manual(values = c("solid","dotted","dashed"), 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Rate of opinion\nmutation (v)") +
    scale_fill_manual(values = colors, 
                      name   = "Rate of opinion\nmutation (v)")  +
    scale_y_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 1)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1) +
    geom_point(aes(y = Mean), size = 1.3, alpha = 1, stroke = 0.5) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 10),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.position = c(0.7, 1), 
          legend.justification = c(0,1))
  
  fig3b_grob <- ggplotGrob(fig3b)
  fig3b_grob <- gtable_filter(
    fig3b_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  fig3b <- as_ggplot(fig3b_grob)
  
  return(fig3b)
}

plot_fig4a <- function(simdata_dist, tag = "A", colors = hamming_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  subdata <- simdata_dist[simdata_dist$Metric %in% c("hamming_all", "hamming_in", "hamming_out") & 
                            simdata_dist$M <= Mmax, ]
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  labels <- c("population-level","within-party","between-party")
  
  fig4a <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 8),
          legend.title = element_text (size = 9),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Party bias (p)",
         y = "Average issue distance",
         tag = tag) + 
    scale_linetype_manual(values = c("solid","dotted","dashed"), 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Rate of opinion\nmutation (v)" ) +
    scale_fill_manual(values = colors,
                      name   = "Rate of opinion\nmutation (v)" )  +
    scale_y_continuous(limits = c(0, 1.5),
                       breaks = seq(0, 1.5, 0.3)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1) +
    geom_point(aes(y = Mean), size = 1.3, alpha = 1, stroke = 0.5) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 10),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.position = c(0.68, 1), 
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
# Plotting functions with normalization
###########################################
plot_fig3b_normed <- function(simdata_dist, tag = "B", colors = cityblock_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  subdata <- simdata_dist[
    simdata_dist$Metric %in% c("cityblock_all_normed", "cityblock_in_normed", "cityblock_out_normed") & 
      simdata_dist$M <= Mmax, ]
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  labels <- c("population-level","within-party","between-party")
  
  fig3b <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 8),
          legend.title = element_text (size = 9),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Party bias (p)",
         y = "Average opinion distance",
         tag = tag) +
    scale_linetype_manual(values = c("solid","dotted","dashed"), 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Rate of opinion\nmutation (v)") +
    scale_fill_manual(values = colors, 
                      name   = "Rate of opinion\nmutation (v)")  +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1) +
    geom_point(aes(y = Mean), size = 1.3, alpha = 1, stroke = 0.5) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 10),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.position = c(0.68, 1), 
          legend.justification = c(0,1))
  
  fig3b_grob <- ggplotGrob(fig3b)
  fig3b_grob <- gtable_filter(
    fig3b_grob,
    "axis-b-2-1|axis-b-3-[12]|axis-b-4-[123]|axis-b-5-[1234]|axis-l-4-5|axis-l-3-[45]|axis-l-2-[345]|axis-l-1-[2345]", 
    trim = FALSE, invert = TRUE)
  
  fig3b <- as_ggplot(fig3b_grob)
  
  return(fig3b)
}

plot_fig4a_normed <- function(simdata_dist, tag = "A", colors = hamming_colors){
  
  simdata_dist$v <- factor(simdata_dist$v)
  
  subdata <- simdata_dist[
    simdata_dist$Metric %in% c("hamming_all_normed", "hamming_in_normed", "hamming_out_normed") & 
      simdata_dist$M <= Mmax, ]
  subdata <- subdata %>% arrange(desc(K))
  subdata$Metric <- factor(subdata$Metric)
  
  labels <- c("population-level","within-party","between-party")
  
  fig4a <- ggplot(data = subdata, 
                  # aes(x = q, colour = v, fill = v, lty = Metric)) +
                  aes(x = p, colour = v, fill = v, lty = Metric)) +
    theme_classic() +
    theme(legend.text = element_text (size = 8),
          legend.title = element_text (size = 9),
          # legend.key.size = unit(0.03, "npc"),
          panel.spacing = unit(0.2,  "lines"),
          # legend.margin = margin(t = 0, unit="npc"),
          # axis.text.y = element_text(size=12)
          title = element_text(size = 18)
    ) +
    labs(x = "Party bias (p)",
         y = "Average issue distance",
         tag = tag) + 
    scale_linetype_manual(values = c("solid","dotted","dashed"), 
                          labels = labels,
                          name   = "Distance type") + 
    scale_color_manual(values = colors,
                       name   = "Rate of opinion\nmutation (v)" ) +
    scale_fill_manual(values = colors,
                      name   = "Rate of opinion\nmutation (v)" )  +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(limits = c(qmin, qmax), 
                       breaks = seq(qmin, qmax, qinc)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
    geom_line(aes(y = Mean), size = 0.4, alpha = 1) +
    geom_point(aes(y = Mean), size = 1.3, alpha = 1, stroke = 0.5) + 
    facet_rep_grid( M2 ~ K2, repeat.tick.labels = TRUE,
                    # switch = "both", 
                    drop = TRUE) +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 10),
          axis.line = element_line() ) +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 8), 
          # strip.background = element_blank(),
          legend.position = c(0.68, 1), 
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
# Figure 3 without normalization, without opinion dispersion
###########################################
# change plot size
if(Mmax > 3){figW <- figW*1.5}

# save multiplot
fig3b <- plot_fig3b(simdata_dist, "A")
fig3c <- plot_fig4a(simdata_dist, "B")

if(saveplots == 1){
  
  plottype <- paste0("fig3_unnormed")
  png(filename = paste0("plots/figs/", plottype, "_", vlabel, "_",
                        format(Sys.Date(), format="%y%m%d"), ".png"),
      width = figW*1.9/1.5, height = figW*ratio*1.75*2, units = "in", res = 300)
  multiplot(fig3b, fig3c, cols = 1)
  dev.off()
  
}

###########################################
# Figure 3 with normalization, without opinion dispersion
###########################################
# save multiplot
fig3b_normed <- plot_fig3b_normed(simdata_dist, "A")
fig3c_normed <- plot_fig4a_normed(simdata_dist, "B")

if(saveplots == 1){
  
  plottype <- paste0("fig3_normed_horizontal")
  png(filename = paste0("plots/figs/", plottype, "_", vlabel, "_",
                        format(Sys.Date(), format="%y%m%d"), ".png"),
      width = figW*3.8/1.5, height = figW*ratio*1.75, units = "in", res = 300)
  multiplot(fig3b_normed, fig3c_normed, cols = 2)
  dev.off()
  
  plottype <- paste0("fig3_normed")
  png(filename = paste0("plots/figs/", plottype, "_", vlabel, "_",
                        format(Sys.Date(), format="%y%m%d"), ".png"),
      width = figW*1.9/1.5, height = figW*ratio*1.75*2, units = "in", res = 300)
  multiplot(fig3b_normed, fig3c_normed, cols = 1)
  dev.off()
  
}

