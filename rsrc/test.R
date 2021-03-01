# test.R
# scratchwork

# select_cols2   <- c("id","u","v",
#                     # "p","beta","epsilon",
#                    # "CC_final","CD_final","DC_final","DD_final",
#                    "y","z","g","h"#,"sia_sid","sia_sjd"
#                    )
# 
# simdata_bymeasure_comp <- threshdata_all[ threshdata_all$p == 0, select_cols2 ] %>% 
#   gather("variable","value"-id,-u,-v) %>%
# simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)
# 
# simdata_comp <- 
#   simdata_bymeasure_all[which(simdata_bymeasure_comp$variable %in% measure_vals),] %>%
#   rename( Fraction = value, Metric = variable ) %>%
#   group_by( M, K, u, v, p, beta, epsilon, Metric ) %>%
#   summarise(
#     Mean = mean(Fraction),
#     SD = sd(Fraction),
#     SE = sd(Fraction) / sqrt(length(Fraction)),
#     numCases = length(Fraction)
#   )

############################################
# plotting strategy distribution across p
############################################
# melt data
select_cols <- c("id","M","K","u","v","p","beta","epsilon","CC","CD","DC","DD",
                 "CC_final","CD_final","DC_final","DD_final","topinion_mean","topinion_var",
                 "cooperation_all", "cooperation_in", "cooperation_out"#, # cooperation
                 # "opn_simpson_mean","sets_simpson_mean", # simpsons indices
                 # "cityblock_all", "cityblock_in", "cityblock_out", # city block distances
                 # "hamming_all",   "hamming_in",   "hamming_out"    # hamming distances
)
id_vars       <- c("M","K","u","v","p","beta","epsilon")

measure_strat  <- c("CC","CD","DC","DD")
measure_strat2 <- c("CC_final","CD_final","DC_final","DD_final")
measure_coop   <- c("cooperation_all", "cooperation_in", "cooperation_out")

simdata_bymeasure_all <- simdata[ select_cols ] %>% 
  gather("variable","value",-M,-K,-id,-p,-u,-v,-beta,-epsilon)
simdata_bymeasure_all$value = as.numeric(simdata_bymeasure_all$value)

# assume simdata is loaded as vsweep
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

plot_fig2x <- function(simdata_strat, v, M, tag = "B", labeled = TRUE, wlegend = TRUE){
  
  subdata <- simdata_strat[simdata_strat$v == v & simdata_strat$M == M, ]
  
  ylabel <- if(labeled){"Relative abundance"}else{""}
  
  fig2x <- ggplot(subdata,
                  aes(x = p, y = Mean, color = Strategy)) +
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
    labs(x = "",
         y = ylabel,
         tag = tag) +
    ggtitle( paste0("v = ", v) ) +
    scale_color_manual(values = c("#0571b0","#92c5de","#f4a582","#ca0020")) +
    scale_y_continuous(limits = c(0.15, 0.35),
                       breaks = seq(0.0, 0.5, 0.05)) +
    # scale_x_continuous(limits = c(qmin, qmax), 
    #                    breaks = seq(qmin, qmax, qinc)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, size = 0.3) +
    geom_line(aes(group = Strategy), size = 0.4, alpha = 1, lty = 1) +
    geom_point(stat="identity", size = 1.1, alpha = 1, stroke = 0.5) # + #, shape = 21, fill = "white") + 
    # facet_grid( ~ M2, space="free_x",
    #             switch = "x", scales="free_x") +
    # theme(strip.placement = "outside") +
    # theme(axis.title.x = element_blank(), 
    #       strip.background = element_blank(),
    #       strip.text.x.bottom = element_text(angle = 45))
  
  return(fig2x)
}

# plot
# v <- 0.001
figTemp1 <- plot_fig2x(simdata_strat_all, 0.001, 1, "A")
figTemp2 <- plot_fig2x(simdata_strat_all, 0.005, 1, "B")
figTemp3 <- plot_fig2x(simdata_strat_all, 0.025, 1, "C")
figTemp4 <- plot_fig2x(simdata_strat_all, 0.125, 1, "D")
figTemp5 <- plot_fig2x(simdata_strat_all, 0.625, 1, "E")

# save plot
plottype <- paste0("figStrat_p_TEMP")

png(filename = paste0("plots/figs/", plottype, "_", 
                      format(Sys.Date(), format="%y%m%d"), "_SD.png"), # !!! change !!!
    width = figW*1.75, height = figW*ratio*2.7, units = "in", res = 300)
multiplot(figTemp1, figTemp2, figTemp3, figTemp4, figTemp5, cols = 2)
          #layout = matrix(c(1,3,2,4), ncol = 2, byrow = TRUE))
# layout = matrix(c(1,2,3,4,5,6), ncol = 3, byrow = FALSE))
dev.off()


############################################
# export simulated y, z, g, h values
############################################
require(tibble) # assume it's loaded

# assume simdata_plot is computed using FigZ.R from neutral data
simdata_neutral <- data.frame(matrix(ncol = 6, nrow = 9))
colnames(simdata_neutral) <- c("u","v","y","z","g","h")
simdata_neutral$u <- 0.001
simdata_neutral$v <- simdata_plot[simdata_plot$Metric %in% c("y") & simdata_plot$p == 0,]$v
simdata_neutral$y <- simdata_plot[simdata_plot$Metric %in% c("y") & simdata_plot$p == 0,]$Mean
simdata_neutral$z <- simdata_plot[simdata_plot$Metric %in% c("z") & simdata_plot$p == 0,]$Mean
simdata_neutral$g <- simdata_plot[simdata_plot$Metric %in% c("g") & simdata_plot$p == 0,]$Mean
simdata_neutral$h <- simdata_plot[simdata_plot$Metric %in% c("h") & simdata_plot$p == 0,]$Mean

file_name <- "data/calc_data/sim_data_yzgh.csv"
write.csv(simdata_neutral, file_name, row.names = FALSE)





