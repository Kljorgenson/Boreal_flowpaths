## Compile and plot summary statistics from mixing models

library(tidyverse)
library(RColorBrewer)
library(R.utils)
library(here)
library(ggpubr)
library(MixSIAR)
library(gridExtra)
library(googledrive)
library(grid)

## Load model data from individual model files
# Set output options for MixSIAR models
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = TRUE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "posterior_density",
                       sup_pairs = TRUE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = F,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       return_obj = TRUE) 

# Create list of all RData files from MixSIAR models
filelist = list.files(path = "MixSIAR_model_output/", pattern = "*_snw.RData", full.names = TRUE)
filelist
mods <- gsub("_com_snw.RData", "", as.character(filelist), perl = TRUE)
models <- gsub("MixSIAR_model_output/", "", as.character(mods), perl = TRUE)
models 


## Save posteriors, statistics, and diagnostics for each model in lists
# Posteriors
posts <- lapply(models, function(x){
  load(paste0("MixSIAR_model_output/", x, "_com_snw.RData"))
  
  posts <- output_posteriors(jags.1, mix, source, output_options)
 
   })

# Summary statistics
stats <- lapply(models, function(x){
  load(paste0("MixSIAR_model_output/", x, "_com_snw.RData"))
  
  stats <- data.frame(output_stats(jags.1, mix, source, output_options))
})

# Diagnostics
diags <- lapply(models, function(x){
  load(paste0("MixSIAR_model_output/", x, "_com_snw.RData"))

  diags <- do.call(cbind, output_diagnostics(jags.1, mix, source, output_options))
  })


# Name list items with site and year
names(posts) <- models
names(stats) <- models
names(diags) <- models


### Check diagnostics for all models
# Combine diagnostic data into dataframe
diag.dat <- do.call(rbind, diags)
diag.dat$name <- rownames(diag.dat)
diag.dat <- diag.dat %>% separate(name, into = c("SY", "v1", "v2", "v3"), sep = '[.]') %>% 
  separate(SY, into = c("Site", "Year"), sep = '[_]')
rownames(diag.dat) <- NULL
names(diag.dat) <- c("gelman.point", "gelman.upperCI", "geweke.chain1", "geweke.chain2", "geweke.chain3", "Site", "Year", "v1", "v2", "v3")

# Summarize diagnostic data results
diag.results <- diag.dat %>% group_by(Site, Year) %>% summarise(geweke.chain1 = sum(abs(geweke.chain1)>1.96),
                                                                geweke.chain2 = sum(abs(geweke.chain2)>1.96),
                                                                geweke.chain3 = sum(abs(geweke.chain3)>1.96),
                                                                n = length(Site),
                                                                gelman1.05 = sum(gelman.point >= 1.05),
                                                                gelman1.10 = sum(gelman.point >= 1.10))

# Export diagnostic data summary
write.csv(diag.results, "Tables/Diagnostics.csv", row.names = F)



### Create plot of all posterior distributions (Figure S3)
# Select and modify individual year+site plots
F.2015 <- posts$FRCH_2015[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2015") + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=30), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(name = NULL, values = c("#66c2a5", "#fc8d62", "#8da0cb","#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater")) + theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
F.2018 <- posts$FRCH_2018[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2018") + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=30), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(name = NULL, values = c("#66c2a5", "#fc8d62", "#8da0cb","#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater")) + theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
F.2019 <- posts$FRCH_2019[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2019") + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=30), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(name = NULL, values = c("#66c2a5", "#fc8d62", "#8da0cb","#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater")) + theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
F.2020 <- posts$FRCH_2020[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2020") + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
F.2021 <- posts$FRCH_2021[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2021") + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
F.2022 <- posts$FRCH_2022[[1]] + xlab(NULL) + ylab(NULL) + ggtitle("2022") + theme(axis.text.y = element_text(size = 20, face = 'bold'), axis.text.x=element_blank(), axis.title=element_text(size=30), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
P.2019 <- posts$POKE_2019[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
P.2020 <- posts$POKE_2020[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), "#e78ac3", labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
P.2021 <- posts$POKE_2021[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
P.2022 <- posts$POKE_2022[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.y = element_text(size = 20, face = 'bold'), axis.text.x=element_blank(), axis.title=element_text(size=30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
V.2019 <- posts$VAUL_2019[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x = element_text(size = 20, face = 'bold'), axis.title=element_text(size=30), axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Groundwater", "Vault Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
V.2020 <- posts$VAUL_2020[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x = element_text(size = 20, face = 'bold'), axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
V.2021 <- posts$VAUL_2021[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x = element_text(size = 20, face = 'bold'), axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
V.2022 <- posts$VAUL_2022[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.y = element_text(size = 20, face = 'bold'), axis.text.x = element_text(size = 20, face = 'bold'), axis.title=element_text(size=30), plot.title = element_text(hjust = 0.5, size = 30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater", "Vault Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2015 <- posts$MOOS_2015[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2018 <- posts$MOOS_2018[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2019 <- posts$MOOS_2019[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2020 <- posts$MOOS_2020[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2021 <- posts$MOOS_2021[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
M.2022 <- posts$MOOS_2022[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.y = element_text(size = 20, face = 'bold'), axis.text.x=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
S.2019 <- posts$STRT_2019[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=30), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
S.2020 <- posts$STRT_2020[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), "#e78ac3", labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
S.2021 <- posts$STRT_2021[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))
S.2022 <- posts$STRT_2022[[1]] + xlab(NULL) + ylab(NULL) + ggtitle(NULL) + theme(axis.text.y = element_text(size = 20, face = 'bold'), axis.text.x=element_blank(), legend.position="none", axis.title=element_text(size=30)) + scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"), labels = c("Precipitation", "Soil water", "Groundwater"))+ theme(plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))


# Blank grids for plot organization
blank <- grid.rect(gp=gpar(col="white"))

# List of all plots
plot.list <- list(F.2022, F.2021, F.2020, F.2019, F.2018, F.2015, M.2022, M.2021, M.2020, M.2019, M.2018, M.2015, P.2022, P.2021, P.2020, P.2019, blank, blank, S.2022, S.2021, S.2020, S.2019, blank, blank, V.2022, V.2021, V.2020, V.2019, blank, blank)

# Arrange plots in grid
p.post <- ggarrange(plotlist = plot.list,
                    common.legend = F, legend = "none", ncol = 6, nrow = 5,
                    widths = c(70,50,50,50,50,50), heights = c(130,100,100,100,125))

# Annotate plots with axis labels and sites
p.post.1 <- annotate_figure(p.post, bottom = text_grob("       Proportion", size = 30))
p.post.3 <- annotate_figure(p.post.1, right = text_grob("French", size = 30, rot = 270, hjust = 4.7, vjust = -0.1))
p.post.3 <- annotate_figure(p.post.3, right = text_grob("Moose", size = 30, rot = 270, hjust = 2.8, vjust = 1.3))
p.post.3 <- annotate_figure(p.post.3, right = text_grob("Poker", size = 30, rot = 270, hjust = .5, vjust = 20.5))
p.post.3 <- annotate_figure(p.post.3, right = text_grob("Stuart", size = 30, rot = 270, hjust = -1.8, vjust = 21.9))
p.post.3 <- annotate_figure(p.post.3, right = text_grob("Vault", size = 30, rot = 270, hjust = -4.8, vjust = 23.1))
p.post.4 <- annotate_figure(p.post.3, left = text_grob("Scaled Posterior Density", size = 30, rot = 90, vjust = -0.05))

# Export plot
ggexport(plot = p.post.4, filename = "Figures/Post plot window.png", height = 1300, width = 1600)




### Summary statistics
## Compile summary statistics
summary.dat.1 <- do.call(rbind, stats)
summary.dat.1$Site <- rownames(summary.dat.1)
row.names(summary.dat.1) <- NULL
summary.dat <- summary.dat.1 %>% separate(Site, into = c("SY", NA, NA, "source"), sep = '[.]') %>% 
  separate(SY, into = c("Site", "Year"), sep = '[_]') %>% na.omit()
names(summary.dat) <- c("Mean", "SD", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "Site", "Year", "Source")

## Plot of source proportions by site and year with error bars for standard deviations (Figure 3)
site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")
names(summary.dat) <- c("Mean", "SD", "p2.5", "p5", "p25", "p50", "p75", "p95", "p97.5", "Site", "Year", "Source")
                                                   
p_wrap <- summary.dat %>% ggplot(aes(Mean, as.factor(Year), color = Source)) + geom_point(cex = 2) +
  facet_grid(rows = vars(Site), scales = 'free_y', space='free_y', 
             labeller = labeller(Site = site_labs)) + 
  geom_errorbar(aes(xmin=p2.5, xmax=p97.5), width=.5) +
  scale_color_manual(breaks=c("ppt/snow","soil_water",   "well",   "upwelling" ),
                     labels = c("Precipitation", "Soil water","Groundwater", "Vault Groundwater"), 
                     values = c("#66c2a5", "#fc8d62","#8da0cb", "#e78ac3")) +
  labs(color = NULL, y = "Year", x = "Proportion") + theme_bw() +
  theme(legend.position = "top", text = element_text(size = 15)) + # c(0,0) bottom left, c(1,1) top-right.
  theme(strip.background =element_rect(fill="white"),
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "inches")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p_wrap


ggsave(here("Figures/Mean and CI window long.png"), width = 4.5, height = 5.5)


# Calculate difference among years
summary.dat %>% group_by(Site, Source) %>% summarise(diff = max(Mean) - min(Mean))

# Export summary statistics data
write.csv(summary.dat, "Output_data/All_summary_data.csv", row.names = F)
