### Plot raw data that will be used in the mixing models
# This code creates Figure 2

library(MixSIAR)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpubr)

## Load data select for mixing models
all_dat <- read.csv("Raw_data/All_mixing_data.csv")
mix_d <- all_dat %>% select(Days_melt, Site, Chloride_uM, Magnesium_uM) %>% na.omit()

## Create data frame with mixture data
mix_dat <- all_dat  %>% group_by(Days_melt, Type, year,Site,Date) %>% 
  summarise(Chloride_uM = mean(Chloride_uM, na.rm = TRUE),
            Magnesium_uM = mean(Magnesium_uM, na.rm = TRUE)) %>% na.omit()


## Create dataframe with just source data
sources_dat <- all_dat %>% filter(Type %in% c("well", "ppt/snow", "soil_water")) %>% select(Type, Chloride_uM, Magnesium_uM)

sources_VAUL <- all_dat %>% filter(Type %in% c("upwelling", "ppt/snow", "soil_water")) %>% select(Type, Chloride_uM, Magnesium_uM)

# Create one dataframe with sources matched to sites
length(sources_dat[,c(1,2,3)]$Chloride_uM)
length(sources_VAUL$Chloride_uM)
s_dat <- data.frame(Site = c(rep("MOOS", 20), rep("FRCH", 20), rep("STRT", 20), rep("POKE", 20), rep("VAUL", 21)), 
                    Chloride_uM = c(rep(sources_dat$Chloride_uM,4), sources_VAUL$Chloride_uM), Magnesium_uM = c(rep(sources_dat$Magnesium_uM,4), sources_VAUL$Magnesium_uM), Type = c(rep(sources_dat$Type,4), sources_VAUL$Type))

### Plots
# Function to calculate standard error function
SE <- function(x) sd(x, na.rm=TRUE)/sqrt(length(na.omit(x)))

# Calculate means and standard errors of source groups
s_m <- s_dat %>% filter(!is.na(Magnesium_uM)) %>% group_by(Type, Site) %>% dplyr::summarise(
                 Chloride = mean(Chloride_uM), 
                 SEC = SE(Chloride_uM), 
                 Magnesium = mean(Magnesium_uM), 
                 SEM = SE(Magnesium_uM),
                 n = length(Chloride_uM)) %>% filter(Type != "well" | Site != "VAUL")
s_m
names(s_m) <- c("Type", "Site", "Chloride_uM", "SEC", "Magnesium_uM", "SEM", "n")

s_m$Type <- as.factor(s_m$Type)
levels(s_m$Type) <- c("Precipitation", "Soil water", "Groundwater", "Groundwater")

## Plot with two panels: One with source means, standard errors and stream samples from Vault Creek, and one colored by site for the other four sites
site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")

s_vaul <- s_m %>% filter(Site == "VAUL") # Select Vault sources

cols <- c('slategray3', 'palegreen2', '#086687', 'khaki', 'lightskyblue')
cols <- c('darkgoldenrod1', '#085908', '#47A8D9', 'ivory4', '#730AB8')

# Plot without Vault
p11 <- mix_d %>% filter(Site != "VAUL") %>% ggplot(aes(Chloride_uM, Magnesium_uM)) + geom_point(aes(color = Site), alpha = 0.3) + 
  geom_point(data = s_m, aes(Chloride_uM, Magnesium_uM)) +
  geom_errorbar(data = s_m,                     # add x-axis error bars
                mapping = aes(x = Chloride_uM,
                              ymin = Magnesium_uM- SEM, 
                              ymax = Magnesium_uM+SEM), width = 3) +
  geom_errorbar(data = s_m,                     # add y-axis error bars
                mapping = aes(y = Magnesium_uM,
                              xmin = Chloride_uM - SEC,
                              xmax = Chloride_uM + SEC), width = 20)  +
  ylim(-8,615) +
  theme_bw() +
  labs(color = "", y = expression(paste("Magnesium (", mu, "M)", sep = "")), x = NULL) +
  theme(text = element_text(size=16)) +
  theme(strip.background =element_rect(fill="white"), axis.title.y = element_text(size = 14)) +
  scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart", "Vault")) +
  guides(color = guide_legend(nrow = 1, override.aes=list(size = 3, alpha = 0.8))) +
  geom_text(data = s_m, aes(label=Type), hjust= -0.2, vjust= 0.3, size = 3.5, fontface = "bold") +
  geom_point(data = s_vaul[1,], aes(Chloride_uM, Magnesium_uM, color = Site), alpha = 0.5)
p11  

# Plot of Vault
p12 <- mix_d %>% filter(Site == "VAUL") %>% ggplot(aes(Chloride_uM, Magnesium_uM)) + geom_point(aes(color = Site), col = "#730AB8", alpha = 0.5) + 
  geom_point(data = s_vaul, aes(Chloride_uM, Magnesium_uM)) +
  geom_errorbar(data = s_vaul,                     # add x-axis error bars
                mapping = aes(x = Chloride_uM,
                              ymin = Magnesium_uM- SEM, 
                              ymax = Magnesium_uM+SEM), width = 3) +
  geom_errorbar(data = s_vaul,                     # add y-axis error bars
                mapping = aes(y = Magnesium_uM,
                              xmin = Chloride_uM - SEC,
                              xmax = Chloride_uM + SEC), width = 200)  +
  theme_bw() +
  labs(color = "", y = "", x = NULL) +
  theme(text = element_text(size=16)) +
  theme(strip.background =element_rect(fill="white")) +
  guides(color = guide_legend(nrow = 3, override.aes=list(size = 5, alpha = 1))) +
  geom_text(aes(x=21, y=0, label = "Precipitation"), size = 3.5, fontface = "bold") +
  geom_text(aes(x = 13, y = 2900, label = 'Groundwater'), size = 3.5, fontface = "bold")+
  geom_text(aes(x = 60, y = 300, label = 'Soil water'), size = 3.5, fontface = "bold")

p12

# Combine frames
p <- ggarrange(p11,p12, common.legend = T)
p
annotate_figure(p, bottom = text_grob(expression(paste("Chloride (", mu, "M)", sep = "")), size = 14)) + theme(plot.background = element_rect(fill = 'white', color = 'white'))

ggsave('Figures/Biplot small.png', width = 7.2, height = 3.8)
