### Plot discharge data
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)

Q.2015 <- read.csv("Raw_data/Discharge/Q_2015.csv")
names(Q.2015) <- c(NA, "Site", "Q", "DateTime")
Q.2018 <- read.csv("Raw_data/Discharge/Q_2018.csv")
Q.2019 <- read.csv("Raw_data/Discharge/Q_2019.csv")
Q.2020 <- read.csv("Raw_data/Discharge/Q_2020.csv")
Q.2021 <- read.csv("Raw_data/Discharge/Q_2021.csv")
Q.2022 <- read.csv("Raw_data/Discharge/Q_2022.csv")
names(Q.2022) <- c("Site", "Q", "DateTime")


Q.dat.1 <- rbind(Q.2015[,-1],Q.2019, Q.2020, Q.2022)
Q.dat.2 <- rbind(Q.2018, Q.2021)
Q.dat <- full_join(Q.dat.1, Q.dat.2) %>% mutate(Site2 = case_when(Site == "POKE.Q.int" ~ "POKE",
                                                Site == "MOOS.Q.int" ~ "MOOS",
                                                Site == "FRCH.Q.int" ~ "FRCH",
                                                Site == "VAUL.Q.int" ~ "VAUL",
                                                Site == "STRT.Q.int" ~ "STRT",
                                                Site == "POKE" ~ "POKE",
                                                Site == "MOOS" ~ "MOOS",
                                                Site == "FRCH" ~ "FRCH",
                                                Site == "VAUL" ~ "VAUL",
                                                Site == "STRT" ~ "STRT"
)) %>% select(DateTime, Q, Site2)
names(Q.dat) <- c("DateTime", "Q", "Site")
Q.dat$DateTime <- as.POSIXct(Q.dat$DateTime, tzone = "America/Anchorage")
Q.dat$Year <- Q.dat$DateTime %>% format('%Y') %>% as.numeric()
Q.dat$Julian <- as.POSIXlt(Q.dat$DateTime)$yday %>% as.numeric()


## Model start and end dates
window <- data.frame(Year = c(2015, 2018:2022), start = c(158,192,165,174,178,172), 
                     end = c(263,274,246,271,261,270))

# Average discharge during analysis windows by site 2019-2022
Q.dat %>% filter(!Year %in% c(2015,2018),
                   Year == 2019 & Julian >= 165 & Julian <=246 |
                   Year == 2020 & Julian >= 174 & Julian <=271 |
                   Year == 2021 & Julian >= 178 & Julian <=261 |
                   Year == 2022 & Julian >= 172 & Julian <=270) %>% group_by(Site) %>% summarise(median = median(Q, na.rm = T), sd = sd(Q, na.rm = T),
                                                                                                 cv = sd/mean(Q,na.rm = T)*100)

## Q plots
cols <- c('darkgoldenrod1', '#085908', '#47A8D9', 'ivory4', '#730AB8')

# Raw Q
Q.dat %>% na.omit() %>% ggplot(aes(Julian, Q, color = Site)) + geom_point(size = 0.2) + 
  facet_wrap(~Year, ncol = 2, nrow = 3) +
  xlab("Julian day") + ylab("Discharge (L/s)") + theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=2)))  +
  theme (strip.background =element_rect (fill="white"), text = element_text(size =18), legend.text=element_text(size=18)) +
  geom_vline(data = window, aes(xintercept = start), linetype = 2) + xlim(130,300) +
  geom_vline(data = window, aes(xintercept = end), linetype = 2) + scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart", "Vault")) + 
  scale_y_log10() + labs(col = NULL)

ggsave("Figures/All Q by year.png", width = 8, height = 6)


### Q and chems comparison
## Join Q
q2019 <- read.csv("Raw_data/Discharge/Q.daily.2019.csv") 
q2020 <- read.csv("Raw_data/Discharge/Q.daily.2020.csv") 
q2021 <- read.csv("Raw_data/Discharge/Q.daily.2021.csv") %>% select(Site,Day,Q)
q2022 <- read.csv("Raw_data/Discharge/Q.daily.2022.csv") 

q <- rbind(q2019,q2020,q2021,q2022) %>% rename(Date =Day)
q$Date = as.Date(q$Date)
head(q)

# Join chems
chems.d <- read.csv("Raw_data/All_mixing_data.csv") %>% filter(Type == "stream") %>% select(Date, year,Site, Chloride_uM, Magnesium_uM) %>% na.omit() %>%
                                                               group_by(Date, year,Site) %>% 
                                                               summarise(Chloride_uM = mean(Chloride_uM, na.rm = TRUE),
                                                                         Magnesium_uM = mean(Magnesium_uM, na.rm = TRUE)) %>% na.omit()
chems.d$Date <- as.Date(chems.d$Date)

                                                             
# Join chems and Q                                                             
chems.q <- left_join(q,chems.d) %>% filter(is.na(Magnesium_uM)==F, is.na(Chloride_uM)==F) %>%
  select(Site,Date,Magnesium_uM, Chloride_uM, Q) %>% na.omit()

# Plot
p1 <- chems.q %>% filter(Site != "VAUL") %>% ggplot(aes(Q,Magnesium_uM,col = Site)) + geom_point(size = 0.5) +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5) +
  scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart", "Vault")) +
  labs(color = "", y = expression(paste("Magnesium (", mu, "M)", sep = "")), x = "Discharge (L/s)", tag = "a)") +
  scale_x_continuous(trans='log10', limits = c(100,15000)) + ylim(0,405) + 
  geom_point(data = chems.q[chems.q$Magnesium_uM > 600,], aes(Chloride_uM, Magnesium_uM, color = Site),size = 0.5) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), plot.margin = margin(4,0,4,0))

p2 <- chems.q %>% filter(Site != "VAUL") %>% ggplot(aes(Q,Chloride_uM,col = Site)) + geom_point(size = 0.5) +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5)+
  scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart")) +
  labs(color = "", y = expression(paste("Chloride (", mu, "M)", sep = "")), x = "Discharge (L/s)", tag = "b)") + ylim(0,30.02) +
  scale_x_continuous(trans='log10') +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), plot.margin = margin(4,0,4,0))

p3 <- chems.q %>% filter(Site == "VAUL") %>% ggplot(aes(Q,Magnesium_uM,col = Site)) + geom_point(size = 0.5, col = "#730AB8") +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5) +
  labs(color = "", y = expression(paste("Magnesium (", mu, "M)", sep = "")), x = NULL, , tag = "c)") +
  scale_x_continuous(trans='log10') + theme(text=element_text(size=9), axis.text.x = element_blank(), plot.margin = margin(10,10,4,4),
                                            axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))

p4 <- chems.q %>% filter(Site == "VAUL") %>% ggplot(aes(Q,Chloride_uM,col = Site)) + geom_point(size = 0.5, col = "#730AB8") +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5)+
  labs(color = "", y = expression(paste("Chloride (", mu, "M)", sep = "")), x = "Discharge (L/s)", , tag = "d)") +
  scale_x_continuous(trans='log10') + theme(text=element_text(size=9),
                                            axis.title.y = element_text(margin = margin(t = 0, r = 10.5, b = 0, l = 10)),
                                            axis.title.x = element_text(size = 11, margin = margin(t = 6, r = 0, b = 0, l = 0)), 
                                            plot.margin = margin(2,10,4,4))

p <- ggarrange(p3,p4, common.legend = T, ncol = 1, heights = c(0.83,1))
ggarrange(p1,p2, p, common.legend = T, nrow = 1, widths = c(1,1,0.7))

ggsave("Figures/Q and chems 2.png", width = 6.5, height = 2.8)


## Just 4 plots
p1 <- chems.q %>% filter(Site != "VAUL") %>% ggplot(aes(Q,Magnesium_uM,col = Site)) + geom_point(size = 0.5) +
  theme_bw() +
  scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart", "Vault")) +
  labs(color = "", y = expression(paste("Magnesium (", mu, "M)", sep = "")), x = NULL) +
  scale_x_continuous(trans='log10', limits = c(100,15000)) + ylim(0,405) +
  geom_point(data = chems.q[chems.q$Magnesium_uM > 600,], aes(Chloride_uM, Magnesium_uM, color = Site),size = 0.5)+
  theme(plot.margin = margin(4,4,4,4), text = element_text(size =10), legend.text=element_text(size=10)) +
  guides(col = guide_legend(override.aes = list(size = 2)))

p2 <- chems.q  %>% filter(Site != "VAUL") %>% ggplot(aes(Q,Chloride_uM,col = Site)) + geom_point(size = 0.5) +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5)+
  scale_color_manual(values = cols, labels = c("French", "Moose", "Poker", "Stuart")) +
  labs(color = "", y = expression(paste("Chloride (", mu, "M)", sep = "")), x = "Discharge (L/s)") +
  scale_x_continuous(trans='log10') + ylim(0,30.02)+ theme(plot.margin = margin(4,4,4,8), text = element_text(size =10), 
                                                           legend.text=element_text(size=10))

p3 <- chems.q %>% filter(Site == "VAUL") %>% ggplot(aes(Q,Magnesium_uM,col = Site)) + geom_point(size = 0.5, col = "#730AB8") +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5) +
  labs(color = "", y = NULL, x = NULL) +
  scale_x_continuous(trans='log10') + theme(plot.margin = margin(4,4,4,4), text = element_text(size =10), 
                                            legend.text=element_text(size=10))

p4 <- chems.q %>% filter(Site == "VAUL") %>% ggplot(aes(Q,Chloride_uM,col = Site)) + geom_point(size = 0.5, col = "#730AB8") +
  theme_bw() + #geom_smooth(method = 'lm', alpha = 0.3, size = 1.5)+
  labs(color = "", y = NULL, x = "Discharge (L/s)") +
  scale_x_continuous(trans='log10') + theme(plot.margin = margin(4,4,4,12), text = element_text(size =10), legend.text=element_text(size=10))

ggarrange(p1,p3,p2,p4, common.legend = T, widths = c(1,0.93), heights = c(0.93,1))+ bgcolor("white") + border('white')       

ggsave("Figures/Q and chems.png", width = 4.3, height = 4)

