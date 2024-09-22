### Correlations among catchment attributes ###
# This code creates Figure S1

library(here)
library(tidyverse)
library(psych)

### Data ###
dat <- read.csv(here("Raw_data/Site_characteristics.csv"), header = F) %>% t() %>% as.data.frame()
names(dat) <- dat[1,]
dat <- dat %>% slice(-1) %>% rename(Area = Catchment_area) %>% mutate(across(Area:Yedoma, as.numeric))
dat$Site <- c("French", "Moose", "Poker", "Stuart", "Vault")


### Plot ###
pdf("Figures/FigS1.pdf")
jpeg("Figures/FigS1.jpeg", width = 1000, height = 1000, res = 160)

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "spearman"), digits = 2)
  txt <- paste0("R = ", r)
  #cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.2)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, 
         col = ccols[as.factor(dat$Site)],
         pch = 19, cex = 1.4)
}

# Create the plots
pairs(dat[,2:9], 
      lower.panel = panel.cor,
      upper.panel = upper.panel, cex.labels=0.8)             

dev.off()
