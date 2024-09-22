## Model relationships between source water proportions estimated by mixing models and catchments characteristics and precipitation

# Load packages
library(tidyverse)
library(R.utils)
library(ggpubr)
library(gridExtra)
library(googledrive)
library(grid)
library(R2jags)


### Load summary and catchment characteristics data
# Load summary data from mixing model output
summary.dat <- read.csv("Output_data/All_summary_data.csv")

# Load catchment characteristics
Site_char <- read.csv("Raw_data/Site_characteristics.csv") %>% pivot_longer(cols = 2:6, names_to = "Site") %>% 
  pivot_wider(names_from = Attribute, values_from = value)
Site_char$Catchment_area <- Site_char$Catchment_area/100 # m^2 to km^2

# Join summary data and catchment characteristics
dat.char <- full_join(summary.dat, Site_char)


                                                                  

############## Linear models of source proportions and site characteristics ##############

## JAGS model with random intercept by site
lm_jags.int <- function(){
  
  # Likelihood
  for (i in 1:n) {
    mu[i] <- alpha[year[i]]+beta*car[i]
    N[i] ~ dnorm(mu[i], tau_res)
    y[i]~dnorm(N[i], pow(Ysd[i], -2)) # Encorporate standard deviation of water source proportions
  }
  
  # Priors
  mu_int~dnorm(0.5, 1) # Mean hyperparameter for random intercepts
  sigma_int~dunif(0, 100) # SD hyperparameter for random intercepts
  tau_int <- 1/(sigma_int*sigma_int)
  for (i in 1:Nyears) {
    alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
  }
  beta ~ dnorm(0, 0.2) # slope
  sigma_res~dunif(0, 100) # Residual standard deviation
  tau_res <- 1/(sigma_res*sigma_res) # Residual precision
  
} # End model

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "sigma_res")

# Chains and iterations
n.chains = 3
n.iter = 20000 
n.burnin = 5000



## Models of well proportions
# Select well proportions data
well <- dat.char %>% filter(Source == "well", Year %in% c("2019", "2020", "2021", "2022")) # Select only years that we have data from all 5 sites
Nyears <- length(levels(as.factor(well$Year)))
well$Year <- rep(1:4,4) # Change year into a level for JAGS model

# % Deciduous model
jagsdata <- with(well, list(y = Mean, car = Deciduous, year = Year, Ysd = SD,
                               n = length(Mean), Nyears = Nyears))
fit_decid <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_decid

print(fit_decid, intervals=c(0.005, 0.5, 0.995))


# Slope model
jagsdata <- with(well, list(y = Mean, car = Slope, year = Year, Ysd = SD,
                            n = length(Mean), Nyears = Nyears))
fit_slope <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope
print(fit_slope, intervals=c(0.005, 0.5, 0.995))

# % Permafrost model
jagsdata <- with(well, list(y = Mean, car = PF_percent, year = Year, Ysd = SD,
                            n = length(Mean), Nyears = Nyears))
fit_PF <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_PF
print(fit_PF, intervals=c(0.005, 0.5, 0.995))

# Catchment size model
jagsdata <- with(well, list(y = Mean, car = Catchment_area, year = Year, Ysd = SD,
                            n = length(Mean), Nyears = Nyears))
fit_size <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_size



## Models of precipitation proportions
# Select precipitation proportions data
ppt <- dat.char %>% filter(Source == "ppt/snow", Year %in% c("2019", "2020", "2021", "2022"), Site != "VAUL")
Nyears <- length(levels(as.factor(ppt$Year)))
ppt$Year <- rep(1:4,4)

# % Decid model
jagsdata <- with(ppt, list(y = Mean, car = Deciduous, year = Year, Ysd = SD,
                           n = length(Mean), Nyears = Nyears))
fit_decid.p <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_decid.p
print(fit_decid.p, intervals=c(0.005, 0.5, 0.995))


# Slope model
jagsdata <- with(ppt, list(y = Mean, car = Slope, year = Year, Ysd = SD,
                           n = length(Mean), Nyears = Nyears))
fit_slope.p <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.p
print(fit_slope.p, intervals=c(0.005, 0.5, 0.995))

# % Permafrost model
jagsdata <- with(ppt, list(y = Mean, car = PF_percent, year = Year, Ysd = SD,
                           n = length(Mean), Nyears = Nyears))
fit_PF.p <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
               n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_PF.p

# Catchment size model
jagsdata <- with(ppt, list(y = Mean, car = Catchment_area, year = Year, Ysd = SD,
                           n = length(Mean), Nyears = Nyears))
fit_size.p <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.int,
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_size.p

print(fit_size.p, intervals=c(0.005, 0.5, 0.995))


## Make table of model output from all above models
# Combine output from all models
catchment_mods <- rbind(fit_decid$BUGSoutput$summary, fit_decid.p$BUGSoutput$summary,
      fit_slope$BUGSoutput$summary, fit_slope.p$BUGSoutput$summary, fit_PF$BUGSoutput$summary, fit_PF.p$BUGSoutput$summary, 
      fit_size$BUGSoutput$summary, fit_size.p$BUGSoutput$summary)

catchment_mods <- as.data.frame(catchment_mods)
catchment_mods$source <- rep(c(rep("well", 9), rep("ppt/snow", 9)), 4)
catchment_mods$characteristic <- c(rep("Deciduous", 18), rep("Slope", 18), rep("Permafrost", 18), rep("Size", 18))
catchment_mods$parameter = rep(c("alpha_1", "alpha_2", "alpha_3", "alpha_4", "beta", "mu_int", "sigma_int", "sigma_res", "deviance"), 8)
head(catchment_mods)

# Export model results
write.csv(catchment_mods, "Tables/Catchment characteristics models.csv", row.names = F)





### Plot model results for relationships between source water proportions and catchment characteristics
# Extract draws for intercepts and slopes for CI

## % deciduous model
# Groundwater proportion model
decid.out.w <- as.data.frame(as.matrix(as.mcmc(fit_decid)))
decid.out.w$source <- "well"

# Precipitation proportion model
decid.out.p <- as.data.frame(as.matrix(as.mcmc(fit_decid.p)))
decid.out.p$source <- "ppt/snow"

# Combine groundwater and precipitation model outputs
decid.out <- rbind(decid.out.w, decid.out.p)
names(decid.out) <- c("a1", "a2", "a3", "a4", "beta", "deviance", "mu_int", "sigma_int", "sigma_res", "source")
decid.out$var <- "Deciduous"

## Slope model
# Groundwater proportion model
slope.out.w <- as.data.frame(as.matrix(as.mcmc(fit_slope)))
slope.out.w$source <- "well"

# Precipitation proportion model
slope.out.p <- as.data.frame(as.matrix(as.mcmc(fit_slope.p)))
slope.out.p$source <- "ppt/snow"

# Combine groundwater and precipitation model outputs
slope.out <- rbind(slope.out.w, slope.out.p)
names(slope.out) <- c("a1", "a2", "a3", "a4", "beta", "deviance", "mu_int", "sigma_int", "sigma_res", "source")
slope.out$var <- "Slope"

## Permafrost model
# Groundwater proportion model
PF.out.w <- as.data.frame(as.matrix(as.mcmc(fit_PF)))
PF.out.w$source <- "well"

# Precipitation proportion model
PF.out.p <- as.data.frame(as.matrix(as.mcmc(fit_PF.p)))
PF.out.p$source <- "ppt/snow"

# Combine groundwater and precipitation model outputs
PF.out <- rbind(PF.out.w, PF.out.p)
names(PF.out) <- c("a1", "a2", "a3", "a4", "beta", "deviance", "mu_int", "sigma_int", "sigma_res", "source")
PF.out$var <- "PF_percent"

## Catchment size model
# Precipitation proportion model
s.out.p <- as.data.frame(as.matrix(as.mcmc(fit_size.p)))
s.out.p$source <- "ppt/snow"

# Groundwater proportion model
s.out.w <- as.data.frame(as.matrix(as.mcmc(fit_size)))
s.out.w$source <- "well"

# Combine groundwater and precipitation model outputs
s.out <- rbind(s.out.w, s.out.p)
names(s.out) <- c("a1", "a2", "a3", "a4", "beta", "deviance", "mu_int", "sigma_int", "sigma_res", "source")
s.out$var <- "Catchment_area"


## Combine all model outputs
out.all <- rbind(PF.out, slope.out, decid.out, s.out)

## Specify endpoints for line segments
ends <- dat.char[,c(1,2,10, 11,12,13,14,15,18)] %>% pivot_longer(cols = 6:9, names_to = "var") %>% 
  group_by(var, Source) %>% 
  summarise(x1 = min(value), x2 = max(value)) %>% na.omit()
ends


## Code for confidence intervals
# Function for creating predicted lines from slopes and intercepts
calc <- function(v, s, i){
  x <- seq(ends[ends$var == v & ends$Source == s,]$x1, ends[ends$var == v & ends$Source == s,]$x2, by = 0.1)
  ints <- out.all[out.all$v == v,c(7,10)]
  slps <- out.all[out.all$v == v,c(5,10)] 
  
  
  y = ints[ints$source == s,][i,1] + slps[slps$source == s,][i,1]*x
  df <- data.frame(x=x, y=y, rep = rep(paste(i), length(y)), source = rep(paste(s), length(y)),
                   var = rep(paste(v), length(y)))
  return(df)
}

# Make a list of all element combinations to calculate predicted lines for
args <- expand.grid( c = unique(out.all$var), b = c("well", "ppt/snow"), a = 1:(length(decid.out$beta)/2))

# Apply function to output points for predicted lines
preds <- mapply(calc, v = args$c, s = args$b, i = args$a)

# Big messy reordering of output
arr <- array2DF(preds)
x <- arr %>% filter(Var1 == "x") %>% select(Value)
y <- arr %>% filter(Var1 == "y") %>% select(Value)
var <- arr %>% filter(Var1 == "var") %>% select(Value)
source <- arr %>% filter(Var1 == "source") %>% select(Value)
rep <- arr %>% filter(Var1 == "rep") %>% select(Value)
pred <- data.frame(x=x, y=y, source = source, rep = rep, var = var)
names(pred) <- c("x", "y", "source", "rep", "var")
pred$x <- as.numeric(pred$x)
pred$y <- as.numeric(pred$y)
head(pred)

# Calculate quantiles to get lines for CIs
pred.int <- pred %>% group_by(x, source, var) %>% summarise(upper = quantile(y, probs = .975),
                                                     lower = quantile(y, probs = .025),
                                                     median = quantile(y, probs = .50))

### Plot water source proportions and catchment characteristics
# Change order of legend                                                    
dat.char$Source <- factor(dat.char$Source, levels=c("ppt/snow", "soil_water", "well", "upwelling"))

# % Deciduous plot
p1 <- dat.char %>% filter(Site != "VAUL", Source != "Upwelling") %>% 
  ggplot(aes(Deciduous, Mean, color = Source)) + geom_point(cex = 2) +
  theme_bw() + xlab("Deciduous cover (%)")  + ylab(NULL)  + scale_color_manual(labels = c("Precipitation", "Soil water", "Groundwater"), 
                                                                 values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  geom_line(data = pred.int[pred.int$var == "Deciduous",], aes(x, median, color = source), stat="smooth") +
  geom_line(data = pred.int[pred.int$var == "Deciduous",], aes(x,upper, color = source), stat="smooth", alpha = 0.5) +
  geom_line(data = pred.int[pred.int$var == "Deciduous",], aes(x,lower, color = source), stat="smooth", alpha = 0.5) +
  ylim(0,0.9)

# % Slope plot
p2 <- dat.char %>% filter(Site != "VAUL") %>% ggplot(aes(Slope, Mean, color = Source)) + geom_point(cex = 2) +
  theme_bw() + xlab("Slope (%)")  + ylab(NULL) + scale_color_manual(labels = c("Precipitation", "Soil water", "Groundwater"), 
                                                                  values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  geom_line(data = pred.int[pred.int$var == "Slope" & pred.int$source == "well",], aes(x, median, color = source), stat="smooth") +
  geom_line(data = pred.int[pred.int$var == "Slope" & pred.int$source == "well",], aes(x,upper, color = source), stat="smooth", alpha = 0.5) +
  geom_line(data = pred.int[pred.int$var == "Slope" & pred.int$source == "well",], aes(x,lower, color = source), stat="smooth", alpha = 0.5)  +
  ylim(0,0.9)


# % Permafrost
p3 <- dat.char %>% filter(Site != "VAUL") %>% ggplot(aes(PF_percent, Mean, color = Source)) + geom_point(cex = 2) +
  theme_bw() + xlab("Permafrost extent (%)")  + ylab(NULL)   + scale_color_manual(labels = c("Precipitation", "Soil water", "Groundwater"), 
                                                                                values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  geom_line(data = pred.int[pred.int$var == "PF_percent" & pred.int$source == "well",], aes(x, median, color = source), stat="smooth") +
  geom_line(data = pred.int[pred.int$var == "PF_percent" & pred.int$source == "well",], aes(x,upper, color = source), stat="smooth", alpha = 0.5) +
  geom_line(data = pred.int[pred.int$var == "PF_percent" & pred.int$source == "well",], aes(x,lower, color = source), stat="smooth", alpha = 0.5)  +
  ylim(0,0.9)


# Catchment size
p4 <- dat.char %>% filter(Site != "VAUL") %>% ggplot(aes(Catchment_area, Mean, color = Source)) + geom_point(cex = 2) +
  theme_bw() + xlab(bquote("Catchment area " (km^2)))  + ylab(NULL)   + scale_color_manual(labels = c("Precipitation", "Soil water", "Groundwater"), 
                                                                                         values = c("#66c2a5", "#fc8d62", "#8da0cb"))+
  geom_line(data = pred.int[pred.int$var == "Catchment_area" & pred.int$source == "ppt/snow",], aes(x, median, color = source), stat="smooth") +
  geom_line(data = pred.int[pred.int$var == "Catchment_area" & pred.int$source == "ppt/snow",], aes(x,upper, color = source), stat="smooth", alpha = 0.5) +
  geom_line(data = pred.int[pred.int$var == "Catchment_area" & pred.int$source == "ppt/snow",], aes(x,lower, color = source), stat="smooth", alpha = 0.5) +
  ylim(0,0.9)

## Combined plot (Figure 4)
p.quad <- ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "bottom")
p.quad <- annotate_figure(p.quad, left = text_grob("                 Proportion", rot = 90)) # Add y-axis label
p.quad









############## Linear models of source proportions and summer rain ##############
### Wrangle rain data and calculate cumulative summer rain
# Read in data from our rain gauges (at French, Stuart, and Vault creeks), LTER gauges (Poker Creek), or Eielson AFB
Precip.dat <- read.csv("DoD_precip_all.csv")
Precip.dat$datetimeAK <- as.POSIXct(Precip.dat$datetimeAK, tz = "America/Anchorage")
Precip.dat$Julian_day <- as.POSIXlt(Precip.dat$datetimeAK)$yday # Calculate Julian day

# Calculate days since zero snow
Precip.dat <- Precip.dat %>% mutate(Days_melt = case_when(year == "2015" ~ Julian_day - 119,
                                                          year == "2018" ~ Julian_day - 130,
                                                          year == "2019" ~ Julian_day - 115,
                                                          year == "2020" ~ Julian_day - 130,
                                                          year == "2021" ~ Julian_day - 119,
                                                          year == "2022" ~ Julian_day - 131))



# Calculate cumulative precipitation between earliest common start date and end of window (143)
Precip.dat %>% group_by(Site, year) %>% summarise(min = min(Days_melt),
                                                  max = max(Days_melt))

Precip.cum.year <- Precip.dat %>% group_by(Site, year) %>% filter(Days_melt >= 35 & Days_melt <= 143) %>% 
  summarise(Precip = sum(inst_rainfall_mm, na.omit = T))


## Assign what rain data to use for sites or years that are missing data
# Eielson 2021 is much too high and doesn't appear accurate

add_precip <- data.frame(Site = c("FRCH", "FRCH", "FRCH", rep("MOOS",5), "STRT", "STRT"), year = c("2018","2021","2022", "2018", "2019", "2020", "2021", "2022", "2019", "2022"), Precip = c(162, 188, 127, 162, 364, 357, 188, 127, 364, 127))
Precip.cum.year$year <- as.character(Precip.cum.year$year)
Precip.cum.year <- rbind(Precip.cum.year[-c(1:4),], add_precip)
names(Precip.cum.year) <- c("Site", "Year", "Precip")
summary.dat$Year <- as.character(summary.dat$Year)

# Join added data
summary.dat.precip <- unique(left_join(summary.dat, Precip.cum.year))
summary.dat.precip$Precip <- as.numeric(summary.dat.precip$Precip)




## Average cumulative rain across sites for each year
avg <- summary.dat.precip %>% group_by(Year) %>% summarise(ppt = mean(Precip),
                                                    max = max(Precip),
                                                    min = min(Precip),
                                                    sd = sd(Precip))
# Plot
avg %>% ggplot(aes(ppt, Year)) + geom_bar(stat = "identity", fill = "grey") +
  geom_errorbar(aes(xmin=ppt - sd, xmax=ppt + sd), width=.5) + theme_bw() + xlab("Cumulative rain (mm)") +
  ylab("") + theme(text = element_text(size = 10), axis.title.x = element_text(size = 10.5))

ggsave("Figures/Yearly precip avg.png", width = 2.75, height = 1.1)

### Model precipitation proportions and cumulative summer rain
# Select precipitation data
ppt <- summary.dat.precip %>% filter(Source == "ppt/snow")
ppt$Precip_cm <- ppt$Precip/10 # mm -> cm

# Assign levels for model
PPT <- ppt %>% mutate(across(Site, ~case_when(Site== "FRCH" ~ 1,
                                              Site== "MOOS" ~ 2,
                                              Site== "POKE" ~ 3,
                                              Site== "STRT" ~ 4,
                                              Site== "VAUL" ~ 5)))


# JAGS model
## Independent intercepts and slopes
lm_jags.4 <- function(){
# Priors
mu_int~dnorm(0.5, 1) # Mean hyperparameter for random intercepts
sigma_int~dunif(0, 100) # SD hyperparameter for random intercepts
tau_int <- 1/(sigma_int*sigma_int)
mu_slope~dnorm(0, 1) # Mean hyperparameter for random slopes
sigma_slope~dunif(0, 100) # SD hyperparameter for slopes
tau_slope <- 1/(sigma_slope*sigma_slope)
for (i in 1:Nsites) {
    alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
    beta[i]~dnorm(mu_slope, tau_slope) # Random slopes
}
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res) # Residual precision
# Likelihood
for (i in 1:n) {
    mu[i] <- alpha[site[i]]+beta[site[i]]*ppt[i]
    y[i]~dnorm(N[i], pow(Ysd[i], -2))
    N[i] ~ dnorm(mu[i], tau_res)
}
}


# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "mu_slope", "sigma_slope", "sigma_res")

# Chains and iterations
n.chains = 3
n.iter = 20000 
n.burnin = 5000

Nsites <- 5

# Run model
jagsdata.4 <- with(PPT, list(y = Mean, ppt = Precip_cm, site = Site, Ysd = SD,
                           n = length(Mean), Nsites = Nsites))
fit_cumul.2 <- jags(data = jagsdata.4, parameters.to.save = params, model.file = lm_jags.4,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_cumul.2


### Plot model output for precipitation proportion and cumulative summer rain model
# Organize model output
mod.sum <- as.data.frame(fit_cumul.2$BUGSoutput$summary)
mod.sum <- mod.sum[1:10, ]
mod.sum$Site <- rep(c("FRCH", "MOOS", "POKE", "STRT", "VAUL"), 2)
mod.sum$param <- c(rep("alpha", 5), rep("beta", 5))
mod.sum <- mod.sum %>% select(mean, param, Site) %>% pivot_wider(values_from = mean, names_from = param)


# Extract draws for intercepts and slopes for CI
turnout.out <- as.data.frame(as.matrix(as.mcmc(fit_cumul.2)))
names(turnout.out) <- c("a1", "a2", "a3", "a4", "a5", "b1", "b2", "b3", "b4", "b5", "deviance", "mu_int", "mu_slope", "sigma_int", "sigma_res", "sigma_slope")

# Endpoints for segments
ends <- summary.dat.precip %>% filter(Source == "ppt/snow") %>% select(Precip, Site) %>% group_by(Site) %>% 
  summarise(p1 = min(Precip)/10, p2 = max(Precip)/10)

# Combine coefficients and endpoints
mod.sum.ends <- full_join(ends, mod.sum)
mod.sum.ends$y1 <- mod.sum.ends$alpha+mod.sum.ends$beta*mod.sum.ends$p1
mod.sum.ends$y2 <- mod.sum.ends$alpha+mod.sum.ends$beta*mod.sum.ends$p2
mod.sum.ends


## Code for confidence intervals
# Create points from slopes and intercepts
nrow <- nrow(turnout.out$a)
x <- seq(0,60,1)
ints <- turnout.out[,1:5]
slps <- turnout.out[,6:10]
Site <- 1:5

# Function for creating predicted lines from slopes and intercepts
calc <- function(s, i){
y = ints[s][i,] + slps[s][i,]*x
df <- data.frame(x=x, y=y, rep = rep(paste(i), length(y)), Site = rep(paste(s), length(y)))
return(df)
}

# Calculate points for prediction lines for all element combinations
args <- expand.grid( b = 1:5, a = 1:length(turnout.out$a1))
preds <- mapply(calc, s = args$b, i = args$a)
arr <- array2DF(preds)

# Wrangle output
x <- arr %>% filter(Var1 == "x") %>% select(Value)
y <- arr %>% filter(Var1 == "y") %>% select(Value)
Site <- arr %>% filter(Var1 == "Site") %>% select(Value)
rep <- arr %>% filter(Var1 == "rep") %>% select(Value)
pred <- data.frame(x=x, y=y, Site = Site, rep = rep)
names(pred) <- c("x", "y", "Site", "rep")
pred$x <- as.numeric(pred$x)
pred$y <- as.numeric(pred$y)
head(pred)

# Calculate CI at every x-value
pred.int <- pred %>% group_by(x, Site) %>% summarise(upper = quantile(y, probs = .975),
                                   lower = quantile(y, probs = .025))
pred.int$Source <- "ppt/snow"

# Convert levels back into site names
pred.int <- pred.int %>% mutate(across(Site, ~case_when(Site== "1" ~ "FRCH",
                                              Site== "2" ~ "MOOS",
                                              Site== "3" ~ "POKE",
                                              Site== "4" ~ "STRT",
                                              Site== "5" ~ "VAUL")))
     
# Filter CI lines by endpoints
pred.int <- pred.int %>% filter(x >= 12.7 & x <= 57.1,
                                x >= 18 | Site %in% c("FRCH", "MOOS", "STRT"),
                                x >= 19 | Site %in% c("FRCH", "MOOS", "STRT", "POKE"),
                                x <= 42.2 | Site %in% c("VAUL"),
                                x <= 36.7 | Site %in% c("STRT", "VAUL"),
                                x <= 36.4 | Site %in% c("STRT", "VAUL", "POKE"))


## Plot
# Set site names
site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")

# Change order of legend                                                    
summary.dat.precip$Source <- factor(summary.dat.precip$Source, levels=c("ppt/snow", "soil_water", "well", "upwelling"))


p_precip <- summary.dat.precip %>% 
  ggplot(aes(Precip/10, Mean, color = Source)) +
  geom_point(cex = 2) +
  facet_wrap(~ Site, labeller = labeller(Site = site_labs)) + 
  geom_errorbar(aes(ymin=Mean - SD, ymax=Mean + SD), width=4) +
  scale_color_manual(labels = c("Precipitation", "Soil water", "Groundwater", "Vault groundwater"), 
                     values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
  labs(color = "Source", y = "Proportion", x = "Cumulative rain (cm)") + theme_bw() +
  theme(strip.background =element_rect(fill="white"), legend.position = c(0.85,0.22)) +
  geom_segment(data = mod.sum.ends[mod.sum.ends$Site != "POKE",], aes(x = p1, xend = p2, y = y1, yend = y2), col = "#66c2a5") +
  geom_line(data = pred.int[pred.int$Site != "POKE",], aes(x,upper), color = "#66c2a5", stat = "smooth", alpha = 0.5) +
  geom_line(data = pred.int[pred.int$Site != "POKE",], aes(x,lower), color = "#66c2a5", stat = "smooth", alpha= 0.5) # Confidence interval

p_precip

ggsave("Figures/Cumulative precip and proportions window.png", width = 4.5, height = 3.5)
