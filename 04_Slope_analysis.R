# Fit slopes to lines of proportions over time from mixing model output
# Load modified functions from 00_Modified_MixSIAR_functions.R prior to running this code

# Load packages
library(R2jags)
library(tidyr)
library(ggplot2)
library(dplyr)

# Set output options for MixSIAR
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


# Make list of all .RData files with mixing model output
filelist = list.files(path = "MixSIAR_model_output/", pattern = "*_snw.RData", full.names = TRUE)
filelist
mods <- gsub("_com_snw.RData", "", as.character(filelist), perl = TRUE)
models <- gsub("MixSIAR_model_output/", "", as.character(mods), perl = TRUE)
models 



###### Calculate slopes and intercepts for all sources, sites, and years
## For each proportion within each mixing model, extract model draws for each line using a modfied version of the MixSIAR continuous variable plotting function, 
# then fit a linear line to these and extract slope and error terms.
# This takes a few hours to run!
slope.df <- lapply(models, function(x){
  # Load .RData file with model output and input
  load(paste0("MixSIAR_model_output/", x, "_com_snw.RData"))
  
  # Make empty lists for output
  slope.dat<-list()
  slope.quants <- list()

  # Loop over 3 sources  
for(s in 1:3){ 
  
# Extract regression lines from mixing model output  
c <- continuous_var(jags.1, mix, source, output_options)
dim(c)

thin <- seq(from = 1, to = 9000, by = 100) # Thin the number of regression lines for efficiency 
df <- as.data.frame(t(c[,thin,s]))
dim(df)
colnames(df) <- seq(38, 143, length.out = 200) # Set time steps in terms of days since snowmelt
df$reps <- seq(1,dim(df)[1])
df2 <- df %>% pivot_longer(cols = 1:200, values_to = 'v', names_to = 'time')
head(df2)
df2$time <- as.numeric(df2$time)
df2$val <- df2$v*100 # PROPORTIONS -> PERCENTAGES FOR EASIER INTERPRETATION OF SMALL NUMBERS

# Fit JAGS model to regression line points. A line is fit to each of the extracted lines
lm_jags_slopes <- function(){
 # Priors
    mu_int~dunif(-300, 300) # Mean hyperparameter for random intercepts
    s_int~dunif(0, 100) # SD hyperparameter for random intercepts
    tau_int <- 1/(s_int*s_int)
    mu_slope~dnorm(0, 1) # Mean hyperparameter for random slopes
    s_slope~dunif(0, 100) # SD hyperparameter for slopes
    tau_slope <- 1/(s_slope*s_slope)
    for (i in 1:90) { # length that draws were thinned to
      alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
      beta[i]~dnorm(mu_slope, tau_slope) # Random slopes
    }
    s_res~dunif(0, 100) # Residual standard deviation
    tau_res <- 1/(s_res*s_res) # Residual precision
  
  # Likelihood
  for (i in 1:n) {
    mu[i] <- alpha[rep[i]]+beta[rep[i]]*X[i]
    Y[i] ~ dnorm(mu[i], tau_res)
  }
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "s_int", "mu_slope", "s_slope", "s_res")

# Chains and iterations
n.chains = 3
n.iter = 5000
n.burnin = 1000

# Model
jagsdata <- with(df2, list(Y = val, X = time, rep = reps, n = length(df2$time)))
fit.0 <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags_slopes,
              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit.0

slope.dat[[s]] <- as.data.frame(fit.0$BUGSoutput$summary) # Save to list

slope.dat[[s]]$source <- source$source_names[s] # Name the source

}


s.dat <- do.call(rbind, slope.dat) # Combine lists
s.dat$name <- rownames(s.dat)

s.dat
}) # End



## Organize slope dataframe
names(slope.df) <- models
slope.df2 <- bind_rows(slope.df, .id = "column_label")
slope.df2$par <- substr(slope.df2$name, start = 1, stop = 4)
slope.df2 <- slope.df2 %>% separate(column_label, into = c("Site", "Year"), sep = '[_]') # Extract site and year from model label
rownames(slope.df2) <- NULL
head(slope.df2)
names(slope.df2) <- c("Site", "Year","Mean", "SD", "p2.5", "p25", "p50", "p75", "p97.5", "Rhat", "n.eff", "Source", "name", "Par")

# Save all slopes and intercepts
write.csv(slope.df2, "Output_data/Slopes.csv", row.names = F)



# Select just mean slopes
slope.dat <- slope.df2 %>% filter(Par %in% c("mu_s", "mu_i")) %>% select("Site", "Year","Mean", "SD", "p2.5", "p50", "p97.5", "Source", "Par") %>% pivot_wider(values_from = c("Mean", "SD", "p2.5", "p50", "p97.5"), names_from = "Par")
names(slope.dat) <- c("Site", "Year", "Source", "Mean_alpha", "Mean_beta", "SD_alpha", "SD_beta", "p2.5_alpha", "p2.5_beta", "p50_alpha", "p50_beta", "p97.5_alpha", "p97.5_beta")




## Plot mean slope and confidence intervals for all sites and years over curves from mixing model output
# Extract modeled lines from MixSIAR output for source proportions over time using a slightly modified plot function
pred_points <- lapply(models, function(x){
  load(paste0("MixSIAR_model_output/", x, "_com_snw.RData"))
  
  pred_points <- plot_continuous_var_2(jags.1, mix, source, output_options, exclude_sources_below=0)[[5]]
  
})
names(pred_points) <- models

pred.dat.1 <- do.call(rbind, pred_points)
pred.dat.1$Site <- rownames(pred.dat.1)
row.names(pred.dat.1) <- NULL
pred.dat <- pred.dat.1 %>% separate(Site, into = c("SY", NA), sep = '[.]') %>% 
  separate(SY, into = c("Site", "Year"), sep = '[_]') %>% na.omit() # Make columns for site and year


### Make 95% CI lines for calculated slopes
## Code for confidence intervals
# Create points for lines over the model window from calculated slopes and intercepts
x <- seq(38, 143, length.out = 100)
int <- slope.df2 %>% filter(Par == "alph")
slp <- slope.df2 %>% filter(Par == "beta")

calc <- function(S,s,yr, i){
  int2 <- int %>% filter(Site == S & Source == s & Year == yr) %>% select(Mean)
  slp2 <- slp %>% filter(Site == S & Source == s & Year == yr) %>% select(Mean)
  y = int2$Mean[i] + slp2$Mean[i]*x
  df <- data.frame(x=x, y=y, rep = rep(paste(i), length(y)), Site = rep(paste(S), length(y)),Source = rep(paste(s), length(y)), Year = rep(paste(yr), length(y)))
  return(df)
}

# Calculate predicted lines for all sites and years
args <- slope.df2 %>% select(Site, Source, Year) %>% unique()
preds <- mapply(calc, S = args$Site, s = args$Source, yr = args$Year, i = rep(1:30, length(args$Site)))
arr <- array2DF(preds)

x <- arr %>% filter(Var1 == "x") %>% select(Value)
y <- arr %>% filter(Var1 == "y") %>% select(Value)
Site <- arr %>% filter(Var1 == "Site") %>% select(Value)
Source <- arr %>% filter(Var1 == "Source") %>% select(Value)
Year <- arr %>% filter(Var1 == "Year") %>% select(Value)
rep <- arr %>% filter(Var1 == "rep") %>% select(Value)

pred <- data.frame(x=x, y=y, Site = Site, Source = Source, Year = Year, rep = rep)
names(pred) <- c("x", "y", "Site", "source", "Year", "rep")
pred$x <- as.numeric(pred$x)
pred$y <- as.numeric(pred$y)
head(pred)

# Calculate the median and CI for each x-value
pred.int <- pred %>% group_by(x, Site, source, Year) %>% summarise(upper = quantile(y, probs = .975),
                                                                   median = quantile(y, probs = .5),
                                                                   lower = quantile(y, probs = .025))


### Plots

# Indicate which slopes are not significantly different from 0
slope.df2 %>% filter(name %in% c("mu_slope", "mu_slope1", "mu_slope2")) %>% filter(p2.5*p97.5 < 0) # Select CI overlapping zero

pred.dat <- pred.dat %>% mutate(s = case_when(Site == "MOOS" & Year == "2021" & source == "soil_water" ~ "n",
                                    Site == "POKE" & Year == "2022" & source == "well" ~ "n",
                                    Site == "MOOS" & Year == "2022" & source == "well" ~ "n")) %>%
                                mutate(across(s, ~ifelse(is.na(s) == T, "y",.)))

# Plot of source proportions over time (Figure 6)
site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")


pred.dat %>% ggplot(aes(x, median, color = source)) +
  facet_grid(vars(Site), vars(factor(Year, levels=c("2022", "2021", "2020", "2019", "2018", "2015")))) +
  geom_line(aes(linetype = s), size = 0.8) + # Linetype indicates significance
  scale_linetype_manual(values = c(2,1)) +
  geom_ribbon(aes(ymin=low,ymax=high, fill = source), alpha=0.5, colour = NA) +
  scale_color_manual(labels = c("Precipitation", "Soil water", "Vault groundwater", "Groundwater"), 
                     values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb")) +
  scale_fill_manual(labels = c("Precipitation", "Soil water", "Vault groundwater", "Groundwater"), 
                    values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb")) +
  labs(color = NULL, y = "Proportion", x = "Days since snowmelt") + theme_bw() +
 theme(legend.position = "top", strip.background =element_rect(fill="white")) +
  guides(color = "none", fill = guide_legend(override.aes=list(alpha = 1))) +
  labs(fill = "Source") + guides(linetype = "none") + scale_x_continuous(n.breaks = 3) +
  theme(axis.text.x = element_text(hjust = 1))

ggsave("Figures/Pred lines.png", height =5, width = 6)


# Plot with fit linear slopes over mixing model curves (Figure S4)
site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")

pred.dat %>% ggplot(aes(x, median, color = source)) +
  facet_wrap(~ Site + Year, labeller = labeller(Site = site_labs), scales = "free_y") + 
  geom_ribbon(aes(ymin=low,ymax=high, fill = source), alpha=0.3, colour = NA) +
  scale_color_manual(labels = c("Precipitation", "Soil water", "Vault groundwater", "Groundwater"), 
                     values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb")) +
  scale_fill_manual(labels = c("Precipitation", "Soil water", "Vault groundwater", "Groundwater"), 
                     values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb")) +
  labs(color = NULL, y = "Proportion", x = "Days since snowmelt") + theme_bw() +
  geom_line(data = pred.int, aes(x,upper/100), stat = "smooth", linetype = 2) + # Divide by 100 because they were calculated as percentages
  geom_line(data = pred.int, aes(x,lower/100), stat = "smooth", linetype = 2) +
  geom_abline(data = slope.dat, aes(slope = Mean_beta/100, intercept = Mean_alpha/100, color = Source)) + 
  theme(legend.position = "top", strip.background =element_rect(fill="white"), text = element_text(size = 20)) +
  guides(color = "none", fill = guide_legend(override.aes=list(alpha = 1))) +
  labs(fill = "Source") + scale_x_continuous(breaks = c(60, 90, 120))

ggsave("Figures/Slope over pred.png", height =10, width = 10)






############### Model slopes with precipitation #############

# Join slopes with site characteristics
Site_char <- read.csv("Site_characteristics.csv") %>% pivot_longer(cols = 2:6, names_to = "Site") %>% pivot_wider(names_from = Attribute, values_from = value) %>% filter(Site != "CARI")
Site_char$Catchment_area <- Site_char$Catchment_area/100

dat.char.slope <- full_join(Site_char, slope.dat)

# Calculate cumulative precipitation
Precip.dat <- read.csv("DoD_precip_all.csv")
Precip.dat$datetimeAK <- as.POSIXct(Precip.dat$datetimeAK, tz = "America/Anchorage")
Precip.dat$Julian_day <- as.POSIXlt(Precip.dat$datetimeAK)$yday

# Calculate days since zero snow for precip data
Precip.dat <- Precip.dat %>% mutate(Days_melt = case_when(year == "2015" ~ Julian_day - 119,
                                                          year == "2018" ~ Julian_day - 130,
                                                          year == "2019" ~ Julian_day - 115,
                                                          year == "2020" ~ Julian_day - 130,
                                                          year == "2021" ~ Julian_day - 119,
                                                          year == "2022" ~ Julian_day - 131))



# Calculate cumulative precipitation between earliest common start date (35) and end of window (143)
Precip.cum.year <- Precip.dat %>% group_by(Site, year) %>% filter(Days_melt >= 35 & Days_melt <= 143) %>% 
  summarise(Precip = sum(inst_rainfall_mm, na.omit = T))

# Assign where to use precip from for all sites and years. At sites that did not have rain gages
# Eielson 2021 is much too high so we will have to use Stuart.

add_precip <- data.frame(Site = c("FRCH", "FRCH", "FRCH", rep("MOOS",5), "STRT", "STRT"), year = c("2018","2021","2022", "2018", "2019", "2020", "2021", "2022", "2019", "2022"), Precip = c(162, 188, 127, 162, 364, 357, 188, 127, 364, 127))
Precip.cum.year$year <- as.character(Precip.cum.year$year)
Precip.cum.year <- rbind(Precip.cum.year[-c(1:4),], add_precip)
Precip.cum.year$Precip <- Precip.cum.year$Precip/10
names(Precip.cum.year) <- c("Site", "Year", "Precip_cm")



# Calculate proportion fall precip
Precip.cum.fall <- Precip.dat %>% group_by(Site, year) %>% filter(Days_melt >= 90, Days_melt <= 143) %>% 
  summarise(Precip = sum(inst_rainfall_mm, na.omit = T))

add_precip.2 <- data.frame(Site = c("FRCH", "FRCH", "FRCH", rep("MOOS",5), "STRT", "STRT"), year = c("2018","2021","2022", "2018", "2019", "2020", "2021", "2022", "2019", "2022"), Precip = c(70,130,53.2,70,230,69,130,53.2,230,53.2        ))#77.4,121, 38.3, 77.4, 176, 78, 121, 38.3, 176, 38.3))
Precip.cum.fall$year <- as.character(Precip.cum.fall$year)
Precip.cum.fall <- rbind(Precip.cum.fall[-c(1:4),], add_precip.2)
Precip.cum.fall$Precip <- Precip.cum.fall$Precip/10
names(Precip.cum.fall) <- c("Site", "Year", "Precip_fall_cm")


## Join precipitation data
Precip.cum <- full_join(Precip.cum.fall, Precip.cum.year)

Precip.cum %>% ggplot(aes(Year, Precip_fall_cm, color = Site)) + geom_point() + ylab("Fall rain (cm)")

Precip.cum %>% ggplot(aes(Precip_cm, Precip_fall_cm, col = Year)) + geom_point() +
  facet_wrap(~Site)

# Proportion of summer precip in fall
Precip.cum$Prop_fall <- Precip.cum$Precip_fall_cm/Precip.cum$Precip_cm
dat.char.slope$Year <- as.character(dat.char.slope$Year)
slope.dat.char <- full_join(dat.char.slope, Precip.cum) %>% select(Site, Year, Source, Mean_beta, SD_beta, p2.5_beta, p97.5_beta, Precip_cm, Precip_fall_cm, Prop_fall, Catchment_area, Slope, Deciduous,Burned_largest, PF_percent)

write.csv(Precip.cum, "Output_data/Cumulative_precip_calc.csv")





#################### Models with SD
# JAGS model: Independent intercepts and slopes by site
## Model set up
lm_jags.1 <- function(){
  # Priors
  mu_int~dnorm(0, 2) # Mean hyperparameter for random intercepts
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
    mu[i] <- alpha[group[i]]+beta[group[i]]*vari[i]
    y[i]~dnorm(N[i], pow(Ysd[i], -2))
    N[i] ~ dnorm(mu[i], tau_res)
  }
}


# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "mu_slope", "sigma_slope", "sigma_res")

# Chains and iterations
n.chains = 3
n.iter = 5000 
n.burnin = 500


#### Run models the influence of cumulative or fall precipitation on groundwater or precipitation proportions.
## Organize slope data
# Slope data for precipitation proportion
PPT.slope <- slope.dat.char %>% filter(Source == "ppt/snow") %>% select(Site, Year, Source, Prop_fall, Precip_cm, Mean_beta, SD_beta)

PPT.slopes <- PPT.slope %>% mutate(across(Site, ~case_when(Site== "FRCH" ~ 1,
                                                           Site== "MOOS" ~ 2,
                                                           Site== "POKE" ~ 3,
                                                           Site== "STRT" ~ 4,
                                                           Site== "VAUL" ~ 5)))

PPT.slopes %>% ggplot(aes(Precip_cm, Mean_beta, col = as.character(Site))) + geom_point()

# Slope data for groundwater proportions (without Vault)
well.slope <- slope.dat.char %>% filter(Source %in% c("well")) %>% select(Site, Year, Source, Prop_fall, Precip_cm, Mean_beta, SD_beta)


well.slopes <- well.slope %>% mutate(across(Site, ~case_when(Site== "FRCH" ~ 1,
                                                             Site== "MOOS" ~ 2,
                                                             Site== "POKE" ~ 3,
                                                             Site== "STRT" ~ 4)))

## Run models
# Cumulative summer precipitation and precipitation proportion
Nsites <- 5
jagsdata <- with(PPT.slopes, list(y = Mean_beta, Ysd = SD_beta, vari = Precip_cm, group = Site,
                                  n = length(Mean_beta), Nsites = Nsites))
fit_slope.ppt <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.1,
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.ppt

# Cumulative fall precipitation and precipitation proportion
jagsdata <- with(PPT.slopes, list(y = Mean_beta, Ysd = SD_beta, vari = Prop_fall, group = Site,
                                  n = length(Mean_beta), Nsites = Nsites))
fit_slope.ppt.fall <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.1,
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.ppt.fall


# Cumulative fall precipitation and groundwater proportion
Nsites <- 4
jagsdata <- with(well.slopes, list(y = Mean_beta, Ysd = SD_beta, vari = Precip_cm, group = Site,
                                   n = length(Mean_beta), Nsites = Nsites))
fit_slope.well.cum <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.1,
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.well.cum

summary(lm(well.slopes[well.slopes$Site == 2,]$Mean_beta ~ well.slopes[well.slopes$Site == 2,]$Precip_cm))
plot(well.slopes[well.slopes$Site == 2,]$Mean_beta ~ well.slopes[well.slopes$Site == 2,]$Precip_cm)

# Cumulative fall precipitation and groundwater proportion
Nsites <- 4
jagsdata <- with(well.slopes, list(y = Mean_beta, Ysd = SD_beta, vari = Prop_fall, group = Site,
                                   n = length(Mean_beta), Nsites = Nsites))
fit_slope.well.fall <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.1,
                            n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.well.fall





### Model for just Vault for relationship between Vault groundwater proportion and cumulative/fall precipitation
# Data
up.slope <- slope.dat.char %>% filter(Source %in% c("upwelling")) %>% select(Site, Year, Source, Prop_fall, Precip_cm, Mean_beta, SD_beta)
up.slopes <- up.slope %>% mutate(across(Site, ~case_when(Site== "VAUL" ~ 1)))

# JAGS model: one site
lm_jags.2 <- function(){
  # Priors
  alpha ~dnorm(0, 1)
  beta ~ dnorm(0, 2)
  sigma~dunif(0, 100) # Residual standard deviation
  tau <- 1/(sigma*sigma) # Residual precision
  
  # Likelihood
  for (i in 1:n) {
    mu[i] <- alpha + beta*vari[i]
    y[i]~dnorm(N[i], pow(Ysd[i], -2))
    N[i] ~ dnorm(mu[i], tau)
  }
}

# Parameters to estimate
params <- c("alpha", "beta", "sigma")


# Cumulative fall precipitation and Vault groundwater proportion

jagsdata <- with(up.slopes, list(y = Mean_beta, Ysd = SD_beta, vari = Prop_fall,
                                 n = length(Mean_beta)))
fit_slope.upwelling.fall <- jags(data = jagsdata, parameters.to.save = params, model.file = lm_jags.2,
                                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin)
fit_slope.upwelling.fall




########## Plot
#### Plot models with significant results: relationship between fall precipitation and groundwater proportion
# Extract draws for intercepts and slopes for CI
mod.sum <- as.data.frame(fit_slope.well.fall$BUGSoutput$summary)
mod.sum <- mod.sum[1:8, ]
mod.sum$Site <- rep(c("FRCH", "MOOS", "POKE", "STRT"), 2)
mod.sum$param <- c(rep("alpha", 4), rep("beta", 4))
mod.sum <- mod.sum %>% select(mean, param, Site) %>% pivot_wider(values_from = mean, names_from = param)


out <- as.data.frame(as.matrix(as.mcmc(fit_slope.well.fall)))
names(out) <- c("a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4", "deviance", "sigma_res")

# Endpoints for segments
ends <- slope.dat.char %>% filter(Source %in% c("well")) %>% select(Site,Prop_fall) %>% group_by(Site) %>% 
  summarise(p1 = min(Prop_fall), p2 = max(Prop_fall))
ends
mod.sum
mod.sum.ends <- full_join(ends, mod.sum)
mod.sum.ends$y1 <- mod.sum.ends$alpha+mod.sum.ends$beta*mod.sum.ends$p1
mod.sum.ends$y2 <- mod.sum.ends$alpha+mod.sum.ends$beta*mod.sum.ends$p2
mod.sum.ends


## Code for confidence intervals
# Create points from slopes and intercepts
nrow <- nrow(out$a1)
x <- seq(0,.90,.01)
ints <- out[,1:4]
slps <- out[,5:8]
Site <- 1:4

calc <- function(s, i){
  y = ints[s][i,] + slps[s][i,]*x
  df <- data.frame(x=x, y=y, rep = rep(paste(i), length(y)), Site = rep(paste(s), length(y)))
  return(df)
}
calc(4,1)
args <- expand.grid( b = 1:4, a = 1:length(out$a1))
preds <- mapply(calc, s = args$b, i = args$a)
arr <- array2DF(preds)

x <- arr %>% filter(Var1 == "x") %>% select(Value)
y <- arr %>% filter(Var1 == "y") %>% select(Value)
Site <- arr %>% filter(Var1 == "Site") %>% select(Value)
rep <- arr %>% filter(Var1 == "rep") %>% select(Value)
pred <- data.frame(x=x, y=y, Site = Site, rep = rep)
names(pred) <- c("x", "y", "Site", "rep")
pred$x <- as.numeric(pred$x)
pred$y <- as.numeric(pred$y)
head(pred)


pred.int <- pred %>% group_by(x, Site) %>% summarise(upper = quantile(y, probs = .975),
                                                     lower = quantile(y, probs = .025))
pred.int$Source <- "well"

pred.int <- pred.int %>% mutate(across(Site, ~case_when(Site== "1" ~ "FRCH",
                                                        Site== "2" ~ "MOOS",
                                                        Site== "3" ~ "POKE",
                                                        Site== "4" ~ "STRT")))

# Bounds 
ends
pred.int <- pred.int %>% filter(x >= 0.193 & x <= 0.804,
                                x >= 0.308 | Site %in% c("MOOS", "FRCH"),
                                x >= 0.313 | Site %in% c("MOOS", "FRCH", "POKE"),
                                x <= 0.691 | Site %in% c("POKE"))


# Plot of all sites
mod.sum.ends <- mod.sum.ends %>% filter(Site %in% c("FRCH",  "STRT"))
pred.int <- pred.int %>% filter(Site %in% c("FRCH", "STRT"))


df2 <- data.frame(Site = rep(c("FRCH", "MOOS", "POKE", "STRT", "VAUL"),2), Prop_fall = rep(.5, 10), Mean_beta = c(rep(0.15,4), 0, rep(-0.25,4), -1.5), Source = rep("well",10))

site_labs <- c("French", "Moose", "Poker", "Stuart", "Vault")
names(site_labs) <- c("FRCH", "MOOS", "POKE", "STRT", "VAUL")


slope.dat.char   %>% filter(is.na(Source) ==F, Source %in% c("upwelling", "well")) %>%
  ggplot() + geom_point(aes(Prop_fall, Mean_beta, color = Source)) +
  facet_wrap(~Site, scales = "free_y", labeller = labeller(Site = site_labs)) + 
  xlab("Proportion of rain in fall") + ylab("Groundwater seasonal slope (%/day)") +
  theme_bw() + theme(strip.background =element_rect(fill="white")) +
  geom_segment(data = mod.sum.ends, aes(x = p1, xend = p2, y = y1, yend = y2), color = "#8da0cb") +
  geom_point(data = df2, aes(Prop_fall, Mean_beta), alpha = 0) + # fix size of plot
  scale_color_manual(breaks = c("well", "upwelling"), labels = c("Groundwater", "Vault groundwater"), 
                     values = c("#8da0cb", "#e78ac3")) +
  theme(strip.background =element_rect(fill="white"), legend.position = c(0.85,0.22),
        panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  labs(color = "Source") +
  geom_ribbon(data = pred.int, aes(x = x, ymin=lower,ymax=upper), alpha=0.3, fill = "#8da0cb") +
  geom_errorbar(aes(Prop_fall, Mean_beta, ymin=Mean_beta-SD_beta, ymax=Mean_beta+SD_beta, color = Source), width=.1)


ggsave("Figures/Groundwater slope and fall rain SD.png", width = 4.5, height = 2.9)







