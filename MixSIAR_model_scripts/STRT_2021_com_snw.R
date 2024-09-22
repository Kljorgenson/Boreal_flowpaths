### DoD Mixing Models
library(MixSIAR)
library(dplyr)
library(ggplot2)
library(lubridate)

all_dat <- read.csv("mixing2/All_mixing_data_Fox.csv")

### MixSIAR models 4 sites
k <- "STRT_2021"
mix_d <- all_dat %>% filter(Site == "STRT", year == "2021") %>% select(Days_melt, Chloride_uM, Magnesium_uM) %>% na.omit()



source_dat <- all_dat %>% filter(Type %in% c("well", "ppt/snow", "soil_water")) %>% select(Type, Chloride_uM, Magnesium_uM)

source_dat$Type <- as.factor(source_dat$Type)
source_dat_means <- source_dat %>% group_by(Type) %>% summarise(n = length(Chloride_uM),
                                                                MeanChloride_uM=mean(Chloride_uM, na.rm = TRUE),
                                                                MeanMagnesium_uM = mean(Magnesium_uM, na.rm = TRUE),
                                                                SDChloride_uM =sd(Chloride_uM, na.rm = TRUE),
                                                                SDMagnesium_uM = sd(Magnesium_uM, na.rm = TRUE))



write.csv(mix_d,paste("mixing2/mix_dat_", k, ".csv", sep = ""), 
          row.names = FALSE)
write.csv(source_dat_means,paste("mixing2/source_dat_", k, ".csv", sep = ""), 
          row.names = FALSE)

# Load mixture data
mix <- load_mix_data(filename=paste("mixing2/mix_dat_", k, ".csv", sep = ""),
                     iso_names=c("Chloride_uM", "Magnesium_uM"),
                     factors= NULL,
                     fac_random= NULL,
                     fac_nested= NULL,
                     cont_effects="Days_melt")


# Load source data

source <- load_source_data(filename=paste("mixing2/source_dat_", k, ".csv", sep = ""),
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)

# Discrimination factors
discr <- as.data.frame(matrix(rep(0, 2*length(source$source_names)*(length(names(source_dat))-1)), 
                              nrow =length(source$source_names))) 
row.names(discr) <- source$source_names
names(discr) <- c(paste("Mean", names(source_dat)[2], sep = ""),paste("Mean", names(source_dat)[3], sep = ""),#paste("Mean", names(source_dat)[4], sep = ""),#,paste("Mean", names(source_dat)[5], sep = ""),
                  paste("SD", names(source_dat)[2], sep = ""),paste("SD", names(source_dat)[3], sep = ""))#,paste("SD", names(source_dat)[4], sep = ""))#, paste("SD", names(source_dat)[5], sep = ""))

write.csv(discr,paste("mixing2/discr_", k, ".csv", sep = ""), 
          row.names = TRUE)
discr <-load_discr_data(filename=paste("mixing2/discr_", k, ".csv", sep = ""), 
                        mix)
discr


# Write the JAGS model file
model_filename <- paste("mixing2/MixSIAR_model_", k, sep = "")   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE # Water may come from one section of the source distribution
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run model
jags.1 <- run_model(run="extreme", mix, source, discr, model_filename,  
                    alpha.prior = 1, resid_err, process_err)


save.image(file = paste("mixing2/",k,"_com_snw.RData", sep = ""))


