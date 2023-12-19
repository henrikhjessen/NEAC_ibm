# Summary --------------------------------------------------------------------------------------
# This is the script for generating plots, to use in my second paper on how climate and fisheries interact to cause evolution
# Plots are organised into subheading for ease of navigating this script. 
# Data import takes 30+ minutes, recommend saving workspace when exiting.
# Plotting dimensions recommended:
#
# PowerPoint: width = 25.4, height = 14.3, unit = "cm". For PowerPoint, use .svg format.
# LaTeX:  width(cm) = 19(full page), 14(1,5 column), 90(single column), 3(minimum)
#         Use .pdf format for best scaling, text size 12 or 13, or 500dpi.
#
#
# Set up the workspace -------------------------------------------------------------------------
library(here)
library(devtools)
library(tidyverse)
library(tidyr)
library(gganimate)
library(readr)
library(dplyr)
library(tibble)
library(reshape2)
library(scales)
library(RNetCDF)
library(abind)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(shiny)
library(esquisse)
library(Hmisc)
library(gifski)
library(transformr)
library(tictoc)
library(paletteer) #https://r-charts.com/color-palettes/

## Define new functions ------------------------------------------------------------------------------------
## N: vector of sizes/weights (number of fish)
## M: vector of means (means after summarizing repetitions)
## S: vector of standard deviations (standard deviations after summarizing repetitions)

#This functions is used to combine standard deviations in the two-step process
#Based on https://www.brainkart.com/article/Combined-Mean-and-Combined-Standard-Deviation_35098/
combine_sd <- function (M = NULL, S = NULL, N = NULL) 
{
  combined_mean <- sum( N * M ) / sum( N )
  mean_dev <- combined_mean - M
  combined_var <- sum( N * ( S^2 + mean_dev^2 ) ) / sum( N )
  combined_sd <- sqrt(combined_var)
  return(combined_sd)
}


my_palette <- c("#000000","#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7")
temp_palette <- c("#25CED1","#697684","#864559","#A3142E")
big_palette <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", 20)
raster_palette <- "ggthemes::Green"

## Load data ------------------------------------------------------------------------------------
#Based on netCDF
#Import the data files
tic()
data_files <- list.files(path = here("output_paper3") , pattern = "*.nc") #create list of data files
setwd(here("output_paper3"))
for (i in 1:length(data_files)) { #Read in the data files
  ncdf_temp <- open.nc(data_files[i])
  print(data_files[i])
  #print.nc(ncdf_temp)
  data_holder2 <- read.nc(ncdf_temp)
  data_holder <- data_holder2$PopResults
  if (i == 1) {
    ncdf_dataframe <- data_holder
  } else {
    ncdf_dataframe <- abind(ncdf_dataframe, data_holder, along = 4)
  }
  parameter_holder <- data_holder2$Parameters
  if (i == 1) {
    parameters <- parameter_holder
  } else {
    parameters <- cbind(parameters, parameter_holder)
  }
  heritability_holder <- data_holder2$Heritabilities
  if (i == 1) {
    heritabilities <- heritability_holder
  } else {
    heritabilities <- abind(heritabilities, heritability_holder,along = 1)
  }
  close.nc(ncdf_temp)
  rm(ncdf_temp, data_holder2, data_holder, parameter_holder, heritability_holder)
}
setwd(here())
parameter_names <- c("Horizon","MaxAge","TargetInd","MaxInd","k","EnDensGon","EnDensSom","HarvestDuration","HarvestStart",
                     "beta","c_R","c_phi","c_SMR","c_SDA","D_M","c_COT","c_u","b_2","b_3","b_4","v1","v2","v3","v4","InvestDecay",
                     "c_predation","c_foraging","c_respiration","M_fixed","pred_exp","for_risk","repro_exp","resp_exp","tau",
                     "eta","EggWeight","ClimScen","GearType","PropTrawl","Fmax","Lmax","FlatMort","FishSelWidth",
                     "density_dependance", "Psi", "Theta","phenotypic_deviance","inheritance_deviance","FoodEnv_deviance","temp_deviance",
                     "recruitment_deviance")
if (length(data_files)>1){colnames(parameters)[1:length(data_files)] <- seq(1,length(data_files),1)}
parameters <- cbind(parameter_names, parameters)
rm(i, parameter_names) #cleanup

array_runs <- length(data_files)
if (length(dim(ncdf_dataframe)) == 4) { #Only do this chunk for multiple runs
  ncdf_dataframe <- aperm(ncdf_dataframe, c(4,1,2,3)) #Place "run" dimension first
  #Create a column with run-number to splice onto the dataframe
  yearsperrun <- dim(ncdf_dataframe)[2]
  array_agegroups <- dim(ncdf_dataframe)[3]
  run_col <- vector(mode='numeric', length = array_runs*yearsperrun*array_agegroups) #Allocate the column
  for (i in 1:array_runs) { #Fill run_col with numbers according to run
    start <- (1+ (i-1)*yearsperrun*array_agegroups)
    stop <- (yearsperrun*array_agegroups+ (i-1)*yearsperrun*array_agegroups)
    run_col[start:stop] = i
    rm(start, stop)
  }
  rm(i) #cleanup
  for (i in 1:array_runs) { #stack runs on top of each other in the "year" column
    data_holder <- ncdf_dataframe[i,,,]
    data_holder_length <- length(data_holder[1])
    if (i == 1) {
      temp_PopResults <- data_holder
    } else {
      temp_PopResults <- abind(temp_PopResults, data_holder, along = 1)
    }
    rm(data_holder, data_holder_length)
  }
  rm(i) #cleanup
  #Do the same for heritabilities
  run_col2 <- vector(mode='numeric', length = array_runs*yearsperrun) #Allocate the column
  for (i in 1:array_runs) { #Fill run_col with numbers according to run
    start <- (1+ (i-1)*yearsperrun)
    stop <- (yearsperrun+ (i-1)*yearsperrun)
    run_col2[start:stop] = i
    rm(start, stop)
  }
  rm(yearsperrun, array_agegroups, i)
} else {
  temp_PopResults <- ncdf_dataframe
  run_col <- rep(1,dim(ncdf_dataframe)[1]*dim(ncdf_dataframe)[2])
  run_col2 <- rep(1,dim(ncdf_dataframe)[1])
}

array_years <- dim(temp_PopResults)[1]
for (i in 1:array_years) { #Stack years on top of each other in the age-group column, and splice the run_col onto this
  data_holder <- temp_PopResults[i,,]
  data_holder_length <- nrow(data_holder)
  if (i == 1) {
    PopResults <- data_holder
  } else {
    PopResults <- rbind(PopResults, data_holder)
  }
  if (i == array_years) {
    PopResults <- cbind(run_col, PopResults)
    heritabilities <- cbind(run_col2, heritabilities)
    rm(run_col, run_col2)
  } 
}
PopResults <- as.data.frame(PopResults) #convert to data frame
heritabilities <- as.data.frame(heritabilities) #convert to data frame
heritabilities$run <- as.factor(heritabilities$run)
parameters <- as.data.frame(parameters) #convert to data frame
colnames(PopResults)[1:48] <-c("run","year","age","temp","FoodEnv","n_recr","n_survived","n_died","SomWeight","SomWeight_SD",
                               "GonWeight","GonWeight_SD","length","length_SD","M_predation","Mpred_SD","M_foraging",
                               "Mfor_SD","M_reproduction","Mrepr_SD","M_respiration","Mresp_SD","FishMort","FishMort_SD",
                               "f_int","fint_SD","n_mature","n_spawnskip","pop_biomass","som_growth", "somgro_SD",
                               "gon_growth", "gongro_SD","appetite","appetite_SD","som_allocation","allocation_SD",
                               "intercept","intercept_SD","ageatmat","ageatmat_SD","lenatmat","lenatmat_SD","CumGonWeight",
                               "CumGonWeight_SD","n_fished","yield","SSB") #name columns
colnames(heritabilities)[1:4] <- c("run","year","ageatmat","lengthatmat")
PopResults$run <- as.factor(PopResults$run)
n <- PopResults$n_survived + PopResults$n_died
PopResults <- add_column(PopResults, n, .after = "n_died")
pct_mature <- (PopResults$n_mature / PopResults$n) * 100
PopResults <- add_column(PopResults, pct_mature, .after = "n_mature")
pct_skippers <- (PopResults$n_spawnskip / PopResults$n) * 100
PopResults <- add_column(PopResults, pct_skippers, .after = "n_spawnskip")
tot_growth <- (PopResults$som_growth + PopResults$gon_growth)
PopResults <- add_column(PopResults, tot_growth, .after = "gongro_SD")
gon_allocation <- as.vector(1-PopResults$som_allocation)
PopResults <- add_column(PopResults, gon_allocation, .after = "som_allocation")
Z <- PopResults$M_predation +PopResults$M_foraging+PopResults$M_reproduction+PopResults$M_respiration+PopResults$FishMort + 0.07
PopResults <- add_column(PopResults, Z, .after = "FishMort_SD")
PopResults$NatMort <- 0.07 + PopResults$M_foraging + PopResults$M_predation + PopResults$M_reproduction + PopResults$M_respiration
rm(data_holder, array_years,array_runs, data_holder_length, i, n, pct_mature, pct_skippers, 
   tot_growth, gon_allocation, temp_PopResults, Z) #cleanup

PopResults <- PopResults %>% drop_na(year)
heritabilities <- subset(heritabilities, select = -c(run))
rm(ncdf_dataframe, data_files)
toc()


## Name the runs ----------------------------------------------------------------------------------------------------------
run_names_holder1 <- subset(parameters, parameter_names=="GearType")
run_names_holder2 <- subset(parameters, parameter_names=="ClimScen")
run_names_holder3 <- subset(parameters, parameter_names=="Fmax")
run_names_holder4 <- subset(parameters, parameter_names=="PropTrawl")
run_names_holder5 <- subset(parameters, parameter_names=="Lmax")
run_names1 <- vector(mode="character",length = ncol(run_names_holder1)-1)
run_names2 <- vector(mode="character",length = ncol(run_names_holder2)-1)
run_names3 <- vector(mode="character",length = ncol(run_names_holder3)-1)
run_names4 <- vector(mode="character",length = ncol(run_names_holder4)-1)
run_names5 <- vector(mode="character",length = ncol(run_names_holder5)-1)


for (i in 1:(ncol(run_names_holder1)-1)){
  run_names1[i] <- run_names_holder1[1,i+1]
}
for (i in 1:(ncol(run_names_holder2)-1)){
  run_names2[i] <- run_names_holder2[1,i+1]
}
for (i in 1:(ncol(run_names_holder3)-1)){
  run_names3[i] <- run_names_holder3[1,i+1]
}
for (i in 1:(ncol(run_names_holder4)-1)){
  run_names4[i] <- run_names_holder4[1,i+1]
}
for (i in 1:(ncol(run_names_holder5)-1)){
  run_names5[i] <- run_names_holder5[1,i+1]
}
run_names3 <- as.double(run_names3) %>%
  round(digits = 2) %>%
  as.character()

PopResults$GearType <- as.character(PopResults$run)
PopResults$ClimScen <- as.character(PopResults$run)
PopResults$Fmax <- as.character(PopResults$run)
PopResults$PropTrawl <- as.character(PopResults$run)
PopResults$Lmax <- as.character(PopResults$run)

for (i in 1:length(run_names1)) {
  PopResults$GearType[PopResults$GearType == i] <- run_names1[i]
}
for (i in 1:length(run_names2)) {
  PopResults$ClimScen[PopResults$ClimScen == i] <- run_names2[i]
}
for (i in 1:length(run_names3)) {
  PopResults$Fmax[PopResults$Fmax == i] <- run_names3[i]
}
for (i in 1:length(run_names3)) {
  PopResults$PropTrawl[PopResults$PropTrawl == i] <- run_names4[i]
}
for (i in 1:length(run_names3)) {
  PopResults$Lmax[PopResults$Lmax == i] <- run_names5[i]
}

rm(run_names_holder1, run_names_holder2, run_names_holder3,run_names_holder4,
   run_names_holder5,i, run_names1, run_names2, run_names3, run_names4, run_names5)

PopResults$PropTrawl <- PopResults$PropTrawl %>%
  as.double() %>% round(digits=2) %>% as.character()

popresults_copy <- PopResults
## Summarize the data -----------------------------------------------------------------------------------------------------
PopResults$GearType[PopResults$GearType == "0"] <- "Unfished"
PopResults$GearType[PopResults$GearType == "4"] <- "Fished"
PopResults$Fmax[PopResults$GearType == "Unfished"] <- "0.0"
PopResults$ClimScen[PopResults$ClimScen == "0"] <- "No warming"
PopResults$ClimScen[PopResults$ClimScen == "1"] <- "SSP1"
PopResults$ClimScen[PopResults$ClimScen == "2"] <- "SSP2"
PopResults$ClimScen[PopResults$ClimScen == "3"] <- "SSP3"

PopResults$PropTrawl[PopResults$PropTrawl == "0"] <- "0%"
PopResults$PropTrawl[PopResults$PropTrawl == "0.3"] <- "30%"
PopResults$PropTrawl[PopResults$PropTrawl == "0.5"] <- "50%"
PopResults$PropTrawl[PopResults$PropTrawl == "0.7"] <- "70%"
PopResults$PropTrawl[PopResults$PropTrawl == "1"] <- "100%"

PopResults$ClimScen[PopResults$ClimScen == "SSP1"] <- "SSP1 (+0.8 °C)"
PopResults$ClimScen[PopResults$ClimScen == "SSP2"] <- "SSP2 (+3.0 °C)"
PopResults$ClimScen[PopResults$ClimScen == "SSP3"] <- "SSP3 (+8.4 °C)"



age_data <- PopResults %>%
  group_by(GearType, ClimScen, Fmax,PropTrawl,Lmax, year, age) %>%
  summarise(temperature = mean(temp),
            avg_appetite = wtd.mean(appetite, n), sd_app = sqrt(wtd.var(appetite, n)), se_app = sqrt(wtd.var(appetite, n))/sqrt(20), #Appetite
            avg_intercept = wtd.mean(intercept, n), sd_int = sqrt(wtd.var(intercept, n)),se_int = sqrt(wtd.var(intercept, n))/sqrt(20), #Intercept
            avg_gonall = wtd.mean(gon_allocation, n), sd_gonall = sqrt(wtd.var(gon_allocation, n)), se_gonall = sqrt(wtd.var(gon_allocation, n))/sqrt(20), #Allocation
            avg_lenatmat = wtd.mean(lenatmat, n_mature), sd_lenatmat= sqrt(wtd.var(lenatmat, n_mature)), #Lenatmat
            avg_ageatmat = wtd.mean(ageatmat, n_mature), sd_ageatmat= sqrt(wtd.var(ageatmat, n_mature)), #Ageatmat
            avg_n = mean(n), sd_n = sd(n), avg_mature = mean(n_mature), sd_mature = sd(n_mature),         #Population dynamics
            avg_recr = mean(n_recr), sd_recr = sd(n_recr), avg_died = mean(n_died), sd_died = sd(n_died), #Population dynamics
            avg_SomWeight = mean(SomWeight), sd_SomWeight = sd(SomWeight), #Somatic weight
            avg_GonWeight = mean(GonWeight), sd_GonWeight = sd(GonWeight), #Gonad weight
            avg_length = mean(length), sd_length=sd(length), #Length
            avg_SomGro = mean(som_growth), sd_SomGro = sd(som_growth), #Somatic growth
            avg_GonGro = mean(gon_growth), sd_GonGro = sd(gon_growth), #Gonadic growth
            avg_Mpred = mean(M_predation), sd_Mpred = sd(M_predation), #Predation mortality
            avg_Mfor = mean(M_foraging), sd_Mfor = sd(M_foraging), #Foraging mortality
            avg_Mrepro = mean(M_reproduction), sd_Mrepro = sd(M_reproduction), #Reproductive mortality
            avg_Mresp = mean(M_respiration), sd_Mresp = sd(M_respiration),#Respiration mortality
            avg_Fishmort = mean(FishMort), sd_Fishmort = sd(FishMort),#Fisheries mortality
            avg_Z = mean(Z), sd_Z = sd(Z), #Total mortality
            avg_somgro = mean(som_growth), sd_somgro = sd(som_growth), #Somatic growth
            avg_gongro = mean(gon_growth), sd_gongro = sd(gon_growth),  #Gonadic growth
            avg_yield = mean(yield), sd_yield = sd(yield), #Yield
            avg_caught = mean(n_fished), sd_caught = sd(n_fished) #Number caught
  )
age_data$avg_appetite <- age_data$avg_appetite/1000
age_data$sd_app <- age_data$sd_app/1000
age_data$avg_gonall <- age_data$avg_gonall*100
age_data$sd_gonall <- age_data$sd_gonall*100
age_data <- subset(age_data, year > 0)
age_data$Lmax <- relevel(age_data$Lmax, ref=c(5,1,2,3,4))

sum_data <- PopResults %>%
  group_by(GearType, ClimScen, Fmax,PropTrawl,Lmax, year) %>%
  summarise(temperature = mean(temp),sd_temperature = sd(temp),
            avg_appetite = wtd.mean(appetite, n), sd_app = sqrt(wtd.var(appetite, n)), se_app = sqrt(wtd.var(appetite, n))/sqrt(20), #Appetite
            avg_intercept = wtd.mean(intercept, n), sd_int = sqrt(wtd.var(intercept, n)),se_int = sqrt(wtd.var(intercept, n))/sqrt(20), #Intercept
            avg_gonall = wtd.mean(gon_allocation, n), sd_gonall = sqrt(wtd.var(gon_allocation, n)), se_gonall = sqrt(wtd.var(gon_allocation, n))/sqrt(20), #Allocation
            avg_lenatmat = wtd.mean(lenatmat, n_mature), sd_lenatmat= sqrt(wtd.var(lenatmat, n_mature)), #Lenatmat
            avg_ageatmat = wtd.mean(ageatmat, n_mature), sd_ageatmat= sqrt(wtd.var(ageatmat, n_mature)) #Ageatmat
  )
sum_data$avg_appetite <- sum_data$avg_appetite/1000
sum_data$sd_app <- sum_data$sd_app/1000
sum_data$avg_gonall <- sum_data$avg_gonall*100
sum_data$sd_gonall <- sum_data$sd_gonall*100
sum_data <- subset(sum_data, year > 0)

gc()

## Evolving traits plots ------------------------------------------------------------------------------------------------
evo_data <- subset(sum_data,GearType=="Fished" & year =="2000")
frames <- subset(evo_data, Lmax == "110" & PropTrawl == "70%")

evo1 <- ggplot(subset(sum_data,GearType=="Fished")) +
  aes(x = year, y = avg_appetite, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_appetite - sd_app, ymax = avg_appetite + sd_app, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = bquote("Appetite "(KJ~kg^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_grid(factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme_bw() +
  theme(text = element_text(size = 11), legend.position = "top")
evo1
ggsave(here("R/Rplots3/evolving_appetite.pdf"), plot=evo1, width = 20, height = 22, units = "cm")
evo2 <- ggplot(subset(sum_data,GearType=="Fished")) +
  aes(x = year, y = avg_intercept, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "PMRN Intercept (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  facet_grid(factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme_bw() +
  theme(text = element_text(size = 11), legend.position = "top")
evo2
ggsave(here("R/Rplots3/evolving_intercept.pdf"), plot=evo2, width = 20, height = 22, units = "cm")
evo_arr1 <- ggarrange(evo1,evo2, ncol=1, nrow=2, common.legend = TRUE, legend = "top")
ggsave(here("R/Rplots3/evolving_traits.pdf"), plot=evo_arr1, width = 20, height = 40, units = "cm")

evo3 <- ggplot(data=evo_data) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_appetite)) +
  geom_tile(data=frames,aes(fill=avg_appetite),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average\nappetite (J)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "right")
evo4 <- ggplot(data=evo_data) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_intercept)) +
  geom_tile(data=frames,aes(fill=avg_intercept),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average PMRN\nintercept (cm)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "right")

evo_arr2 <- ggarrange(evo3,evo4, ncol=1, nrow=2, common.legend = F, legend = "right", align = "v")
ggsave(here("R/Rplots3/evo_traits_raster.pdf"), plot=evo_arr2, width = 20, height = 10, units = "cm")

rm(evo1,evo2, evo_arr1)
gc()

## Population dynamics ----------------------------------------------------------------------------------------------
popdynamics_data <- subset(age_data,GearType=="Fished") %>%
  group_by(GearType, ClimScen, Fmax,PropTrawl,Lmax, year) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            n_mat = sum(avg_mature), mat_sd = combine_sd(avg_mature, sd_mature, avg_n),
            pop_biomass = sum(avg_SomWeight*avg_n)/1000, biomass_sd = sum(sd_SomWeight * avg_n)/1000
  )
popdynamics_data <- subset(popdynamics_data, year >0)
popdynamics_data$rel_n <- popdynamics_data$n/121590.00
popdynamics_data$rel_biomass <- popdynamics_data$pop_biomass/296.3876
popdynamics_data$rel_n_sd <- popdynamics_data$n_sd/121590.00
popdynamics_data$rel_biomass_sd <- popdynamics_data$biomass_sd/296.3876
raster_dynamics <- subset(popdynamics_data, year == 2000)
frames <- subset(raster_dynamics, Lmax == "110" & PropTrawl == "70%")
popdynamics_data$n_sd <- popdynamics_data$n_sd/1000
popdynamics_data$n <- popdynamics_data$n/1000

pop3 <- ggplot(popdynamics_data) +
  aes(x = year, y = n, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=n - n_sd, ymax = n + n_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average number of individuals (thousands)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_grid(factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme(text = element_text(size = 10), legend.position = "top")
ggsave(here("R/Rplots3/total_individuals.pdf"), plot=pop3, width = 20, height = 22, units = "cm")

pop4 <- ggplot(popdynamics_data) +
  aes(x = year, y = pop_biomass, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin= pop_biomass - biomass_sd, ymax = pop_biomass + biomass_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average population biomass (t)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_grid(factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10), legend.position = "top")
ggsave(here("R/Rplots3/total_biomass.pdf"), plot=pop4, width = 20, height = 22, units = "cm")

pop_arr1 <- ggarrange(pop3,pop4, ncol=1, nrow=2, common.legend = TRUE, legend = "top")
ggsave(here("R/Rplots3/evolving_popdynamics.pdf"), plot=pop_arr1, width = 20, height = 30, units = "cm")

#Raster plots ----
pop5 <- ggplot(data=raster_dynamics) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=rel_n)) +
  geom_tile(data=frames,aes(fill=rel_n),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Relative number of individuals") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 12), legend.position = "top")

pop6 <- ggplot(data=raster_dynamics) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=rel_biomass)) +
  geom_tile(data=frames,aes(fill=rel_biomass),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Relative population biomass") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 12), legend.position = "top")

pop7 <- ggplot(data=raster_dynamics) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=pop_biomass)) +
  geom_tile(data=frames,aes(fill=pop_biomass),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average population\nbiomass (t)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "right")
pop8 <- ggplot(data=raster_dynamics) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=n)) +
  geom_tile(data=frames,aes(fill=n),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average number\nof individuals") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "right")

#pop_arr1 <- pop5/pop6
pop_arr1 <- ggarrange(pop7,pop8, ncol=1, nrow=2, common.legend = F, legend = "right", align = "v")

ggsave(here("R/Rplots3/popdynamics_raster.pdf"), plot=pop_arr1, width = 20, height = 10, units = "cm")

aov_biomass <- aov(data=popdynamics_data, pop_biomass~ClimScen+PropTrawl+Lmax+
                PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_biomass)

sink(here("R/rplots3/anova_biomass.txt"))
capture.output(summary(aov_biomass, file= here("R/rplots3/anova_biomass.txt")))
sink()
aov_individuals <- aov(data=popdynamics_data, n~ClimScen+PropTrawl+Lmax+
                     PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_individuals)

sink(here("R/rplots3/anova_individuals.txt"))
capture.output(summary(aov_individuals, file= here("R/rplots3/anova_individuals.txt")))
sink()

rm(popdynamics_data,raster_dynamics,frames, pop3,pop4,pop5,pop6,pop_arr1, pop7,pop8) 
gc()



## Life history ----------------------------------------------------------------------------------------------
lh_dat2 <- subset(sum_data, year==2000 & GearType=="Fished")
lh_dat2$PropTrawl[lh_dat2$PropTrawl == "0"] <- "0%"
lh_dat2$PropTrawl[lh_dat2$PropTrawl == "0.3"] <- "30%"
lh_dat2$PropTrawl[lh_dat2$PropTrawl == "0.5"] <- "50%"
lh_dat2$PropTrawl[lh_dat2$PropTrawl == "0.7"] <- "70%"
lh_dat2$PropTrawl[lh_dat2$PropTrawl == "1"] <- "100%"
frames = subset(lh_dat2, Lmax == "110" & PropTrawl == "70%")

lh1 <- ggplot(data = subset(sum_data,GearType=="Fished")) +
  aes(x = year, y = avg_lenatmat, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_lenatmat - (sd_lenatmat/sqrt(20)), ymax = avg_lenatmat + (sd_lenatmat/sqrt(20)), 
                  fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average length at maturity (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_grid(PropTrawl~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
ggsave(here("R/Rplots3/length_at_mat_time.pdf"), plot=lh1, width = 30, height = 20, units = "cm")
lh2 <- ggplot(data = subset(sum_data,GearType=="Fished")) +
  aes(x = year, y = avg_ageatmat, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_ageatmat - (sd_ageatmat/sqrt(20)), ymax = avg_ageatmat + (sd_ageatmat/sqrt(20)),
                  fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average age at maturity (y)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw()+
  facet_grid(PropTrawl~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
ggsave(here("R/Rplots3/age_at_mat_time.pdf"), plot=lh2, width = 30, height = 20, units = "cm")


lh3 <- ggplot(data = lh_dat2) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = avg_lenatmat, color = PropTrawl, fill = PropTrawl) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_lenatmat-sd_lenatmat, ymax = avg_lenatmat+sd_lenatmat),position=position_dodge()) +
  #scale_color_manual(values = temp_palette) +
  #scale_fill_manual(values = temp_palette) +
  labs(x = "Size at maximum selectivity (cm)", y = "Average length at maturity (cm)", color = "Proportion trawled", fill = "Proportion trawled") +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10)) 
ggsave(here("R/Rplots3/length_at_mat_lmax.pdf"), plot=lh3, width = 30, height = 20, units = "cm")

lh4 <- ggplot(data = lh_dat2) +
  aes(x = Lmax, y = avg_ageatmat, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymin = avg_ageatmat-sd_ageatmat, ymax = avg_ageatmat+sd_ageatmat),position=position_dodge()) +
  scale_color_manual(values = temp_palette) +
  #scale_fill_manual(values = temp_palette) +
  labs(x = "Size at maximum selectivity (cm)", y = "Average age at maturity (y)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
ggsave(here("R/Rplots3/age_at_mat_lmax.pdf"), plot=lh4, width = 30, height = 20, units = "cm")

lh5 <- ggplot(data = lh_dat2) +
  aes(x = PropTrawl, y = avg_lenatmat, color = factor(Lmax, levels=c("90","100","110","120","130")), 
      fill = factor(Lmax, levels=c("90","100","110","120","130"))) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_lenatmat-sd_lenatmat, ymax = avg_lenatmat+sd_lenatmat),position=position_dodge()) +
  #scale_color_manual(values = temp_palette) +
  #scale_fill_manual(values = temp_palette) +
  labs(x = "Proportion trawled", y = "Average length at maturity (cm)", color = "Proportion trawled", fill = "Proportion trawled") +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10)) 
ggsave(here("R/Rplots3/length_at_mat_proptrawl.pdf"), plot=lh5, width = 30, height = 20, units = "cm")

lh6 <- ggplot(data = lh_dat2) +
  aes(x = PropTrawl, y = avg_ageatmat, color = factor(Lmax, levels=c("90","100","110","120","130")), 
      fill = factor(Lmax, levels=c("90","100","110","120","130"))) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_ageatmat-sd_ageatmat, ymax = avg_ageatmat+sd_ageatmat),position=position_dodge()) +
  #scale_color_manual(values = temp_palette) +
  #scale_fill_manual(values = temp_palette) +
  labs(x = "Proportion trawled", y = "Average age at maturity (y)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
ggsave(here("R/Rplots3/age_at_mat_proptrawl.pdf"), plot=lh6, width = 30, height = 20, units = "cm")

#Raster plots-----
lh7 <- ggplot(data=lh_dat2) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_lenatmat)) +
  geom_tile(data=frames,aes(fill=avg_lenatmat),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average length at\nmaturity (cm)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "top")

lh8 <- ggplot(data=lh_dat2) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_ageatmat)) +
  geom_tile(data=frames,aes(fill=avg_ageatmat),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average age at\nmaturity (y)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 10), legend.position = "top")
  
lh_arr1 <- ggarrange(lh7,lh8, ncol=1, nrow=2, common.legend = F, legend = "right", align="v")
ggsave(here("R/Rplots3/maturation_raster.pdf"), plot=lh_arr1, width = 20, height = 10, units = "cm")


aov_ageatmat <- aov(data=lh_dat2, avg_ageatmat~ClimScen+PropTrawl+Lmax+
                PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_ageatmat)

sink(here("R/rplots3/anova_ageatmat.txt"))
capture.output(summary(aov_ageatmat, file= here("R/rplots3/anova_ageatmat.txt")))
sink()
aov_lenatmat <- aov(data=lh_dat2, avg_lenatmat~ClimScen+PropTrawl+Lmax+
                      PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_lenatmat)

sink(here("R/rplots3/anova_lenatmat.txt"))
capture.output(summary(aov_lenatmat, file= here("R/rplots3/anova_lenatmat.txt")))
sink()

rm(lh_dat2,frames, lh1,lh2,lh3,lh4,lh5,lh6,lh7,lh8,lh_arr1)
gc()


## Growth ----------------------------------------------------------------------------------------------
gro_dat = subset(age_data, year == 2000 & GearType=="Fished")
gro_dat2 = subset(age_data, year == 2000 & GearType=="Fished" & age == "10")
frames <- subset(gro_dat2, Lmax == "110" & PropTrawl == "70%")

gro3 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_SomWeight, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_SomWeight - sd_SomWeight, ymax = avg_SomWeight + sd_SomWeight, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_grid(PropTrawl~factor(Lmax, levels=c("90","100","110","120","130"))) +
  labs(x= "Age (y)", y = "Somatic weight (kg)") +
  theme_bw()
ggsave(here("R/Rplots3/growth_curve_weight.pdf"), plot=gro3, width = 30, height = 20, units = "cm")
gro4 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_length, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_length - sd_length, ymax = avg_length + sd_length, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_grid(factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))~factor(Lmax, levels=c("90","100","110","120","130"))) +
  labs(x= "Age (y)", y = "Length (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw()+
  theme(legend.position = "top")
ggsave(here("R/Rplots3/growth_curve_length.pdf"), plot=gro4, width = 20, height = 22, units = "cm")


gro5 <- ggplot(data=gro_dat2) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_length)) +
  geom_tile(data=frames,aes(fill=avg_length),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average length at\nage 10 (cm)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 9), legend.position = "right")

ggsave(here("R/Rplots3/growth_raster.pdf"), plot=gro5, width = 20, height = 5, units = "cm")

aov_gro <- aov(data=gro_dat2, avg_length~ClimScen+PropTrawl+Lmax+
                 PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_gro)

sink(here("R/rplots3/anova_lenat10.txt"))
capture.output(summary(aov_gro, file= here("R/rplots3/anova_lenat10.txt")))
sink()

rm(gro_dat,gro3,gro4)
gc()

## Mortality ----------------------------------------------------------------------------------------------
mort_data2 <- subset(age_data, year > 1890 & GearType=="Fished") %>%
  group_by(Fmax,ClimScen,Lmax,PropTrawl, age) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            Mpred = mean(avg_Mpred), Mpred_sd = combine_sd(avg_Mpred,sd_Mpred,avg_n),
            Mfor = mean(avg_Mfor), Mfor_sd = combine_sd(avg_Mfor,sd_Mfor,avg_n),
            Mresp = mean(avg_Mresp), Mresp_sd = combine_sd(avg_Mresp,sd_Mresp,avg_n),
            Mrepro = mean(avg_Mrepro), Mrepro_sd = combine_sd(avg_Mrepro,sd_Mrepro,avg_n),
            Fishmort = mean(avg_Fishmort), Fishmort_sd = combine_sd(avg_Fishmort,sd_Fishmort,avg_n),
            Z = mean(avg_Z), Z_sd = combine_sd(avg_Z, sd_Z, avg_n)
  )

mort1 <- ggplot(data = mort_data2) +
  aes(x = age, y = Z, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Z - Z_sd, ymax = Z + Z_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Total mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_grid(PropTrawl~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme_bw()
ggsave(here("R/Rplots3/mortality_total.pdf"), plot=mort1, width = 30, height = 20, units = "cm")

mort2 <- ggplot(data = mort_data2) +
  aes(x = age, y = Fishmort, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Fishmort - Fishmort_sd, ymax = Fishmort + Fishmort_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Fisheries mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_grid(PropTrawl~factor(Lmax, levels=c("90","100","110","120","130"))) +
  theme_bw()
ggsave(here("R/Rplots3/mortality_fisheries.pdf"), plot=mort2, width = 30, height = 20, units = "cm")

rm(mort_data2,mort1,mort2)
gc()


## Yield ----------------------------------------------------------------------------------------------
yield_data <- subset(age_data,GearType=="Fished" & year == 2000) %>%
  group_by(GearType, ClimScen, Fmax,PropTrawl,Lmax, year) %>%
  summarise(tot_yield = sum(avg_yield)/1000, caught = sum(avg_caught))
yield_data$avg_size <- (yield_data$tot_yield*1000) / yield_data$caught
frames <- subset(yield_data, Lmax == "110" & PropTrawl == "70%")

yie1 <- ggplot(data=yield_data) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=tot_yield)) +
  geom_tile(data=frames,aes(fill=tot_yield),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average yield (t)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 9), legend.position = "right")

yie2 <- ggplot(data=yield_data) +
  aes(x = factor(Lmax, levels=c("90","100","110","120","130")), y = factor(PropTrawl, levels = c("0%","30%","50%","70%","100%"))) +
  geom_raster(aes(fill=avg_size)) +
  geom_tile(data=frames,aes(fill=avg_size),color="black", alpha=0.2, size=0.5) +
  labs(x= "Target length (cm)", y= "Proportion trawled", fill = "Average weight\nof caught\nfish (kg)") +
  scale_fill_paletteer_c(raster_palette) +
  theme_bw() +
  facet_wrap(~ClimScen, nrow = 1) +
  theme(text = element_text(size = 9), legend.position = "right")

ggsave(here("R/Rplots3/yield_raster.pdf"), plot=yie1, width = 20, height = 5, units = "cm")

yie_arr1 <- ggarrange(yie1,yie2, ncol=1, nrow=2, common.legend = F, legend = "right", align = "v")
ggsave(here("R/Rplots3/yield_raster2.pdf"), plot=yie_arr1, width = 20, height = 10, units = "cm")

aov_yie <- aov(data=yield_data, tot_yield~ClimScen+PropTrawl+Lmax+
                PropTrawl*Lmax+PropTrawl*ClimScen+Lmax*ClimScen)
summary(aov_yie)

sink(here("R/rplots3/anova_yield.txt"))
capture.output(summary(aov_yie, file= here("R/rplots3/anova_yield.txt")))
sink()




## Combination plots ----------------------------------------------------------------------------------------------
comb1 <- (pop7+pop8)/yie1
comb1 <- ggarrange(pop7,pop8,yie1, ncol=1, nrow=3, common.legend = F, legend = "right")
ggsave(here("R/Rplots3/popdynamics+yield.pdf"), plot=comb1, width = 25, height = 25, units = "cm")
