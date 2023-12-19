# Summary --------------------------------------------------------------------------------------
# This is the script for generating plots, to use in my paper on how climate and fisheries interact to cause evolution
# Plots are organised into subheading for ease of navigating this script. 
# Data import takes 20-30 minutes, recommend saving workspace when exiting.
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

## Load data ------------------------------------------------------------------------------------
#Based on netCDF
#Import the data files
tic()
data_files <- list.files(path = here("output_paper2") , pattern = "*.nc") #create list of data files
setwd(here("output_paper2"))
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
run_names1 <- vector(mode="character",length = ncol(run_names_holder1)-1)
run_names2 <- vector(mode="character",length = ncol(run_names_holder2)-1)
run_names3 <- vector(mode="character",length = ncol(run_names_holder2)-1)

for (i in 1:(ncol(run_names_holder1)-1)){
    run_names1[i] <- run_names_holder1[1,i+1]
}
for (i in 1:(ncol(run_names_holder2)-1)){
  run_names2[i] <- run_names_holder2[1,i+1]
}
for (i in 1:(ncol(run_names_holder3)-1)){
  run_names3[i] <- run_names_holder3[1,i+1]
}
run_names3 <- as.double(run_names3) %>%
  round(digits = 2) %>%
  as.character()

PopResults$GearType <- as.character(PopResults$run)
PopResults$ClimScen <- as.character(PopResults$run)
PopResults$Fmax <- as.character(PopResults$run)

for (i in 1:length(run_names1)) {
  PopResults$GearType[PopResults$GearType == i] <- run_names1[i]
}
for (i in 1:length(run_names2)) {
  PopResults$ClimScen[PopResults$ClimScen == i] <- run_names2[i]
}
for (i in 1:length(run_names3)) {
  PopResults$Fmax[PopResults$Fmax == i] <- run_names3[i]
}

rm(run_names_holder1, run_names_holder2, run_names_holder3,
   i, run_names1, run_names2, run_names3)



## Summarize the data -----------------------------------------------------------------------------------------------------
PopResults$GearType[PopResults$GearType == "0"] <- "Unfished"
PopResults$GearType[PopResults$GearType == "4"] <- "Fished"
PopResults$Fmax[PopResults$GearType == "Unfished"] <- "0.0"
PopResults$ClimScen[PopResults$ClimScen == "0"] <- "No warming"
PopResults$ClimScen[PopResults$ClimScen == "1"] <- "SSP1-2.6"
PopResults$ClimScen[PopResults$ClimScen == "2"] <- "SSP2-4.5"
PopResults$ClimScen[PopResults$ClimScen == "3"] <- "SSP3-7.0"
  


age_data <- PopResults %>%
  group_by(GearType, ClimScen, Fmax, year, age) %>%
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
            avg_gongro = mean(gon_growth), sd_gongro = sd(gon_growth)  #Gonadic growth
            )
age_data$avg_appetite <- age_data$avg_appetite/1000
age_data$sd_app <- age_data$sd_app/1000
age_data$avg_gonall <- age_data$avg_gonall*100
age_data$sd_gonall <- age_data$sd_gonall*100
age_data <- subset(age_data, year > 0)


sum_data <- PopResults %>%
  group_by(GearType, ClimScen, Fmax, year) %>%
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


## Temperature ------------------------------------------------------------------------------------------------------
temp_data <- PopResults %>%
  group_by(ClimScen, year) %>%
  summarise(temperature = mean(temp),sd_temperature = sd(temp))
temp_data <- subset(temp_data, year>0)
            
temp1 <- ggplot(data = temp_data) +
  aes(x = year, y = temperature, color = ClimScen, fill = ClimScen) +
  #geom_point() +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = temperature - sd_temperature, ymax = temperature + sd_temperature, fill=ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Temperature (°C)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() + 
  theme(legend.position = "top",text = element_text(size = 10))
temp1

ggsave(here("R/Rplots2/temperature.pdf"), plot=temp1, width = 14.5, height = 9, units = "cm")
rm(temp_data, temp1) 

## Population dynamics ----------------------------------------------------------------------------------------------
popdynamics_data <- age_data %>%
  group_by(GearType, ClimScen, Fmax, year) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            n_mat = sum(avg_mature), mat_sd = combine_sd(avg_mature, sd_mature, avg_n),
            pop_biomass = sum(avg_SomWeight*avg_n)/1000, biomass_sd = sum(sd_SomWeight * avg_n)/1000
            )
popdynamics_data <- subset(popdynamics_data, year >0)
popdynamics_data$rel_n <- popdynamics_data$n/121590.00
popdynamics_data$rel_biomass <- popdynamics_data$pop_biomass/296.3876
popdynamics_data$rel_n_sd <- popdynamics_data$n_sd/121590.00
popdynamics_data$rel_biomass_sd <- popdynamics_data$biomass_sd/296.3876
popdynamics_data$rel_biomass_sd
popdynamics_data$Fmax <- as.character(popdynamics_data$Fmax)
popdynamics_data$Fmax[popdynamics_data$Fmax == "0.0"] <- "Fmax = 0.0"
popdynamics_data$Fmax[popdynamics_data$Fmax == "0.1"] <- "Fmax = 0.1"
popdynamics_data$Fmax[popdynamics_data$Fmax == "0.2"] <- "Fmax = 0.2"
popdynamics_data$Fmax[popdynamics_data$Fmax == "0.3"] <- "Fmax = 0.3"
popdynamics_data$n <- popdynamics_data$n/1000
popdynamics_data$n_sd <- popdynamics_data$n_sd/1000

pop1 <- ggplot(popdynamics_data) +
  aes(x = year, y = n, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=n - n_sd, ymax = n + n_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Total number of\nindividuals (thousands)", color = "IPCC scenario", fill= "IPCC scenario") +
  theme_bw() +
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
pop2 <- ggplot(popdynamics_data) +
  aes(x = year, y = pop_biomass, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin= pop_biomass - biomass_sd, ymax = pop_biomass + biomass_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Population biomass (t)", color = "IPCC scenario", fill= "IPCC scenario") +
  theme_bw() +
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size =10))


pop3 <- ggplot(popdynamics_data) +
  aes(x = year, y = rel_n, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=rel_n - rel_n_sd, ymax = rel_n + rel_n_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Relative number of individuals", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))
pop4 <- ggplot(popdynamics_data) +
  aes(x = year, y = rel_biomass, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin= rel_biomass - rel_biomass_sd, ymax = rel_biomass + rel_biomass_sd, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Relative population biomass", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(text = element_text(size = 10))

pop7 <- ggplot(data = subset(popdynamics_data, year == "2000")) +
  aes(x = Fmax, y = n, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = n-n_sd, ymax = n+n_sd),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(y = "Average number of individuals (thousands)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pop7
pop8 <- ggplot(data = subset(popdynamics_data, year == "2000")) +
  aes(x = Fmax, y = pop_biomass, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = pop_biomass-biomass_sd, ymax = pop_biomass+biomass_sd),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(y = "Average total biomass (t)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pop8



popdat1 <-subset(PopResults, year > 1890) %>%
  group_by(GearType, ClimScen, Fmax, age) %>%
  summarise(temperature = mean(temp),
            avg_n = mean(n), sd_n = sd(n), avg_mature = mean(n_mature), sd_mature = sd(n_mature),         #Population dynamics
            avg_recr = mean(n_recr), sd_recr = sd(n_recr), avg_died = mean(n_died), sd_died = sd(n_died), #Population dynamics
  )
popdat1$age <- as.factor(popdat1$age)
popdat1$Fmax[popdat1$Fmax == "0.0"] <- "Fmax = 0.0"
popdat1$Fmax[popdat1$Fmax == "0.1"] <- "Fmax = 0.1"
popdat1$Fmax[popdat1$Fmax == "0.2"] <- "Fmax = 0.2"
popdat1$Fmax[popdat1$Fmax == "0.3"] <- "Fmax = 0.3"

pop5 <- ggplot(popdat1) +
  aes(x=Fmax,y=avg_n, color=age, fill=age) +
  geom_col(position="fill") +
  theme_bw() +
  scale_color_manual(values = big_palette) +
  scale_fill_manual(values = big_palette) +
  labs(x = bquote("Maximum fisheries mortality "(y^-1)), y = "Proportion of individuals in age group", color = "Age", fill = "Age") +
  theme(axis.text.x = element_text(angle = 60, hjust=1), axis.title.x = element_blank()) +
  facet_wrap(~ClimScen, nrow = 1)
pop5
pop6 <- ggplot(popdat1) +
  aes(x=ClimScen,y=avg_n, color=age, fill=age) +
  geom_col(position="fill") +
  theme_bw() +
  scale_color_manual(values = big_palette) +
  scale_fill_manual(values = big_palette) +
  labs(x = "IPCC scenario", y = "Proportion of individuals in age group", color = "Age", fill = "Age") +
  theme(axis.text.x = element_text(angle = 60, hjust=1), axis.title.x = element_blank()) +
  facet_wrap(~Fmax, nrow = 1)
pop6


#pop_arr1 <- ggarrange(pop1,pop2, ncol=1, nrow=2, common.legend = TRUE, legend = "top", labels = "AUTO",
                      #label.x = 0.0, font.label = list(size = 16))
pop_arr2 <- ggarrange(pop3,pop4, ncol=1, nrow=2, common.legend = TRUE, legend = "top", labels = "AUTO",
                      label.x = 0.0, font.label = list(size = 12))
pop_arr2
#ggsave(here("R/Rplots2/popdynamics.pdf"), plot=pop_arr1, width = 45.4, height = 30, units = "cm")
ggsave(here("R/Rplots2/relative_popdynamics.pdf"), plot=pop_arr2, width = 19, height = 15, units = "cm")

pop_arr3 <- ggarrange(pop5,pop6, ncol=1, common.legend = TRUE, legend = "right",
                      labels = c("B1","B2"),
                      label.x = 0.0, font.label = list(size = 12))
ggsave(here("R/Rplots2/age_distribution.pdf"), plot=pop_arr3, width = 20, height = 26, units = "cm")


pop_arr4 <- ggarrange(pop7,pop8, ncol=2, nrow=1, common.legend = TRUE, legend = "top",labels = "AUTO",
                      label.x = 0.0, font.label = list(size = 12))
ggsave(here("R/Rplots2/popdynamics2.pdf"), plot=pop_arr4, width = 19, height = 10, units = "cm")



rm(pop1,pop2,pop3,pop4,pop_arr1,pop_arr2, popdynamics_data)

## Evolving traits plots ------------------------------------------------------------------------------------------------
evo_data <- sum_data
evo_data$Fmax <- as.character(evo_data$Fmax)
evo_data$Fmax[evo_data$Fmax == "0.0"] <- "Fmax = 0.0"
evo_data$Fmax[evo_data$Fmax == "0.1"] <- "Fmax = 0.1"
evo_data$Fmax[evo_data$Fmax == "0.2"] <- "Fmax = 0.2"
evo_data$Fmax[evo_data$Fmax == "0.3"] <- "Fmax = 0.3"

evo1 <- ggplot(evo_data) +
  aes(x = year, y = avg_appetite, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_appetite - sd_app, ymax = avg_appetite + sd_app, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = bquote("Appetite "(KJ~kg^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme(text = element_text(size = 10)) +
  theme_bw()
evo1
evo2 <- ggplot(evo_data) +
  aes(x = year, y = avg_intercept, color = ClimScen) +
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int, fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "PMRN Intercept (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme(text = element_text(size = 12)) +
  theme_bw()
evo2

evo3 <- ggplot(data = subset(evo_data, year == "2000")) +
  aes(x = Fmax, y = avg_appetite, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_appetite-sd_app, ymax = avg_appetite+sd_app),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(y = bquote("Average appetite "(KJ~kg^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#evo3
evo4 <- ggplot(data = subset(evo_data, year == "2000")) +
  aes(x = Fmax, y = avg_intercept, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_intercept-sd_int, ymax = avg_intercept+sd_int),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(y = "Average PMRN intercept (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#evo4

# evo_arr1 <- ggarrange(evo1,evo2, ncol=1, nrow=2, common.legend = TRUE, legend = "top")
# ggsave(here("R/Rplots2/evolving_traits.pdf"), plot=evo_arr1, width = 35, height = 20, units = "cm")

evo_arr2 <- ggarrange(evo3,evo4, ncol=2, nrow=1, common.legend = TRUE, legend = "top", labels = "AUTO",
                      label.x = 0.0, font.label = list(size = 12))
ggsave(here("R/Rplots2/evolving_traits2.pdf"), plot=evo_arr2, width = 19, height = 10, units = "cm")

rm(evo_data, evo1,evo2, evo_arr1)

## Life history plots --------------------------------------------------------------------------------------------
lh_dat <- subset(PopResults, year > 1970)
lh_dat <- lh_dat %>%
  group_by(run, GearType, ClimScen, Fmax, year) %>%
  summarise(avg_lenatmat = wtd.mean(lenatmat, n_mature), sd_lenatmat= sqrt(wtd.var(lenatmat, n_mature)), #Lenatmat
            avg_ageatmat = wtd.mean(ageatmat, n_mature), sd_ageatmat= sqrt(wtd.var(ageatmat, n_mature)))
lh_dat2 <- subset(sum_data, year==2000)
lh_dat2$Fmax <- as.character(lh_dat2$Fmax)
lh_dat2$Fmax[lh_dat2$Fmax == "0.0"] <- "Fmax = 0.0"
lh_dat2$Fmax[lh_dat2$Fmax == "0.1"] <- "Fmax = 0.1"
lh_dat2$Fmax[lh_dat2$Fmax == "0.2"] <- "Fmax = 0.2"
lh_dat2$Fmax[lh_dat2$Fmax == "0.3"] <- "Fmax = 0.3"

lh1 <- ggplot(data = sum_data) +
  aes(x = year, y = avg_lenatmat, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_lenatmat - (sd_lenatmat/sqrt(20)), ymax = avg_lenatmat + (sd_lenatmat/sqrt(20)), 
                  fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average length at maturity (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(text = element_text(size = 12)) 
lh2 <- ggplot(data = sum_data) +
  aes(x = year, y = avg_ageatmat, color = ClimScen) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_ageatmat - (sd_ageatmat/sqrt(20)), ymax = avg_ageatmat + (sd_ageatmat/sqrt(20)),
                  fill = ClimScen), alpha = 0.2) +
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Year", y = "Average age at maturity (y)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw()+
  facet_wrap(~Fmax, nrow = 1) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(text = element_text(size = 12))

lh_arr1 <- ggarrange(lh1,lh2, nrow=2,common.legend = TRUE, legend = "top")
ggsave(here("R/Rplots2/maturation_schedule.pdf"), plot=lh_arr1, width = 35.4, height = 30, units = "cm")
#ggsave(here("R/Rplots2/maturation_schedule_poster.pdf"), plot=lh_arr1, width = 75, height = 40, units = "cm")

lh3 <- ggplot(data = lh_dat2) +
  aes(x = Fmax, y = avg_lenatmat, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_lenatmat-sd_lenatmat, ymax = avg_lenatmat+sd_lenatmat),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = bquote("Maximum fisheries mortality "(y^-1)), y = "Average length at maturity (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10),axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
lh3
lh4 <- ggplot(data = lh_dat2) +
  aes(x = Fmax, y = avg_ageatmat, color = ClimScen, fill = ClimScen) +
  geom_point(position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg_ageatmat-sd_ageatmat, ymax = avg_ageatmat+sd_ageatmat),position=position_dodge()) +
  geom_line(aes(group=ClimScen),position=position_dodge(width = 0.9), size=0.5, linetype="dashed")+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = bquote("Maximum fisheries mortality "(y^-1)), y = "Average age at maturity (y)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10),axis.title.x = element_blank()) +
  theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
lh4

lh_arr2 <- ggarrange(lh4,lh3, nrow=1, common.legend = TRUE, legend = "top",labels = "AUTO",
                     label.x = 0.0, font.label = list(size = 12))
ggsave(here("R/Rplots2/maturation_schedule2.pdf"), plot=lh_arr2, width = 19, height = 9, units = "cm")

rm(lh1,lh2,lh_arr1, lh3,lh4,lh_arr2)

## Mortality plots ------------------------------------------------------------------------------------------------
mort_data1 <- subset(age_data, year == 20 | year == 2000) %>%
  group_by(year, Fmax, age) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            Mpred = mean(avg_Mpred), Mpred_sd = combine_sd(avg_Mpred,sd_Mpred,avg_n),
            Mfor = mean(avg_Mfor), Mfor_sd = combine_sd(avg_Mfor,sd_Mfor,avg_n),
            Mresp = mean(avg_Mresp), Mresp_sd = combine_sd(avg_Mresp,sd_Mresp,avg_n),
            Mrepro = mean(avg_Mrepro), Mrepro_sd = combine_sd(avg_Mrepro,sd_Mrepro,avg_n),
            Fishmort = mean(avg_Fishmort), Fishmort_sd = combine_sd(avg_Fishmort,sd_Fishmort,avg_n),
            Z = mean(avg_Z), Z_sd = combine_sd(avg_Z, sd_Z, avg_n)
            )
mort_data1$year <- as.factor(mort_data1$year)


mort1 <- ggplot(data = mort_data1) +
  aes(x = age, y = Z, color = Fmax, fill = Fmax, linetype = year, shape = year) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Z - Z_sd, ymax = Z + Z_sd, fill = Fmax), alpha = 0.2)+
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Age (y)", y = bquote("Total mortality "(y^-1))) +
  #facet_wrap(~year, nrow = 1) +
  theme_bw() +
  theme(text = element_text(size = 12))
mort1
ggsave(here("R/Rplots2/tot_mortality.pdf"), plot=mort1, width = 30, height = 20, units = "cm")


mort_data2 <- subset(age_data, year > 1890) %>%
  group_by(Fmax,ClimScen, age) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            Mpred = mean(avg_Mpred), Mpred_sd = combine_sd(avg_Mpred,sd_Mpred,avg_n),
            Mfor = mean(avg_Mfor), Mfor_sd = combine_sd(avg_Mfor,sd_Mfor,avg_n),
            Mresp = mean(avg_Mresp), Mresp_sd = combine_sd(avg_Mresp,sd_Mresp,avg_n),
            Mrepro = mean(avg_Mrepro), Mrepro_sd = combine_sd(avg_Mrepro,sd_Mrepro,avg_n),
            Fishmort = mean(avg_Fishmort), Fishmort_sd = combine_sd(avg_Fishmort,sd_Fishmort,avg_n),
            Z = mean(avg_Z), Z_sd = combine_sd(avg_Z, sd_Z, avg_n)
  )
mort_data2$Fmax <- as.character(mort_data2$Fmax)
mort_data2$Fmax[mort_data2$Fmax == "0.0"] <- "Fmax = 0.0"
mort_data2$Fmax[mort_data2$Fmax == "0.1"] <- "Fmax = 0.1"
mort_data2$Fmax[mort_data2$Fmax == "0.2"] <- "Fmax = 0.2"
mort_data2$Fmax[mort_data2$Fmax == "0.3"] <- "Fmax = 0.3"

mort2 <- ggplot(data = mort_data2) +
  aes(x = age, y = Z, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Z - Z_sd, ymax = Z + Z_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Total mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()
mort2
ggsave(here("R/Rplots2/tot_mortality2.pdf"), plot=mort2, width = 19, height = 10, units = "cm")

mort3 <- ggplot(data = mort_data2) +
  aes(x = age, y = Mpred, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Mpred - Mpred_sd, ymax = Mpred + Mpred_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Predation mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
mort4 <- ggplot(data = mort_data2) +
  aes(x = age, y = Mfor, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Mfor - Mfor_sd, ymax = Mfor + Mfor_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Foraging mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
mort5 <- ggplot(data = mort_data2) +
  aes(x = age, y = Mresp, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Mresp - Mresp_sd, ymax = Mresp + Mresp_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Respiration mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
mort6 <- ggplot(data = mort_data2) +
  aes(x = age, y = Mrepro, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Mrepro - Mrepro_sd, ymax = Mrepro + Mrepro_sd, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Reproductive mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
mort7 <- ggplot(data = mort_data2) +
  aes(x = age, y = Fishmort, color = ClimScen, fill = ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = Fishmort - Fishmort, ymax = Fishmort + Fishmort, fill = ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  labs(x = "Age (y)", y = bquote("Fisheries mortality "(y^-1)), color = "IPCC scenario", fill = "IPCC scenario") +
  facet_wrap(~Fmax, nrow = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())

mort_arr1 <- ggarrange(mort3,mort4,mort5,mort6,mort7,mort2, nrow=6,common.legend = TRUE,align = "v",
                       legend = "top", labels = c("C1","C2","C3","C4","C5","C6"),
                       label.x = 0.0, font.label = list(size = 12), hjust = -1.6, vjust = 0.1)
#mort_arr1 <- annotate_figure(mort_arr1, top = text_grob("Average of last 100 years runtime"))
ggsave(here("R/Rplots2/mort_breakdown.pdf"), plot=mort_arr1, width = 20, height = 30, units = "cm")

rm(mort_data1, mort_data2, mort1,mort2,mort3,mort4,mort5,mort6,mort_arr1)

## Growth patterns ------------------------------------------------------------------------------------------------
gro_dat = subset(age_data, year == 2000)
gro_dat$Fmax <- as.character(gro_dat$Fmax)
gro_dat$Fmax[gro_dat$Fmax == "0.0"] <- "Fmax = 0.0"
gro_dat$Fmax[gro_dat$Fmax == "0.1"] <- "Fmax = 0.1"
gro_dat$Fmax[gro_dat$Fmax == "0.2"] <- "Fmax = 0.2"
gro_dat$Fmax[gro_dat$Fmax == "0.3"] <- "Fmax = 0.3"

gro1 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_SomGro, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_SomGro - sd_SomGro, ymax = avg_SomGro + sd_SomGro, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_wrap(~Fmax, nrow = 1) +
  labs(x= "Age (y)", y = "Somatic growth (kg)") +
  theme_bw()
gro1
gro2 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_GonGro, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_GonGro - sd_GonGro, ymax = avg_GonGro + sd_GonGro, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_wrap(~Fmax, nrow = 1) +
  labs(x= "Age (y)", y = "Gonadal growth (kg)") +
  theme_bw()
gro2
gro3 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_SomWeight, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_SomWeight - sd_SomWeight, ymax = avg_SomWeight + sd_SomWeight, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_wrap(~Fmax, nrow = 1) +
  labs(x= "Age (y)", y = "Somatic weight (kg)") +
  theme_bw()
gro3
gro4 <- ggplot(data=gro_dat) +
  aes(x = age, y = avg_length, color = ClimScen, fill=ClimScen) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_length - sd_length, ymax = avg_length + sd_length, fill=ClimScen), alpha = 0.2)+
  scale_color_manual(values = temp_palette) +
  scale_fill_manual(values = temp_palette) +
  facet_wrap(~Fmax, nrow = 1) +
  labs(x= "Age (y)", y = "Length (cm)", color = "IPCC scenario", fill = "IPCC scenario") +
  theme_bw() +
  theme(legend.position = "top") 
gro4
ggsave(here("R/Rplots2/growth_curves2.pdf"), plot=gro4, width = 19, height = 10, units = "cm")

gro_arr1 <- ggarrange(gro3,gro4, nrow=2,common.legend = TRUE, legend = "top")
gro_arr1 <- annotate_figure(gro_arr1, top = text_grob("End of 2000 year runtime"))
ggsave(here("R/Rplots2/growth_curves.pdf"), plot=gro_arr1, width = 19, height = 20, units = "cm")
rm(gro1,gro2,gro3,gro4,gro_arr1, gro_dat)

## Combined plots ------------------------------------------------------------------------------------------------
comb1 <- ggarrange(evo2,lh_arr2, nrow=2,common.legend = TRUE, legend = "top", labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 13))
comb1
ggsave(here("R/Rplots2/maturation_plot.pdf"), plot=comb1, width = 19, height = 15, units = "cm")

comb2 <- ggarrange(evo1,gro1,gro4, nrow=3,common.legend = TRUE, legend = "top", labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 13))
comb2
ggsave(here("R/Rplots2/growth_plot.pdf"), plot=comb2, width = 19, height = 15, units = "cm")

comb3
comb4 <- ggarrange(evo1,evo2,pop1,pop2, nrow=4,common.legend = TRUE, legend = "top", align = "v",
                   labels = c("A1","A2","A3","A4"),
                   label.x = 0.0, font.label = list(size = 12))
comb4
ggsave(here("R/Rplots2/time_appendix.pdf"), plot=comb4, width = 19, height = 25, units = "cm")
