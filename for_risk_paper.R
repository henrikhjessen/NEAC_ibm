# Summary --------------------------------------------------------------------------------------
# This is the script for generating plots, to use in my paper on how foraging risk impacts life-history evolution.
# Plots are organised into subheading for ease of navigating this script. Animations are found at the end of the script.
#Data import takes 20-30 minutes, recommend saving workspace when exiting.
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
library(tictoc)
library(Hmisc)
library(gifski)
library(transformr)

# Define new functions ------------------------------------------------------------------------------------
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

# Load data ------------------------------------------------------------------------------------
#Based on netCDF
#Import the data files
tic()
data_files <- list.files(path = here("output - paper1") , pattern = "*.nc") #create list of data files
setwd(here("output - paper1"))
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
colnames(PopResults)[1:45] <-c("run","year","age","temp","FoodEnv","n_recr","n_survived","n_died","SomWeight","SomWeight_SD",
                               "GonWeight","GonWeight_SD","length","length_SD","M_predation","Mpred_SD","M_foraging",
                               "Mfor_SD","M_reproduction","Mrepr_SD","M_respiration","Mresp_SD","FishMort","FishMort_SD",
                               "f_int","fint_SD","n_mature","n_spawnskip","pop_biomass","som_growth", "somgro_SD",
                               "gon_growth", "gongro_SD","appetite","appetite_SD","som_allocation","allocation_SD",
                               "intercept","intercept_SD","ageatmat","ageatmat_SD","lenatmat","lenatmat_SD","CumGonWeight",
                               "CumGonWeight_SD") #name columns
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
Z <- PopResults$M_predation +PopResults$M_foraging+PopResults$M_reproduction+PopResults$M_respiration+PopResults$FishMort
PopResults <- add_column(PopResults, Z, .after = "FishMort_SD")
rm(data_holder, array_years,array_runs, data_holder_length, i, n, pct_mature, pct_skippers, 
   tot_growth, gon_allocation, temp_PopResults, Z) #cleanup

PopResults <- PopResults %>% drop_na(year)
heritabilities <- subset(heritabilities, select = -c(run))
rm(ncdf_dataframe, data_files)
toc()
#Name the runs  --------------------------------------------------------------------------------------------------------
run_names_holder <- subset(parameters, parameter_names=="for_risk") #Just edit the parameter_name here to the testing 
run_names <- vector(mode="character",length = ncol(run_names_holder)-1)  #variable to know the values for each run
run_names_holder <- as.matrix(run_names_holder[1,])
for (i in 1:(ncol(run_names_holder)-1)){
  run_names[i] <- run_names_holder[1,i+1]
}
run_names <- as.double(run_names)
run_names <- round(run_names, digits = 2)
run_names <- as.character(run_names)
PopResults$run <- as.character(PopResults$run)
for (i in 1:length(run_names)) {
  PopResults$run[PopResults$run == i] <- run_names[i]
}
rm(run_names, i, run_names_holder)
gc()

#Import ICES data  --------------------------------------------------------------------------------------------------------
ices_l_data <- read_delim(here("R/ices_dat/ices_length_data.csv"), 
                          ";", escape_double = FALSE, trim_ws = TRUE)
ices_length <- colMeans(subset(ices_l_data, year > "2004"),na.rm = T) 
ices_length <- ices_length[2:15]
ices_l_age <- c(1:14)
ices_length <- as.data.frame(cbind(ices_l_age,ices_length))
colnames(ices_length)[1:2] <- c("age","length")

ices_w_data <- read_delim(here("R/ices_dat/ices_weight_data.csv"), 
                          ";", escape_double = FALSE, trim_ws = TRUE)
ices_weight <- colMeans(subset(ices_w_data, Year_age > "2004"))
ices_weight <- ices_weight[2:13]  
#ices_age <- c(3:14)
#ices_weight <- as.data.frame(cbind(ices_age,ices_weight))
ices_weight <- c(NA,NA,ices_weight)
ices_dat <- cbind(ices_length,ices_weight)
colnames(ices_dat)[1:3] <- c("age","length","weight")

rm(ices_w_data,ices_l_data,ices_l_age,ices_length,ices_weight)

#Select runs  --------------------------------------------------------------------------------------------------------
select_runs <- c(1,1.4,1.8,2.2,2.6)
PopResults <- subset(PopResults, (run %in% select_runs))
rm(select_runs)
#Summarize the data -----------------------------------------------------------------------------------------------------
## First step - Summarize for the repetitions, maintaining age groups, reps are not weighted differently
age_data <- PopResults %>%
  group_by(run, year, age) %>%
  summarise(avg_temp = mean(temp), sd_temp = sd(temp), avg_FoodEnv = mean(FoodEnv), sd_FoodEnv = sd(FoodEnv),
            avg_fint = mean(f_int), sd_fint = sd(f_int), avg_appetite = mean(appetite), sd_appetite=sd(appetite),
            avg_intercept = mean(intercept), sd_intercept=sd(intercept), avg_gonall = mean(gon_allocation),
            sd_gonall = sd(gon_allocation), avg_ageatmat = mean(ageatmat), sd_ageatmat = sd(ageatmat),
            avg_lenatmat = mean(lenatmat), sd_lenatmat = sd(lenatmat), avg_somweight = mean(SomWeight),
            sd_somweight = sd(SomWeight), avg_gonweight = mean(GonWeight), sd_gonweight = sd(GonWeight),
            avg_n = mean(n), sd_n = sd(n), avg_mature = mean(n_mature), sd_mature = sd(n_mature), avg_recr = mean(n_recr),
            sd_recr = sd(n_recr), avg_died = mean(n_died), sd_died = sd(n_died), avg_biomass = mean(pop_biomass),
            sd_biomass = sd(pop_biomass), avg_length = mean(length), sd_length = sd(length), avg_Mpred = mean(M_predation),
            sd_Mpred = sd(M_predation), avg_Mfor = mean(M_foraging), sd_Mfor = sd(M_foraging), avg_Mrepro = mean(M_reproduction),
            sd_Mrepro = sd(M_reproduction), avg_Mresp = mean(M_respiration), sd_Mresp = sd(M_respiration), avg_Z = mean(Z),
            sd_Z = sd(Z), avg_cumgonweight = CumGonWeight, sd_cumgonweight = CumGonWeight_SD)

age_data$year <- as.integer(age_data$year)
#Evolving traits ---------------------------------------------------------------------------------------------------------
## One-step process, pooling Rep-Age groups
evo_data <- PopResults %>%
  group_by(run, year) %>%
  summarise(avg_appetite = wtd.mean(appetite, n), sd_app = sqrt(wtd.var(appetite, n)),se_app = sqrt(wtd.var(appetite, n))/sqrt(20),
  avg_intercept = wtd.mean(intercept, n), sd_int = sqrt(wtd.var(intercept, n)),se_int = sqrt(wtd.var(intercept, n))/sqrt(20),
  avg_gonall = wtd.mean(gon_allocation, n), sd_gonall = sqrt(wtd.var(gon_allocation, n)),se_gonall = sqrt(wtd.var(gon_allocation, n))/sqrt(20))

evo_data$avg_appetite <- evo_data$avg_appetite/1000
evo_data$sd_app <- evo_data$sd_app/1000
evo_data$avg_gonall <- evo_data$avg_gonall*100
evo_data$sd_gonall <- evo_data$sd_gonall*100

evo1 <- ggplot(evo_data) +
  aes(x = year, y = avg_appetite, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_appetite - sd_app, ymax = avg_appetite + sd_app, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Appetite (KJ kg^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
evo1
evo2 <- ggplot(evo_data) +
  aes(x = year, y = avg_intercept, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "PMRN Intercept (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
evo2
evo3 <- ggplot(evo_data) +
  aes(x = year, y = avg_gonall, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_gonall - sd_gonall, ymax = avg_gonall + sd_gonall, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Initial gonadic allocation (%)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
evo3
evo4 <- ggplot(subset(evo_data, year == 2000)) +
  aes(x = run, y = avg_intercept, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "PMRN Intercept (cm)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
evo4
evo5 <- ggplot(subset(evo_data, year == 2000)) +
  aes(x = run, y = avg_gonall, color=run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=avg_gonall - sd_gonall, ymax = avg_gonall + sd_gonall)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Initial gonadic allocation (%)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position = "none")
evo5
evo_data$run <- as.character(evo_data$run)

evo0 <- ggarrange(evo1,evo2,evo3, ncol=1, common.legend = TRUE, legend = "top", labels = "AUTO",
                label.x = 0.0, font.label = list(size = 16))
ggsave(here("paper1_runs/figures/evolving_traits.pdf"), plot=evo0, width = 19, height = 30, units = "cm")
ggsave(here("paper1_runs/figures/evolving_appetite.svg"), plot=evo1, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/evolving_intercept.svg"), plot=evo2, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/evolving_allocation.svg"), plot=evo3, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/intercept_foraging.png"), plot=evo4, width = 25.4, height = 14.3, units = "cm")
rm(evo0, evo1,evo2,evo3, evo_data)

#Maturation schedule ---------------------------------------------------------------------------------------------------------
## One-step process, pooling Rep-Age groups
mat_data <- PopResults %>%
  group_by(run,year) %>%
  summarise(avg_lenatmat = wtd.mean(lenatmat, n_mature), sd_lenatmat = sqrt(wtd.var(lenatmat, n_mature)), 
            se_lenatmat = sqrt(wtd.var(lenatmat, n_mature))/sqrt(20),
            avg_ageatmat = wtd.mean(ageatmat, n_mature), sd_ageatmat = sqrt(wtd.var(ageatmat, n_mature)), 
            se_ageatmat = sqrt(wtd.var(ageatmat, n_mature))/sqrt(20))
## Two-step process, unweighted rep averaging, now weighted age averaging to follow
# mat_data2 <- age_data %>%
#   group_by(run, year) %>%
#   summarise(new_lenatmat = wtd.mean(avg_lenatmat, avg_mature), lenatmat_sd = combine_sd(avg_lenatmat,sd_lenatmat,avg_mature),
#             se_lenatmat = combine_sd(avg_lenatmat,sd_lenatmat,avg_mature)/sqrt(20),
#             new_ageatmat = wtd.mean(avg_ageatmat, avg_mature), ageatmat_sd = combine_sd(avg_ageatmat,sd_ageatmat,avg_mature),
#             se_ageatmat = combine_sd(avg_ageatmat,sd_ageatmat,avg_mature)/sqrt(20))
# colnames(mat_data2) <- colnames(mat_data)

mat1 <- ggplot(mat_data) +
  aes(x = year, y = avg_lenatmat, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_lenatmat - se_lenatmat, ymax = avg_lenatmat + se_lenatmat, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Length at maturation (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mat1
mat2 <- ggplot(mat_data) +
  aes(x = year, y = avg_ageatmat, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_ageatmat - se_ageatmat, ymax = avg_ageatmat + se_ageatmat, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Age at maturation (y)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mat2
mat0 <- ggarrange(mat1,mat2, ncol=2,nrow=1, common.legend = TRUE, legend = "top", labels = "AUTO",
                label.x = 0.0, font.label = list(size = 16))
mat3 <- ggarrange(mat1,mat2, ncol=2, common.legend = TRUE, legend = "bottom", labels = "AUTO",
                label.x = 0.0, font.label = list(size = 16))


mat4 <- ggplot(subset(mat_data, year == 2000)) +
  aes(x = run, y = avg_lenatmat, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=avg_lenatmat - se_lenatmat, ymax = avg_lenatmat + se_lenatmat)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Length at maturation (cm)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
mat4
mat5 <- ggplot(subset(mat_data, year == 2000)) +
  aes(x = run, y = avg_ageatmat, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=avg_ageatmat - se_ageatmat, ymax = avg_ageatmat + se_ageatmat)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Age at maturation (y)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
mat5
mat6 <- ggarrange(mat4,mat5,ncol=2, labels = "AUTO",label.x = 0.0, font.label = list(size = 16))

ggsave(here("paper1_runs/figures/maturation.pdf"), plot=mat0, width = 19, height = 9.6, units = "cm")
ggsave(here("paper1_runs/figures/length_at_maturation.svg"), plot=mat1, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/age_at_maturation.pdf"), plot=mat2, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/maturation_wide.svg"), plot=mat3, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/for_length_at_maturation.pdf"), plot=mat4, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/for_age_at_maturation.pdf"), plot=mat5, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/foraging_maturation.png"), plot=mat6, width = 25.4, height = 14.3, units = "cm")
rm(mat0,mat1,mat2,mat_data)


#Mortality ---------------------------------------------------------------------------------------------------------
mort_data <- age_data %>%
  group_by(run, year) %>%
  summarise(
    Z = wtd.mean(avg_Z, avg_n), Z_sd = combine_sd(avg_Z,sd_Z,avg_n), Z_se = combine_sd(avg_Z,sd_Z,avg_n)/sqrt(20),
    M_for = wtd.mean(avg_Mfor,avg_n), Mfor_sd = combine_sd(avg_Mfor,sd_Mfor,avg_n),
    Mfor_se = combine_sd(avg_Mfor,sd_Mfor,avg_n)/sqrt(20),
    M_pred = wtd.mean(avg_Mpred,avg_n), Mpred_sd = combine_sd(avg_Mpred,sd_Mpred,avg_n),
    Mpred_se = combine_sd(avg_Mpred,sd_Mpred,avg_n)/sqrt(20),
    M_resp = wtd.mean(avg_Mresp,avg_n), Mresp_sd = combine_sd(avg_Mresp,sd_Mresp,avg_n),
    Mresp_se = combine_sd(avg_Mresp,sd_Mresp,avg_n)/sqrt(20),
    M_repro = wtd.mean(avg_Mrepro,avg_n), Mrepro_sd = combine_sd(avg_Mrepro,sd_Mrepro,avg_n),
    Mrepro_se = combine_sd(avg_Mrepro,sd_Mrepro,avg_n)/sqrt(20)
  )

mort1 <- ggplot(subset(age_data, year == 2000)) +
  aes(x = age, y = avg_Mfor, color = run, fill = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_Mfor - (sd_Mfor/sqrt(20)), ymax = avg_Mfor + (sd_Mfor/sqrt(20)), fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Age", y = "Foraging mortality (y^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mort1
mort_after <- subset(age_data, year==2000)
mort2 <- ggplot(subset(mort_after, run==1 | run==1.8 | run==2.6)) +
  aes(x = age, y = avg_Z, color = run, fill = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_Z - (sd_Z/sqrt(20)), ymax = avg_Z + (sd_Z/sqrt(20)), fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette[c(1,3,5)]) +
  scale_fill_manual(values = my_palette[c(1,3,5)]) +
  labs(x = "Age", y = bquote("Total mortality"~(y^-1)), color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mort2

mort0 <- ggarrange(mort1,mort2, ncol=2, common.legend = TRUE, legend = "bottom", labels = "AUTO",
                label.x = 0.0, font.label = list(size = 16))
ggsave(here("paper1_runs/figures/foraging_mortality.svg"), plot=mort1, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/total_mortality.svg"), plot=mort2, width = 25.4, height = 14.3, units = "cm")
ggsave(here("paper1_runs/figures/mortality.pdf"), plot=mort0, width = 25.4, height = 14.3, units = "cm")


mort3 <- ggplot(subset(age_data, year == 1)) +
  aes(x = age, y = avg_Mfor, color = run, fill = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_Mfor - (sd_Mfor/sqrt(20)), ymax = avg_Mfor + (sd_Mfor/sqrt(20)), fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  ylim(0,0.9) +
  labs(x = "Age", y = "Foraging mortality (y^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mort3

mort4 <- ggplot(subset(age_data, year == 1)) +
  aes(x = age, y = avg_Z, color = run, fill = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_Z - (sd_Z/sqrt(20)), ymax = avg_Z + (sd_Z/sqrt(20)), fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  ylim(0,0.9) +
  labs(x = "Age", y = "Total mortality (y^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
mort4

mort5 <- ggarrange(mort4,mort2 + rremove("ylab")+rremove("y.text")+rremove("y.ticks"), ncol=2, common.legend = TRUE, legend = "bottom")
ggsave(here("paper1_runs/figures/mortality_beforenafter.svg"), plot=mort5, width = 25.4, height = 14.3, units = "cm")
rm(mort0,mort1,mort2,mort3,mort4,mort5,mort_data)

# Growth rates ---------------------------------------------------------------------------------------------------------

grow_dat <- subset(age_data, year != 1)
grow_dat1 <- subset(age_data, year == 1)
grow_dat1$age <- grow_dat1$age+1
grow_dat <- rbind(subset(grow_dat1,age!=0),grow_dat)

grow1 <- ggplot(subset(age_data, year == 1)) +
  aes(x = age, y = avg_somweight) +
  geom_point(aes(color = run, fill = run)) +
  geom_point(data=ices_dat, aes(x=age, y=weight)) +
  geom_line() +
  geom_ribbon(aes(ymin = avg_somweight - sd_somweight, ymax = avg_somweight + sd_somweight, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme(text = element_text(size = 12)) +
  labs(x = "Age", y = "Somatic Weight (kg)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw()
grow1
grow2 <- ggplot(subset(age_data, year == 2000)) +
  aes(x = age, y = avg_somweight, color = run, fill = run) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_somweight - sd_somweight, ymax = avg_somweight + sd_somweight, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme(text = element_text(size = 12)) +
  labs(x = "Age", y = "Somatic Weight (kg)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw()
grow2
grow3 <- ggplot(subset(age_data_after, run == 1 | run==1.8 | run==2.6)) +
  aes(x = age, y = avg_length, color = run, fill = run) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_length - sd_length, ymax = avg_length + sd_length, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette[c(1,3,5)]) +
  scale_fill_manual(values = my_palette[c(1,3,5)]) +
  theme(text = element_text(size = 12)) +
  labs(x = "Age", y = "Length (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw()
grow3
grow4 <- ggplot(subset(grow_dat, age == 7)) +
  aes(x = year, y = avg_length, color = run, fill = run) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_length - sd_length, ymax = avg_length + sd_length, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme(text = element_text(size = 12)) +
  labs(x = "Year", y = "Length (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw()
grow4

ggsave(here("paper1_runs/figures/length_at_age.svg"), plot=grow3, width = 25.4, height = 14.3, units = "cm")
rm(grow1,grow2,grow3,grow4,grow_dat, grow_dat1)

# Population plots --------------------------------------------------------------------------------------------------
population_dat <- subset(age_data, run == 1 | run == 1.4 | run == 1.8 | run == 2.2 | run == 2.6)%>%
  group_by(run, year) %>%
  summarise(n = sum(avg_n), n_sd = combine_sd(avg_n,sd_n,avg_n),
            n_mat = sum(avg_mature), mat_sd = combine_sd(avg_mature, sd_mature, avg_n),
            pop_biomass = sum(avg_somweight*avg_n)/1000, biomass_sd = sum(sd_somweight * avg_n)/1000,
            CUMGON = mean(avg_cumgonweight), CUMGON_sd = combine_sd(avg_cumgonweight,sd_cumgonweight, avg_n))
dummy_data_after <- subset(age_data, year ==2000 & age == 10) %>% 
  group_by(run, year) %>%
  summarise(CUMGON = mean(avg_cumgonweight), CUMGON_sd = combine_sd(avg_cumgonweight,sd_cumgonweight, avg_n))
dummy_data_after <- subset(dummy_data_after, run == 1 | run == 1.4 | run == 1.8 | run == 2.2 | run == 2.6)

pop1 <- ggplot(subset(population_dat, year == 2000)) +
  aes(x = run, y = n, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=n - n_sd, ymax = n + n_sd)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Number of individuals") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
pop1
pop2 <- ggplot(subset(population_dat, year == 2000)) +
  aes(x = run, y = n_mat, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=n_mat - mat_sd, ymax = n_mat + mat_sd)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Number of mature individuals") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
pop2
pop3 <- ggplot(subset(population_dat, year == 2000)) +
  aes(x = run, y = pop_biomass, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=pop_biomass - biomass_sd, ymax = pop_biomass + biomass_sd)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Population biomass (t)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
pop3
pop4 <- ggplot(subset(dummy_data_after)) +
  aes(x = run, y = CUMGON, color = run) +
  geom_point(shape = "circle", size = 2.5) +
  #geom_line() +
  geom_errorbar(aes(ymin=CUMGON - CUMGON_sd, ymax = CUMGON + CUMGON_sd)) +
  scale_color_manual(values = my_palette) +
  #scale_fill_manual(values = my_palette) +
  labs(x = "Foraging risk", y = "Average cumulative gonad weight \n of 10-year old fish (kg)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
pop4
pop5 <- ggplot(population_dat) +
  aes(x = year, y = pop_biomass, color = run, fill=run) +
  geom_point(shape = "circle", size = 2.5) +
  geom_line() +
  geom_ribbon(aes(ymin=pop_biomass - biomass_sd, ymax = pop_biomass + biomass_sd), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Population biomass (t)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw()
pop5

# Combination plots --------------------------------------------------------------------------------------------------
comb1 <- ggarrange(mort2,grow3, ncol=2,nrow=1, common.legend = TRUE, legend = "top", labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 16))
comb1
ggsave(here("paper1_runs/figures/combined_plot1.pdf"), plot=comb1, width = 19, height = 9.5, units = "cm")
comb2 <- ggarrange(evo1,grow3,mort2, ncol=1,nrow=3, common.legend = TRUE, legend = "top", labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 16))
comb2
ggsave(here("paper1_runs/figures/combined_plot2.pdf"), plot=comb2, width = 9.5, height = 27.5, units = "cm")
comb3 <- ggarrange(evo2, evo3, mat1, mat2, ncol=2,nrow=2, common.legend = TRUE, legend = "top", labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 16))
comb3
ggsave(here("paper1_runs/figures/combined_plot3.pdf"), plot=comb3, width = 19, height = 19, units = "cm")
comb4 <- ggarrange(evo4,evo5, mat4, mat5, nrow=2, ncol=2, labels = "AUTO",
                   label.x = 0.0, font.label = list(size = 16))
comb4 <- annotate_figure(comb4, bottom = "Foraging risk exponent")
ggsave(here("paper1_runs/figures/combined_plot4.pdf"), plot=comb4, width = 19, height = 19, units = "cm")
comb5 <- ggarrange(pop1,pop2,pop3,pop4, nrow=2,ncol=2, legend = "none") %>% 
  annotate_figure(bottom = "Foraging risk exponent")
comb6 <- ggarrange(pop5, comb5, ncol=1, common.legend = TRUE, heights = c(1,2))
ggsave(here("paper1_runs/figures/combined_plot6.pdf"), plot=comb6, width = 19, height = 27, units = "cm")
# Animations ---------------------------------------------------------------------------------------------------------
 ## Mortality
p1 <- ggplot(age_data) +
  aes(x = age, y = avg_Z, color = run, fill = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_Z - (sd_Z/sqrt(20)), ymax = avg_Z + (sd_Z/sqrt(20)), fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(title = 'Year: {frame_time}', x = "Age", y = "Total mortality (y^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw() + 
  transition_time(year) +
  ease_aes('linear')

# p1 <- ggplot(age_data) +
#   aes(x = age, y = avg_Z, color = run, fill = run) +
#   geom_point(shape = "circle", size = 1.5) +
#   geom_line() +
#   geom_ribbon(aes(ymin=avg_Z - (sd_Z), ymax = avg_Z + (sd_Z), fill = run), alpha = 0.2) +
#   scale_color_manual(values = my_palette) +
#   scale_fill_manual(values = my_palette) +
#   labs(title = 'Year: {frame_time}', x = "Age", y = "Total mortality (y^-1)", color = "Foraging risk", fill = "Foraging risk") +
#   theme(text = element_text(size = 12)) +
#   theme_bw() + 
#   transition_time(year) +
#   ease_aes('linear')

anim_save(here("paper1_runs/figures/anim_mortality.gif"), p1, duration = 20, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))

 ## Evolving traits
evo_data <- PopResults %>%
  group_by(run, year) %>%
  summarise(avg_appetite = wtd.mean(appetite, n), sd_app = sqrt(wtd.var(appetite, n)),se_app = sqrt(wtd.var(appetite, n))/sqrt(20),
            avg_intercept = wtd.mean(intercept, n), sd_int = sqrt(wtd.var(intercept, n)),se_int = sqrt(wtd.var(intercept, n))/sqrt(20),
            avg_gonall = wtd.mean(gon_allocation, n), sd_gonall = sqrt(wtd.var(gon_allocation, n)),se_gonall = sqrt(wtd.var(gon_allocation, n))/sqrt(20))
evo_data$avg_appetite <- evo_data$avg_appetite/1000
evo_data$sd_app <- evo_data$sd_app/1000
evo_data$avg_gonall <- evo_data$avg_gonall*100
evo_data$sd_gonall <- evo_data$sd_gonall*100

p2 <- ggplot(evo_data) +
  aes(x = year, y = avg_appetite, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_appetite - sd_app, ymax = avg_appetite + sd_app, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Appetite (KJ kg^-1)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  transition_reveal(along = year)

p3 <- ggplot(evo_data) +
  aes(x = year, y = avg_intercept, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "PMRN Intercept (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  transition_reveal(along = year)
p4 <- ggplot(evo_data) +
  aes(x = year, y = avg_gonall, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_gonall - sd_gonall, ymax = avg_gonall + sd_gonall, fill = run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Year", y = "Gonad allocation (%)", color = "Foraging risk", fill = "Foraging risk") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  transition_reveal(along = year)


anim_save(here("paper1_runs/figures/anim_appetite.gif"), p2, duration = 5, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))
anim_save(here("paper1_runs/figures/anim_intercept.gif"), p3, duration = 5, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))
anim_save(here("paper1_runs/figures/anim_gonall.gif"), p4, duration = 5, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))

 ## Growth rate
p5 <- ggplot(age_data) +
  aes(x = age, y = avg_somweight, color = run, fill = run) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_somweight - sd_somweight, ymax = avg_somweight + sd_somweight, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme(text = element_text(size = 12)) +
  labs(title = "Year: {frame_time}", x = "Age", y = "Somatic Weight (kg)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw() +
  transition_time(year)
p6 <- ggplot(age_data) +
  aes(x = age, y = avg_length, color = run, fill = run) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = avg_length - sd_length, ymax = avg_length + sd_length, fill= run), alpha = 0.2) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme(text = element_text(size = 12)) +
  labs(title = "Year: {frame_time}", x = "Age", y = "Length (cm)", color = "Foraging risk", fill = "Foraging risk") +
  theme_bw() +
  transition_time(year)

anim_save(here("paper1_runs/figures/anim_growth.gif"), p5, duration = 20, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))
anim_save(here("paper1_runs/figures/anim_growth_length.gif"), p6, duration = 20, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))

rm(p1,p2,p3,p4,p5,p6,evo_data)
