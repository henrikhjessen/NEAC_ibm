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


my_palette <- c("#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7")

# Load data ------------------------------------------------------------------------------------
data_files <- list.files(path = here("saved_pops") , pattern = "*.csv") #create list of data files
setwd(here("saved_pops"))
for (i in 1:length(data_files)) {
  data_holder <- read_csv(data_files[i],col_names = FALSE)
  run_holder <- rep(i, length(data_holder))
  data_holder <- cbind(run_holder, data_holder)
  if (i == 1) {
    dataframe <- data_holder
  } else {
    dataframe <- abind(dataframe, data_holder, along = 1)
  }
  rm(data_holder, run_holder)
}
dataframe <- as.data.frame(dataframe)
colnames(dataframe)[1:21] <-c("source_file","status","age","length","som_wt","stored_energy","maturity","gon_wt","cum_gon_wt",
                              "pmrn_intercept","pmrn_slope","pmrn_width","init_invest","appetite","midparent_ageatmat",
                              "phen_intercept","phen_slope","phen_initinvest","phen_appetite","lenatmat","midparent_lenatmat")
for (i in 1:length(data_files)) {
  dataframe["source_file"][dataframe["source_file"] == i] <- data_files[i]
}
init_pops <- subset(dataframe, status != 0)
init_pops <- subset(init_pops, age != 0)
rm(data_files, i, dataframe)
setwd(here())

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
ices_weight <- c(NA,NA,ices_weight)
ices_dat <- cbind(ices_length,ices_weight)
colnames(ices_dat)[1:3] <- c("age","length","weight")

ices_m_data <- read_delim(here("R/ices_dat/ices_maturity_data.csv"), 
                           ";", escape_double = FALSE, trim_ws = TRUE)
ices_mat_46 <- subset(ices_m_data, Year_age == "1946")
ices_mat_21 <- colMeans(subset(ices_m_data, Year_age > "2014"))
ices_mat_46 <- as.double(ices_mat_46[2:13])
ices_mat_21 <- as.double(ices_mat_21[2:13])
ices_mat_46 <- c(NA,NA,ices_mat_46)
ices_mat_21 <- c(NA,NA,ices_mat_21)

ices_dat <- cbind(ices_dat, ices_mat_46,ices_mat_21)

rm(ices_w_data,ices_l_data,ices_l_age,ices_length,ices_weight, ices_m_data, ices_mat)
gc()

# Generate plots ------------------------------------------------------------------------------------  
#Sizes
size_dat <- subset(init_pops, select = c(source_file, age, length, som_wt, maturity))
size_dat$maturity <- ifelse(size_dat$maturity>0,1,0)
size_dat <- size_dat %>%
  group_by(source_file, age) %>%
  summarise(mean_length=mean(length), sd_length=sd(length),
            mean_wt=mean(som_wt), sd_wt=sd(som_wt),
            sum_mat = sum(maturity), n_tot =n())
size_dat$prop_mat <- (size_dat$sum_mat / size_dat$n_tot) * 100

size_dat$ices_w <- 0
size_dat$ices_l <- 0
#size_dat$ices_m <- 0
for (i in 1:nrow(size_dat)) {
  if (size_dat[i,"age"] <= 14) {
    size_dat[i,"ices_w"] = ices_dat[i,"weight"]
    size_dat[i,"ices_l"] = ices_dat[i,"length"]
  } else {
    size_dat[i,"ices_w"] = NA
    size_dat[i,"ices_l"] = NA
  }
}

p1 <- ggplot(data = size_dat, aes(x = age, y =mean_wt)) +
  ylab("Somatic weight (kg)") +
  geom_line(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(color=source_file)) +
  geom_line(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(color=source_file)) +
  geom_point(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(color=source_file)) +
  geom_point(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(color=source_file)) +
  geom_ribbon(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(ymin=mean_wt-sd_wt, ymax=mean_wt+sd_wt, fill=source_file), alpha=0.25) +
  geom_ribbon(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(ymin=mean_wt-sd_wt, ymax=mean_wt+sd_wt, fill=source_file), alpha=0.25) +
  geom_point(aes(y=ices_w), shape=15) +
  scale_color_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  scale_fill_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p1
p2 <- ggplot(data = size_dat, aes(x = age, y =mean_length)) +
  ylab("Length (cm)") +
  geom_line(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(color=source_file)) +
  geom_line(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(color=source_file)) +
  geom_point(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(color=source_file)) +
  geom_point(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(color=source_file)) +
  geom_ribbon(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length, fill=source_file), alpha=0.25) +
  geom_ribbon(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length, fill=source_file), alpha=0.25) +
  scale_color_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  scale_fill_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  geom_point(aes(y=ices_l), shape=15) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p2
p3 <- ggplot(data = size_dat, aes(x = age, y =prop_mat)) +
  ylab("Proportion mature (%)") +
  geom_line(data = subset(size_dat, source_file == "20230621_fished_finalpop.csv"), aes(color=source_file), size=1.2) +
  geom_line(data = subset(size_dat, source_file == "20230621_unfished_finalpop.csv"), aes(color=source_file), size=1.2) +
  geom_point(data = ices_dat, aes(y=ices_mat_21*100, x = age), shape=15) +
  geom_point(data = ices_dat, aes(y=ices_mat_46*100, x = age), shape=0) +
  scale_color_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  scale_fill_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p3
p <- ggarrange(p2,p3, ncol=2, labels = "AUTO", common.legend = TRUE, label.x = 0.0, font.label = list(size = 16))
p <- annotate_figure(p, bottom = "Age (y)")
p
#ggsave(here("saved_pops/growth_trajectory.pdf"), plot=p, width = 40, height = 20, units = "cm")
ggsave(here("R/ices_calibration.pdf"), plot=p, width = 19, height = 10, units = "cm")
rm(p, p1, p2, size_dat, i)


#Inherited values
inherited_vals <- subset(init_pops, select = c(source_file, appetite, pmrn_intercept, init_invest))
inherited_vals$init_invest <- inherited_vals$init_invest*100
mean_vals <- inherited_vals %>%
  group_by(source_file) %>%
  summarise(mean_app=mean(appetite), sd_app=sd(appetite),
            mean_intercept=mean(pmrn_intercept), sd_intercept=sd(pmrn_intercept),
            mean_invest=mean(init_invest), sd_invest=sd(init_invest))
p1 <- ggplot(data=inherited_vals, aes(fill=source_file)) +
  ylab("Appetite (J kg^-1)") +
  geom_violin(aes(x=source_file, y=appetite)) +
  geom_text(data=mean_vals, aes(x=source_file, y = mean_app, label=mean_app), parse=T) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p1
p2 <- ggplot(data=inherited_vals, aes(fill=source_file)) +
  ylab("PMRN Intercept (cm)") +
  geom_violin(aes(x=source_file, y=pmrn_intercept)) +
  geom_text(data=mean_vals, aes(x=source_file, y = mean_intercept, label=mean_intercept), parse=T) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p2
p3 <- ggplot(data=inherited_vals, aes(fill=source_file)) +
  ylab("Initial somatic investment (%)") +
  geom_violin(aes(x=source_file, y=init_invest)) +
  geom_text(data=mean_vals, aes(x=source_file, y = mean_invest, label=mean_invest), parse=T) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p3
p <- ggarrange(p1,p2,p3, ncol=3, nrow=1, common.legend = TRUE, legend = "top")
p
ggsave(here("saved_pops/inherited_vals.pdf"), plot=p, width = 40, height = 20, units = "cm")
rm(p, p1, p2, p3, mean_vals, inherited_vals)
gc()



















# Time series
data_files <- list.files(path = here("saved_pops") , pattern = "*.nc") #create list of data files
setwd(here("saved_pops"))
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

#Evo-plots
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
  aes(x = year, y = avg_appetite, color = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_appetite - sd_app, ymax = avg_appetite + sd_app, color = run, fill = run), alpha = 0.2) +
  #scale_color_discrete(name=NULL) +
  #scale_fill_discrete(name=NULL) +
  labs(x = "Year", y = "Appetite (KJ kg^-1)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
evo1
evo2 <- ggplot(evo_data) +
  aes(x = year, y = avg_intercept, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line() +
  geom_ribbon(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int,colour = run, fill = run), alpha = 0.2) +
  #scale_color_discrete(name=NULL,labels = c("Fished","Unfished")) +
  #scale_fill_discrete(name=NULL,labels = c("Fished","Unfished")) +
  labs(x = "Year", y = "PMRN Intercept (cm)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
evo2
evo3 <- ggplot(evo_data) +
  aes(x = year, y = avg_gonall, colour = run) +
  geom_point(shape = "circle", size = 1.5,colour = run) +
  geom_line(colour = run) +
  geom_ribbon(aes(ymin=avg_gonall - sd_gonall, ymax = avg_gonall + sd_gonall,colour = run, fill = run), alpha = 0.2) +
  #scale_color_discrete(name=NULL,labels = c("Fished","Unfished")) +
  #scale_fill_discrete(name=NULL,labels = c("Fished","Unfished")) +
  labs(x = "Year", y = "Initial gonadic allocation (%)") +
  theme(text = element_text(size = 12)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none")
evo3
# evo4 <- ggplot(subset(evo_data, year == 2000)) +
#   aes(x = run, y = avg_intercept, color = run) +
#   geom_point(shape = "circle", size = 2.5) +
#   #geom_line() +
#   geom_errorbar(aes(ymin=avg_intercept - sd_int, ymax = avg_intercept + sd_int)) +
#   scale_color_discrete(name=NULL,labels = c("Fished","Unfished")) +
#   #scale_fill_discrete(name=NULL,labels = c("Fished","Unfished")) +
#   labs(x = "Foraging risk", y = "PMRN Intercept (cm)") +
#   theme(text = element_text(size = 12)) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
# evo4
# evo5 <- ggplot(subset(evo_data, year == 2000)) +
#   aes(x = run, y = avg_gonall, color=run) +
#   geom_point(shape = "circle", size = 2.5) +
#   #geom_line() +
#   geom_errorbar(aes(ymin=avg_gonall - sd_gonall, ymax = avg_gonall + sd_gonall)) +
#   scale_color_discrete(name=NULL,labels = c("Fished","Unfished")) +
#   #scale_fill_discrete(name=NULL,labels = c("Fished","Unfished")) +
#   labs(x = "Foraging risk", y = "Initial gonadic allocation (%)") +
#   theme(text = element_text(size = 12)) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position = "none")
# evo5
evo_data$run <- as.character(evo_data$run)


combine <- ggarrange(evo1,evo2, ncol=1, nrow=2, common.legend = T)
combine <- annotate_figure(combine, bottom = "Year (y)")
combine
ggsave(here("R/rplots2/init_pop_equilibrium.pdf"), plot=combine, width = 19, height = 20, units = "cm")





## Extra plot for ICES promat calibration
propmat_dat <- subset(PopResults, year >= 4800)
propmat_dat <- propmat_dat %>%
  group_by(run,age) %>%
  summarise(propmat = mean(pct_mature), sd_propmat=sd(pct_mature))

propmat_dat$ices_m46 <- 0
propmat_dat$ices_m21 <- 0
for (i in 1:nrow(propmat_dat)) {
  if (propmat_dat[i,"age"] <= 14) {
    propmat_dat[i,"ices_m46"] = ices_dat[i,"ices_mat_46"] *100
    propmat_dat[i,"ices_m21"] = ices_dat[i,"ices_mat_21"] *100
  } else {
    propmat_dat[i,"ices_m46"] = NA
    propmat_dat[i,"ices_m21"] = NA
  }
}


pm1 <- ggplot(data = propmat_dat, aes(x = age, y =propmat)) +
  ylab("Proportion mature (%)") +
  geom_point(aes(color=run)) +
  geom_line(aes(color=run)) +
  geom_ribbon(aes(ymin=propmat - sd_propmat, ymax = propmat + sd_propmat, fill = run), alpha = 0.2) +
  geom_point(aes(y=ices_m46), shape=0) +
  geom_point(aes(y=ices_m21), shape=15) +
  scale_color_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  scale_fill_manual(values=my_palette, name=NULL,labels = c("Fished","Unfished")) +
  theme_bw() +
  theme(axis.title.x = element_blank())
pm1
p <- ggarrange(p1,pm1, ncol=2, nrow=1, common.legend = TRUE, legend = "top", labels = "AUTO",
               label.x = 0.0, font.label = list(size = 16))
p <- annotate_figure(p, bottom = "Age (y)")
p
#ggsave(here("saved_pops/growth_trajectory.pdf"), plot=p, width = 40, height = 20, units = "cm")
ggsave(here("R/rplots/ices_calibration.pdf"), plot=p, width = 19, height = 10, units = "cm")



