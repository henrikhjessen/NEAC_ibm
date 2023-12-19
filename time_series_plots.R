# Set up the workspace -------------------------------------------------------------------------
library(here)
library(tidyverse)
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
library(plotly)
library(htmlwidgets)
library(shiny)
library(esquisse)
library(gifski)
library(transformr)
library(gganimate)

# Load data ------------------------------------------------------------------------------------
#Based on netCDF
#Import the data files
data_files <- list.files(path = here("output") , pattern = "*.nc") #create list of data files
setwd(here("output"))
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

min_year <- min(PopResults$year[PopResults$year>0])
#min_year <- min(PopResults$year)
max_year <- max(PopResults$year)
pop_start <- subset(PopResults, year==min_year)
pop_start$year <- as.factor(pop_start$year)
pop_end <- subset(PopResults, year==max_year)
pop_end$year <- as.factor(pop_end$year)
before_n_after <- rbind(pop_start,pop_end)

rm(min_year,max_year,pop_start,pop_end)


#Import data from ICES -------------------------------------------------------------------------------------------
#Import ICES data
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


#Subset the data by runs [Optional] -------------------------------------------------------------------------------------------
select_runs <- c(18:23)
PopResults <- subset(PopResults, (run %in% select_runs))
before_n_after <- subset(before_n_after, (run %in% select_runs))
rm(select_runs)


#Name the runs [Optional] -------------------------------------------------------------------------------------------
run_names_holder <- subset(parameters, parameter_names=="InvestDecay") #Just edit the parameter_name here to the testing 
run_names <- vector(mode="character",length = ncol(run_names_holder)-1)  #variable to know the values for each run
run_names_holder <- as.matrix(run_names_holder[1,])
for (i in 1:(ncol(run_names_holder)-1)){
  run_names[i] <- run_names_holder[1,i+1]
}
rm(run_names_holder,i)


#Cumulative fecundity plots-----------------------------------------------------------------------------
CumFec <- subset(PopResults, select = c(run, year, age, n, CumGonWeight, CumGonWeight_SD))
CumFec <- CumFec %>%
  group_by(run, year) %>%
  summarise(cum_gon_wt = weighted.mean(CumGonWeight, n))

p <- ggplot(CumFec) +
 aes(x = year, y = cum_gon_wt, colour = run) +
 geom_point(shape = "circle", size = 1.5, colour = "#112446") +
 geom_smooth(span = 0.3) +
 xlab("Year") +
 ylab("Cumulative gonad weight (kg)") +
 theme_bw()
p
ggsave(here("R/Rplots/cumulative_fecundity.pdf"), plot=p, width = 20, height = 20, units = "cm")
rm(CumFec, p)

#Population dynamics and Biomass plots-----------------------------------------------------------------------------
PopDynamics <- subset(PopResults, select = c(run, year, n_recr, n_survived,n_mature, n_died, n))
PopDynamics <- PopDynamics %>%
  group_by(run, year) %>%
  summarise(n_total=sum(n), n_recr=mean(n_recr), n_died=sum(n_died), n_survived=sum(n_survived), n_mature=sum(n_mature))
PopDynamics_long <- melt(PopDynamics, id.vars = c("run","year"))

popbiomass <- subset(PopResults, select = c(run, year, pop_biomass))
popbiomass <- popbiomass %>%
  group_by(run, year) %>%
  summarise(biomass=(sum(pop_biomass))/1000)


p1 <- ggplot(popbiomass) +
 aes(x = year, y = biomass, colour = run) +
 geom_point(shape = "circle", size = 1.5) +
 geom_smooth(span = 0.48) +
 xlab("Year") +
 ylab("Population Biomass (t)") +
 scale_color_hue(direction = 1) +
 theme_bw()+
  theme(axis.title.x = element_blank())
p1

p2 <- ggplot(PopDynamics) +
  aes(x = year, y = n_total, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.2) +
  xlab("Year") +
  ylab("Total number of fish") +
  scale_color_hue(direction = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
p2

p3 <- ggplot(PopDynamics) +
 aes(x = year, y = n_mature, colour = run) +
 geom_point(shape = "circle", size = 1.5) +
 geom_smooth(span = 0.2) +
 xlab("Year") +
 ylab("Number of mature fish") +
 scale_color_hue(direction = 1) +
 theme_bw()+
  theme(axis.title.x = element_blank())
p3

p4 <- ggplot(PopDynamics) +
  aes(x = year, y = n_recr, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.2) +
  xlab("Year") +
  ylab("Number of recruits") +
  scale_color_hue(direction = 1) +
  theme_bw()+
  theme(axis.title.x = element_blank())
p4
p0 <- ggarrange(p1,p2,p3,p4, ncol=1, nrow=4, common.legend = TRUE, legend = "right")
p0 <- annotate_figure(p0, bottom = "Years")
p0
ggsave(here("R/Rplots/popdynamics.pdf"), plot=p0, width = 40, height = 50, units = "cm")

#saveWidget(ggplotly(p1),here("R/Rplots/Interactives/biomass_over_time.html"))
#saveWidget(ggplotly(p1),here("R/Rplots/Interactives/ntot_over_time.html"))
#saveWidget(ggplotly(p3),here("R/Rplots/Interactives/nmature_over_time.html"))
#saveWidget(ggplotly(p4),here("R/Rplots/Interactives/nrecr_over_time.html"))

rm(PopDynamics, PopDynamics_long, popbiomass, p1,p2,p3,p4,p0)

#Mortality plot ------------------------------------------------------------------------------------

p1 <- ggplot(before_n_after) +
 aes(x = age, y = M_predation, colour = run) +
 geom_point(shape = "circle", size = 1.5) +
 geom_smooth(span = 0.75) +
 scale_color_hue(direction = 1) +
 labs(x = "Age (y)", y = "Predation mortality (y^-1)") +
 theme_bw() +
 facet_wrap(vars(year)) +
 theme(axis.title.x = element_blank())
p2 <- ggplot(before_n_after) +
  aes(x = age, y = M_foraging, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Foraging mortality (y^-1)") +
  theme_bw() +
  facet_wrap(vars(year)) +
  theme(axis.title.x = element_blank())
p3 <- ggplot(before_n_after) +
  aes(x = age, y = M_reproduction, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Reproductive mortality (y^-1)") +
  theme_bw() +
  facet_wrap(vars(year)) +
  theme(axis.title.x = element_blank())
p4 <- ggplot(before_n_after) +
  aes(x = age, y = M_respiration, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Respiration mortality (y^-1)") +
  theme_bw() +
  facet_wrap(vars(year)) +
  theme(axis.title.x = element_blank())
p5 <- ggplot(subset(PopResults, year == 20)) +
  aes(x = age, y = FishMort, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Fisheries mortality (y^-1)") +
  theme_bw() +
  #facet_wrap(vars(year)) +
  theme(axis.title.x = element_blank())
p5
p6 <- ggplot(before_n_after) +
  aes(x = age, y = Z, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Total mortality, Z (y^-1)") +
  theme_bw() +
  facet_wrap(vars(year)) +
  theme(axis.title.x = element_blank())

p <- ggarrange(p1,p2,p3,p4,p5,p6, ncol=1, nrow=6, common.legend = TRUE, legend = "right")
p <- annotate_figure(p, bottom = "Age (y)")
p
ggsave(here("R/Rplots/mort_overview.pdf"), plot=p, width = 30, height = 70, units = "cm")
rm(p,p1,p2,p3,p4,p5,p6)

p1 <- ggplot(subset(PopResults, run == 3)) +
  aes(x = age, y = FishMort) +
  geom_point(shape = "circle", size = 1.5) +
  geom_line(aes(group=year)) +
  #geom_ribbon(aes(ymin=FishMort - FishMort_SD, ymax = FishMort + FishMort_SD, fill = run), alpha = 0.2) +
  labs(title = 'Year: {frame_time}', x = "Age", y = "Total mortality (y^-1)") +
  theme(text = element_text(size = 12)) +
  theme_bw()  +
  transition_time(year) +
  ease_aes('linear')
anim_save(here("R/Rplots/anim_mortality.gif"), p1, duration = 20, fps = 30, width = 25.4, height = 14.3, units = "cm", 
          res= 200, renderer = gifski_renderer(loop = FALSE))


#Growth plot -------------------------------------------------------------------------------------------
p <-  ggplot(before_n_after) +
  aes(x = age, y = som_growth, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  labs(x = "Age (y)", y = "Somatic growth (kg y^-1)") +
  theme_bw() +
  facet_wrap(vars(year))
ggsave(here("R/Rplots/growth_by_age.pdf"), plot=p,width = 30, height = 20, units = "cm")
rm(p)

#Foraging plot -------------------------------------------------------------------------------------------
unfish <- ggplot(data=subset(PopResults, FishMort == 0)) +
  geom_histogram(aes(x=f_int), binwidth = 0.5) +
  theme_bw()
fish <- ggplot(data=subset(PopResults, FishMort != 0)) +
  geom_histogram(aes(x=f_int), binwidth = 0.5) +
  theme_bw()
unfish/fish
rm(unfish, fish)


#MISC. PLOTS------------------------------------------------------------------------------------
#Growth by biomass
ggplot(data=PopResults, aes(x=pop_biomass, y=tot_growth)) +
  geom_point() +
  theme_classic()
ggsave(here("R/Rplots/growth_by_biomass.pdf"), plot=last_plot())

#Mortality by biomass
tot_natmort <- (PopResults$M_predation + PopResults$M_foraging + PopResults$M_reproduction + PopResults$M_respiration)
PopResults <- add_column(PopResults, tot_natmort, .after = "FishMort")
ggplot(data=PopResults, aes(x=pop_biomass, y=tot_natmort)) +
  geom_point() +
  theme_classic()
ggsave(here("R/Rplots/mortality_by_biomass.pdf"), plot=last_plot())

#Temperature over the years
p1 <-ggplot(PopResults) +
 aes(x = year, y = temp, colour = run) +
 geom_point(shape = "circle", size = 1.5) +
 geom_smooth(span = 0.5) +
 scale_color_hue(direction = 1) +
 xlab("Year") +
 ylab("Temperature (C)") +
 theme_bw()
ggsave(here("R/Rplots/temperature_by_year.pdf"), plot=p1,width = 30, height = 20, units = "cm")
rm(p1)

#Survival plot
surv_df <- cbind(PopResults$year, PopResults$age, PopResults$M_foraging + PopResults$M_predation + 
                   PopResults$M_reproduction + PopResults$M_respiration + PopResults$FishMort ) %>%
  as.data.frame()
surv_df[,4] <- exp(-surv_df[,3])
colnames(surv_df)[1:4] <- c("year","age","z","S")

ggplot(data=surv_df) +
  geom_point(aes(x=age, y=S))+
  theme_bw()
ggsave(here("R/Rplots/survival_by_age.pdf"), plot=last_plot())
rm(surv_df)


#GSI plot 
GSI <- as.data.frame(cbind(PopResults$year, PopResults$age, PopResults$GonWeight/PopResults$SomWeight))
colnames(GSI)[1:3] <- c("year","age","gsi")

ggplot(data=GSI) +
  geom_point(aes(x=age, y=gsi))+
  theme_bw()
ggsave(here("R/Rplots/gsi_by_age.pdf"), plot=last_plot())


#Percentage mature by age
ggplot(data=PopResults) +
  geom_point(aes(x=age, y=pct_mature, color = year)) +
  theme_bw()
ggsave(here("R/Rplots/pct_mature_by_age.pdf"), plot=last_plot())
#Percentage mature by length
ggplot(data=PopResults) +
  geom_point(aes(x=length, y=pct_mature)) +
  theme_bw()
ggsave(here("R/Rplots/pct_mature_by_length.pdf"), plot=last_plot())


#Recruitment by year plot
ggplot(data=PopResults) +
  geom_point(aes(x = year, y= n_recr)) +
  theme_bw()


#Skipping spawners by age
ggplot(data=PopResults) +
  #geom_point(aes(x=age, y=n_spawnskip)) +
  geom_point(aes(x=age, y=pct_skippers)) +
  theme_classic()
ggsave(here("R/Rplots/pct_spawnskipper_by_age.pdf"), plot=last_plot())




#Heritability plots--------------------------------------------------------------------------------------------------

aamat <- ggplot(heritabilities) +
 aes(x = year, y = ageatmat, colour = run) +
 geom_point(shape = "circle", size = 1.5) +
  ggtitle("Age at maturation") +
  xlab("Year") +
  ylab("Heritability") +
 geom_smooth(span = 0.46) +
 scale_color_hue(direction = 1) +
 theme_bw()
lamat <- ggplot(heritabilities) +
  aes(x = year, y = lengthatmat, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  ggtitle("Length at maturation") +
  xlab("Year") +
  ylab("Heritability") +
  geom_smooth(span = 0.46) +
  scale_color_hue(direction = 1) +
  theme_bw()
p <- ggarrange(aamat, lamat, ncol=1, nrow=2, common.legend = TRUE, legend = "right")
p
ggsave(here("R/Rplots/heritabilities.pdf"), plot=p, width = 30, height = 20, units = "cm")
rm(aamat, lamat, p)

#Evolving Plots + Maturity-----------------------------------------------------------------------------------------------------------------------
mat_dat <- subset(PopResults, select = c(run, year,age,n,n_mature, ageatmat, lenatmat, gon_allocation, intercept, appetite))
mat_dat_avg <- mat_dat %>%
  group_by(run, year) %>%
  summarise(tot_n=sum(n), tot_n_mature=sum(n_mature), ageatmat=weighted.mean(ageatmat,n_mature), lenatmat=weighted.mean(lenatmat,n_mature),
            intercept=weighted.mean(intercept,n_mature), allocation=(weighted.mean(gon_allocation,n_mature))*100, 
            appetite=weighted.mean(appetite,n))
p1 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year,y=intercept, color=run), size = 1.4) +
  #scale_color_hue("Run",labels = c("Unfished","Fished")) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("PMRN Intercept (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank())

p2 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year,y=allocation, color=run), size = 1.4) +
  #scale_color_hue("Run",labels = c("Unfished", "Fished")) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("Gonadic allocation (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
  
p3 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year,y=ageatmat, color=run), size = 1.4) +
  #scale_color_hue("Run",labels = c("Unfished","Fished")) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("Age at maturation (y)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p4 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year,y=lenatmat, color=run), size = 1.4) +
  #scale_color_hue("Run",labels = c("Unfished","Fished")) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("Length at maturation (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank())

p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend = "right")
p <- annotate_figure(p, bottom = "Years")
p
ggsave(here("R/Rplots/maturation_evolution.pdf"), plot=p, width = 30, height = 20, units = "cm")

p5 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year, y= appetite, color=run), size =1.4) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ylab("Appetite") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p6 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year, y= allocation, color=run), size =1.4) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("Gonadic allocation (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p7 <- ggplot(data=mat_dat_avg) +
  geom_line(aes(x=year, y= intercept, color=run), size =1.4) +
  #scale_color_brewer("Run",labels = c("Unfished","Fished"),palette="Paired") +
  #scale_color_hue("Run",labels = run_names) +
  ylab("PMRN Intercept (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p0 <- ggarrange(p5,p6,p7, ncol=1, common.legend = TRUE, legend = "right")
p0 <- annotate_figure(p0, bottom="Year")
p0
ggsave(here("R/Rplots/evolving_traits.pdf"), plot=p0, width = 30, height = 20, units = "cm")

rm(mat_dat, mat_dat_avg, p1,p2,p3,p4,p5,p6,p7,p,p0)

#Evolving Growth Plots------------------------------------------------------------------------------------------------------------------------------
grow_dat <- subset(PopResults, year>20, select = c(run, year,age,n,SomWeight, GonWeight, length, som_growth, gon_growth, gon_allocation, appetite))
grow_dat$age <- as.factor(grow_dat$age)


p1 <- grow_dat %>%
  filter(age %in% c("3", "6", "10")) %>%
  ggplot() +
  aes(x = year, y = SomWeight, fill = age, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.8) +
  scale_fill_hue(direction = 1) +
  scale_color_hue(direction = 1) +
  labs(x = "Year", y = "Somatic weight (kg)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p2 <-grow_dat %>%
  filter(age %in% c("3", "6", "10")) %>%
  ggplot() +
  aes(x = year, y = GonWeight, fill = age, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.8) +
  scale_fill_hue(direction = 1) +
  scale_color_hue(direction = 1) +
  labs(x = "Year", y = "Gonad weight (kg)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p3 <- grow_dat %>%
  filter(age %in% c("3", "6", "10")) %>%
  ggplot() +
  aes(x = year, y = length, fill = age, colour = run) +
  geom_point(shape = "circle", size = 1.5) +
  geom_smooth(span = 0.8) +
  scale_fill_hue(direction = 1) +
  scale_color_hue(direction = 1) +
  labs(x = "Year", y = "Length (cm)") +
  theme_bw() +
  theme(axis.title.x = element_blank())
p <- ggarrange(p1,p2,p3, ncol=1, nrow=3, common.legend = TRUE, legend = "right")
p <- annotate_figure(p, bottom = "Year")
p
ggsave(here("R/Rplots/growth_evolution.pdf"), plot=p, width = 30, height = 20, units = "cm")
rm(p,p1,p2,p3,grow_dat, grow_dat_avg)


#Growth trajectory plots--------------------------------------------------------------------------------------------------

before_n_after$ices_w <- 0
before_n_after$ices_l <- 0

for (i in 1:nrow(before_n_after)) {
  if (before_n_after[i,"age"] <= 14) {
    before_n_after[i,"ices_w"] = ices_dat[i,"weight"]
    before_n_after[i,"ices_l"] = ices_dat[i,"length"]
  } else {
    before_n_after[i,"ices_w"] = NA
    before_n_after[i,"ices_l"] = NA    
  }
}

p1 <- ggplot(before_n_after) +
  aes(x = age, y = SomWeight, colour = run, linetype = year) +
  ylab("Somatic weight (kg)") +
  geom_line(size = 1.2) +
  #geom_smooth(span = 0.6) +
  scale_color_hue(direction = 1) +
  geom_point(aes(x=age, y=ices_w),shape=15, color="black",size = 1.5)+
  theme_bw() +
  theme(axis.title.x = element_blank())
p1
p2 <- ggplot(before_n_after) +
  aes(x = age, y = length, colour = run, linetype = year) +
  ylab("Length (cm)") +
  geom_line(size = 1.2) +
  #geom_smooth(span = 0.6) +
  scale_color_hue(direction = 1) +
  geom_point(aes(x=age, y=ices_l),shape=15, color="black",size = 1.5)+
  theme_bw() +
  theme(axis.title.x = element_blank())
p2

p3 <- ggarrange(p1,p2, ncol=2, nrow=1, common.legend = TRUE, legend = "top")
p3 <- annotate_figure(p3, bottom = "Age (y)")
p3
ggsave(here("R/Rplots/growth_trajectory.pdf"), plot=p3, width = 40, height = 20, units = "cm")


#saveWidget(ggplotly(p1),here("R/Rplots/Interactives/growth_trajectories.html"))

rm(p1,p2,p3,i)



#Mortality plot round 2-------------------------------------------------------------------------------------------
min_year <- min(PopResults$year)
max_year <- max(PopResults$year)
mort_start <- subset(PopResults, year==min_year)
mort_start$year <- as.factor(mort_start$year)
mort_end <- subset(PopResults, year==max_year)
mort_end$year <- as.factor(mort_end$year)
mort <- rbind(mort_start,mort_end)
mort <- mort[,c("run","year","age","M_predation","M_foraging","M_reproduction","M_respiration","FishMort")]
mort <- melt(mort, id.vars = c("run","year","age"))






p1 <- ggplot() +
  geom_line(data = subset(mort,year==min_year),size=1.4, aes(x=age,y=value, color=run, linetype=variable)) +
  ylab("Mortality (y^-1)") +
  labs(color = "Run", variable="Mortality") +
  scale_color_brewer("Run",labels = run_names) +
  theme_bw() +
  theme(axis.title.x = element_blank())
p2 <- ggplot() +
  geom_line(data = subset(mort,year==max_year),size=1.4, aes(x=age,y=value, color=run, linetype=variable)) +
  ylab("Mortality (y^-1)") +
  labs(color = "Run", variable="Mortality") +
  scale_color_brewer("Run",labels = run_names) +
  theme_bw() +
  theme(axis.title.x = element_blank())
  
p <- ggarrange(p1,p2, ncol=2, nrow=1, common.legend = TRUE, legend = "right")
p <- annotate_figure(p, bottom = "Age (y)")
p
ggsave(here("R/Rplots/mortality_lines.pdf"), plot=p, width = 40, height = 20, units = "cm")

saveWidget(ggplotly(p1),here("R/Rplots/Interactives/Mort_before.html"))
saveWidget(ggplotly(p2),here("R/Rplots/Interactives/Mort_after.html"))
rm(p,p1,p2,mort,mort_end,mort_start,min_year,max_year)
