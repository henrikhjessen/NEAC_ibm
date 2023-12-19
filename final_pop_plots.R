#Set up the workspace
library(here)
library(tidyverse)
library(readr)
library(dplyr)
setwd(here("R/RPlots"))
data <- read.csv(here("output/final_Pop.csv"))
colnames(data)[1:19] <-c("status", "age", "length", "SomWeight", "stored_energy", "maturity", "GonWeight", "cumfec", "gen_intercept",
                      "gen_slope", "pmrn_width", "gen_allocation", "gen_appetite", "midparental_aam",
                      "phen_intercept", "phen_slope", "phen_allocation", "phen_appetite", "length_at_maturation")
########################################################################################
#Age at maturity frequency distribution
########################################################################################
n_mature <- count(subset(data, maturity != "0"), maturity)
n_mature <- sum(n_mature$n)
ageatmat <- ggplot(data=subset(data, maturity != "0"), aes(x=maturity))
    ageatmat + geom_bar() +
    geom_text(stat="count", aes(label=..count..), vjust=-1) +
    xlab("Age at maturity") +
    ylab("Frequency") +
    scale_x_continuous(breaks = seq(0,12,1)) + 
    annotate("text", label = n_mature, x= 1, y= 20000) +
    theme_classic(base_size = 13)
    ggsave(here("R/Rplots/finalpop_ageatmat.pdf"), plot=last_plot())
    
#Import ICES length data
ices_l_data <- read_delim(here("R/ices_dat/ices_length_data.csv"), 
                 ";", escape_double = FALSE, trim_ws = TRUE)
ices_l_data <- colMeans(subset(ices_l_data, year > "2004"),na.rm = T) 
ices_l_data <- ices_l_data[2:15]
ices_l_age <- c(1:14)
ices_l_data <- as.data.frame(cbind(ices_l_age,ices_l_data))
colnames(ices_l_data)[1:2] <- c("age","length")

########################################################################################    
#Length at age
########################################################################################
n_inds <- sum(data$status)
length_dat <- subset(data, select = c(age, length))
length_dat <- length_dat %>%
  group_by(age) %>%
  summarise(mean_length=mean(length), sd_length=sd(length))

lengthatage <- ggplot(data=length_dat, aes(x=age, y=mean_length))
lengthatage + geom_point() +
  geom_errorbar(data = length_dat, aes(x=age, ymin = mean_length - sd_length, ymax = mean_length + sd_length)) +
  geom_point(data=ices_l_data, shape=1, aes(x=age, y=length)) +
  geom_line(data=ices_l_data, aes(x=age, y=length)) +
  scale_x_continuous(breaks = seq(0,25,1)) +
  xlab("Age (y)") +
  ylab("Length (cm)") +
  annotate("text", label = n_inds, x= 1, y= 40) +
  theme_classic()
  ggsave(here("R/Rplots/finalpop_lengthatage.pdf"), plot=last_plot())  
  

#Import ICES weight data
  ices_w_data <- read_delim(here("R/ices_dat/ices_weight_data.csv"), 
                                   ";", escape_double = FALSE, trim_ws = TRUE)
ices_weight <- colMeans(subset(ices_w_data, Year_age > "2004"))
ices_weight <- ices_weight[2:13]  
ices_age <- c(3:14)
ices_weight <- as.data.frame(cbind(ices_age,ices_weight))
ices_age_plot <- ggplot(data=ices_weight, aes(x=ices_age, y=ices_weight))
ices_age_plot +
  geom_point() +
  theme_classic()

########################################################################################
#Weight at age
########################################################################################
n_inds <- sum(data$status)

weight_dat <- subset(data, select = c(age, SomWeight, GonWeight))
weight_dat[,4] <- weight_dat[,2] +weight_dat[,3]
colnames(weight_dat)[4] <- "TotWeight"
weight_dat <- weight_dat %>%
  group_by(age) %>%
  summarise(mean_weight=mean(TotWeight), sd_weight=sd(TotWeight))

weightatage <- ggplot(data=weight_dat, aes(x=age, y=mean_weight))
weightatage + geom_point() +
  geom_errorbar(data = weight_dat, aes(x=age, ymin = mean_weight - sd_weight, ymax = mean_weight + sd_weight)) +
  geom_point(data = ices_weight, shape=1, size=4, aes(x=ices_age, y=ices_weight)) +
  geom_line(data = ices_weight, aes(x=ices_age, y=ices_weight)) +
  scale_x_continuous(breaks = seq(0,25,1)) +
  xlab("Age (y)") +
  ylab("Weight (kg)") +
  labs(fill = "Data") +
  annotate("text", label = n_inds, x= 1, y= 15) +
  theme_classic()
ggsave(here("R/Rplots/finalpop_weightatage.pdf"), plot=last_plot())

#Import ICES maturity data
ices_m_data <- read_delim(here("R/ices_dat/ices_maturity_data.csv"), 
                          ";", escape_double = FALSE, trim_ws = TRUE)
ices_mat <- colMeans(subset(ices_m_data, Year_age > "2004"))
ices_mat <- ices_mat[2:13]  
ices_age <- c(3:14)
ices_mat <- as.data.frame(cbind(ices_age,ices_mat))
ices_mat_plot <- ggplot(data=ices_mat, aes(x=ices_age, y=ices_mat))
ices_mat_plot +
  geom_point() +
  theme_classic()

########################################################################################
#Proportion mature plot
########################################################################################
n_inds <- sum(data$status)
model_m_data <- as.data.frame(cbind(data$age, data$maturity))
colnames(model_m_data) <- c("age","maturity")
model_m_data$maturity <- ifelse(model_m_data$maturity>0,1,0)
model_m_data <- model_m_data %>%
  group_by(age) %>%
  summarise(n_mature=sum(maturity), n_tot=n())
model_m_data[,4] <- model_m_data[,2]/model_m_data[,3]
colnames(model_m_data)[4] <- "prop_mature"
propmatureplot <- ggplot(data=model_m_data, aes(x=age, y=prop_mature))
propmatureplot + geom_point() + geom_line() +
  geom_point(data=ices_mat, shape=1, aes(x=ices_age, y=ices_mat)) +
  geom_line(data=ices_mat, linetype = "dashed", aes(x=ices_age, y=ices_mat)) +
  xlab("Age (y)") +
  ylab("Proportion mature") +
  theme_classic()
ggsave(here("R/Rplots/finalpop_propmature.pdf"), plot=last_plot())
