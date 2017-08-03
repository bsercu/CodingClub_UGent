#Script for resampling on Irene issue

setwd("/home/lionel/Documents/PostDoc_Ghent/Course/CodingClub/")
#load data
dat <- read.table("CodingClub_UGent/Session_20170803/Irenes_data.txt")

#load libraries
library(plyr)

#a function to fit the models and return their AIC plus a sim number
fit_i <- function(sim,data_resamp){
  m1 <- glm(abun ~ fragmean + total_all_ba,data_resamp,family = poisson)
  m2 <- glm(abun ~ -1 + rel_qrob + rel_qrub + rel_fsyl + fragmean + total_all_ba,data_resamp,family = poisson)
  m3 <- glm(abun ~ -1 + rel_qrob + rel_qrub + rel_fsyl + fragmean + total_all_ba +  rel_qrob:rel_qrub + rel_qrob:rel_fsyl + rel_qrub:rel_fsyl,data_resamp,family = poisson)
  aic_m <- AIC(m1,m2,m3)
  out <- data.frame(Sim=sim,Best=paste0("Model_",which.min(aic_m$AIC)))
  return(out)
}


#a function to sample leaves from a specific tree
sample_i <- function(treeID){
  subst <- subset(dat,tree_id == treeID)
  n_sample <- c(48,24,12)[unique(subst$specrich)]
  return(subst[sample(1:dim(subst)[1],n_sample,replace=FALSE),])
}



#a function putting together these two functions
compute_m <- function(sim){
  data_resamp <- rbind.fill(lapply(tree_name,function(treeID) sample_i(treeID)))
  fits <- fit_i(sim,data_resamp)
  return(fits)
}

#all trees ID
tree_name <- levels(dat$tree_id)
#do this 10 time for now
system.time(pre_res <- rbind.fill(lapply(1:1000,function(sim) compute_m(sim))))
table(pre_res$Best) #Model 2 seems better

#my solution to this issue would be, first to aggregate at plot level
library(dplyr)
dat %>%
  group_by(id_plot,specrich,total_all_ba,rel_fsyl,rel_qrob,rel_qrub,fragmean)%>%
  summarise(Abun=sum(abun),Weight=n()) -> dat_dd

mo <- glm.nb(Abun ~ fragmean + total_all_ba+offset(Weight),dat_dd,family = poisson)



