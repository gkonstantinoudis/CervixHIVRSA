


# Calculate age standardized incidence rate


# Created 01.03.2020


#---------------------------------------------------------------------------


library(INLA)
library(sf)
library(dplyr)
library(patchwork)
library(ggplot2)
library(viridis)


# load data
dat <- readRDS("data/data_DataSouthAfrica")


# load the shapefiles
ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")



# Calculate incidence rate age standardized using the World population weights


WorldPopulation <- data.frame(age.group = c("[0,5)", "[5,10)", "[10,15)", "[15,20)", 
                                            "[20,25)", "[25,30)", "[30,35)", "[35,40)", 
                                            "[40,45)", "[45,50)", "[50,55)", "[55,60)", 
                                            "[60,65)", "[65,70)", "[70,75)", "[75,80)", ">80"), 
                              segi = c(12, 10, 9, 9, 8, 8, 6, 6, 6, 6,5,4,4,3,2,1,1), 
                              WHO = c(8.86, 8.69, 8.60, 8.47, 8.22, 7.93, 7.61, 7.15, 6.59, 6.04, 5.37,
                                      4.55, 3.72, 2.96, 2.21, 1.52, 1.54))

#------- Temporal Incidence Rate




# A function to make the integration wrt space and perform the World age standardization 

TemporalIncidenceRate <- function(X){
  
  temp_mod <- X
  # merge with the population weights
  temp_mod <- merge(temp_mod, WorldPopulation, by.x = "age", by.y = "age.group", all.x = TRUE)
  temp_mod$segi <- as.numeric(temp_mod$segi)/100 
  
  # calculate the rate with the world population weights
  temp_plot <- cbind(temp_mod[,-c(7:13)], 
                     (temp_mod[,c("fitted_median", "fitted_025", "fitted_975")]/temp_mod$HIVcases)*
                       temp_mod$segi)
  
  # sum over the age categories to get it adjusted for age
  temp_mod2 <- aggregate(temp_plot[,c("fitted_median", "fitted_025", "fitted_975")], 
                         by = list(year = temp_plot$year, space = temp_plot$ID), sum)
  
  # take the mean over space
  temp_mod3 <- aggregate(temp_mod2[,c("fitted_median", "fitted_025", "fitted_975")], 
                         by = list(year = temp_mod2$year), mean)
  
  ret <-cbind(temp_mod3$year, temp_mod3[,c("fitted_median", "fitted_025", "fitted_975")]*100000)
  
  return(ret)
  
  
}



listl <- list()

# 1. No correction


mod_1 <- readRDS("results/mod_NoCorrection")

dat$fitted_median <- mod_1$summary.fitted.values$`0.5quant`
dat$fitted_025 <- mod_1$summary.fitted.values$`0.025quant`
dat$fitted_975 <- mod_1$summary.fitted.values$`0.975quant`
dat$offset <- dat$HIVcases

# store the results as the first element of a list
listl[[1]] <- TemporalIncidenceRate(dat)
  



# First we will calculate correction II since it is directly applied to the posteriors we sampled from


# 2. Correction II

mod_3 <- readRDS("results/mod_CorrectionII")

# and apply the correction on the posteriors
dat$fitted_median <- mod_3$summary.fitted.values$`0.5quant`
dat$fitted_025 <- mod_3$summary.fitted.values$`0.025quant`
dat$fitted_975 <- mod_3$summary.fitted.values$`0.975quant`
dat$offset <- dat$HIVcases


listl[[2]] <- TemporalIncidenceRate(dat)


# 3. Correction I

mod_2 <- readRDS("results/mod_CorrectionI")

# transform the marginal to get the fitted
lp_marginal.exp_CorI <- lapply(X = mod_2$marginals.linear.predictor, FUN = function(x) inla.tmarginal(exp, marginal = x))
# get summary statistics from the marginals
sum.stats_CorI <- lapply(X = lp_marginal.exp_CorI, FUN = function(x) inla.zmarginal(marginal = x, silent = TRUE))

dat$fitted_median <- sapply(sum.stats_CorI, function(X) X$quant0.5)
dat$fitted_025 <- sapply(sum.stats_CorI, function(X) X$quant0.025)
dat$fitted_975 <- sapply(sum.stats_CorI, function(X) X$quant0.975)
dat$offset <- dat$HIVcases

listl[[3]] <- TemporalIncidenceRate(dat)



# 4. Correction IV

mod_4 <- readRDS("results/model_FullCorrection")

# transform the marginal to get the fitted
lp_marginal.exp_Corfull <- lapply(X = mod_4$marginals.linear.predictor, FUN = function(x) inla.tmarginal(exp, marginal = x))
# get summary statistics from the marginals
sum.stats_CorFull <- lapply(X = lp_marginal.exp_Corfull, FUN = function(x) inla.zmarginal(marginal = x, silent = TRUE))

dat$fitted_median <- sapply(sum.stats_CorFull, function(X) X$quant0.5)
dat$fitted_025 <- sapply(sum.stats_CorFull, function(X) X$quant0.025)
dat$fitted_975 <- sapply(sum.stats_CorFull, function(X) X$quant0.975)
dat$offset <- dat$HIVcases

listl[[4]] <- TemporalIncidenceRate(dat)



nam <- c("No Correction", "Correction II", "Correction I", "Full Correction")
# and the plot

listl <- lapply(listl, function(X){
  colnames(X) <- c("year", "median", "LL", "UP")
  return(X)
})



listplot <- list()
listval <- list()


for(i in 1:4){
  
  datNoCor_1 <- listl[[i]]
  listval[[i]] <- datNoCor_1
  
  listplot[[i]] <- ggplot(datNoCor_1, aes(x=year)) +  
    geom_line(aes(x = year, y=median), col = viridis(n = 1)) + 
    geom_point(aes(x = year, y=median), size=1.5, shape=19, col = viridis(n = 1))  +
    geom_ribbon(aes(ymin=LL, ymax = UP), alpha = .2, fill = viridis(n = 1)) +
    xlab("Year") + ylab("Incidence Rate") + 
    theme_bw() +  ggtitle(nam[i]) + theme(plot.margin = unit(c(0,0,0,0), "cm")) + ylim(0, 900) + 
    theme(text = element_text(size=9))
  
}





png("results/Fig2.png", width =17, height = 10, units ="cm", res = 300)
print(
  (listplot[[1]]|listplot[[3]])/
    (listplot[[2]]|listplot[[4]])
)
dev.off()






temp.list <- lapply(listl, function(X){
  
  ret <- data.frame(median = round(X$median), CrI = paste0("(", round(X$LL), ", ", round(X$UP), ")"))
  return(ret)
}
)

temp.list <- do.call(cbind, temp.list)
rownames(temp.list) <- 2004:2014
library(xtable)
xtable(temp.list)

xtable(t(as.data.frame(apply(temp.list, 2, median))))



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################






