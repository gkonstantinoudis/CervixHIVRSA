

# Created 10.12.2019


# Fit INLA for South Africa




#-----------------------------------------------------------------------------------------

library(maptools)
library(INLA)
library(spdep)
library(pbapply)
library(sf)
library(dplyr)

# load the data


dat <- readRDS("data/data_DataSouthAfrica")

# load the shapefiles
ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")

# Prepare the neighboring matrix
W.nb <- poly2nb(ShapeSA_mun)
nb2INLA("W.adj", W.nb) 

# Set the IDs for the random effects
# ID for age 
dat$id.age <- as.numeric(dat$age)
# id for spatial units to fit the BYM2
dat$id.space <- as.numeric(dat$ID)
# ID for time units to prepare for the temporal rw
dat$id.time <- as.numeric(c(dat$year - 2003))
# ID for spacetime interaction type I
temp <- expand.grid(ID = unique(dat$id.space), year = c(2004:2014))
temp$spatim <- 1:nrow(temp)
dat <- left_join(dat, temp, by = c("ID" = "ID", "year" = "year"))





# set the priors
# Prior for bym2
hyper.bym <- list(theta1 = list('PCprior', c(1, 0.01)), theta2 = list('PCprior', c(0.5, 0.5)))
# Prior of iid and rws
hyper.iid <- list(theta = list(prior="pc.prec", param=c(1,0.01)))


# Type I interaction full model
form <- CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
    constr = TRUE, hyper = hyper.bym) + 
  f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(spatim, model="iid", hyper=hyper.iid, constr = TRUE)




# INLA.call


in.mod <- 
  inla(formula = as.formula(form),
       data=dat,
       family="poisson",  
       verbose = TRUE, 
       control.compute=list(dic=TRUE, waic = TRUE, cpo=TRUE, config = TRUE), 
       control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
       control.mode=list(restart = T)
  )

in.mod <- inla.rerun(in.mod)

# saveRDS(in.mod, file = "results/mod_NoCorrection")


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


