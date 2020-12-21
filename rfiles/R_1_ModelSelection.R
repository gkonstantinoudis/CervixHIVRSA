

# Created 10.12.2019


# Fit INLA for South Africa




#-----------------------------------------------------------------------------------------

library(maptools)
library(INLA)
library(spdep)
library(pbapply)
library(sf)
library(dplyr)
library(xtable)

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

# define the different models to be fit

formulas <- list(
  
  # 1. Model with an intercept
  CCcases ~ 1 + offset(log(dat$HIVcases)),
  
  # 2. Model with an age component
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE),
  
  # 3. Model with a time component
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE),
  
  # 4. Model with a space component
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym), 
  
  # 5. Model with age and time components
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE),
  
  # 6. Model with age and space components
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym), 
  
  # 7. Model with time and space components
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym), 
  
  # 8. Model for time and age and space components
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym) + 
                                            f(id.age, model='rw2', hyper=hyper.iid, constr = TRUE, scale.model = TRUE), 
  
  # Type I interaction spacetime interaction without the age effect
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym) + 
                                            f(spatim, model="iid", hyper=hyper.iid, constr = TRUE),

  # Type I interaction full model
  CCcases ~ 1 + offset(log(dat$HIVcases)) + f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                              constr = TRUE, hyper = hyper.bym) + 
                                            f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                                            f(spatim, model="iid", hyper=hyper.iid, constr = TRUE)

  )




# INLA.call

mod.list <- lapply(formulas, function(X){
  
  in.mod <- 
  inla(formula = as.formula(X),
       data=dat,
       family="poisson",  
       verbose = TRUE, 
       control.compute=list(dic=TRUE, waic = TRUE, cpo=TRUE, config = TRUE), 
       control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
       control.mode=list(restart = T)
  )
  
  return(inla.rerun(in.mod))
  
})


names(mod.list) <- c("int", "age", "time", "space", "agetime", "agespace", "timespace", "timeagespace", 
                     "timespaceT1", "timeagespaceT1")

# saveRDS(mod.list, file = "ModelSelection")


# For the DIC, WAIC and CPO

options(scipen=999)

xtable(
data.frame(
"DIC" = sapply(mod.list, function(X)X$dic$dic),
"WAIC" = sapply(mod.list, function(X)X$waic$waic),
"CPO" = round(sapply(mod.list, function(X)-mean(log(X$cpo$cpo))), digits = 3)
)
)



#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


