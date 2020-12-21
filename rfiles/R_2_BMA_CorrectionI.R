


# Created 03.03.2020


# Run best fitting models and propagate uncertainty for oversampling


#---------------------------------------------------------------------------


library(sf)
library(INLA)
library(spdep)
library(pbapply)
library(parallel)
library(dplyr)


listOV <- readRDS("data/data_DataSouthAfrica_corI")

# spatial structure
ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")

# Neighboring structure
W.nb <- poly2nb(ShapeSA_mun)
W.adj <- nb2INLA("W.adj", W.nb) 


# priors
hyper.bym <- list(theta1 = list('PCprior', c(1, 0.01)), theta2 = list('PCprior', c(0.5, 0.5)))
hyper.iid <- list(theta = list(prior="pc.prec", param=c(1,0.01)))



# The best fit model

model_CorrectionI <- function(X){
   
  form <- CCcases ~ 1 + offset(log(HIVcases)) +
    
    f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
      constr = TRUE, hyper = hyper.bym) + 
    f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(spatim, model="iid", hyper=hyper.iid, constr = TRUE)
  
  
  in.mod <- 
    inla(formula = form,
         data=listOV[[X]],
         family="poisson",  
         verbose = F, 
         num.threads = 1, 
         control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE), 
         control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
         control.mode=list(restart = T, theta = c(3.5525792, 1.2725765, 0.4515591, 0.2792196, 2.7877063))
    )
  
  return(in.mod)
  
}



# Set up parallel environment

t_0 <- Sys.time()

k <- 1:100
ncores <- 20
cl_inla <- makeCluster(ncores, methods=FALSE)
  
clusterEvalQ(cl_inla, {
    library(INLA)
})
  
par.fun <- function(k){model_CorrectionI(X = k)}

clusterExport(cl_inla, c("model_CorrectionI", "listOV", "W.adj", 
                         "hyper.bym", "hyper.iid", "k", "par.fun"))
  
  
# For 1 parameter use parSapply
outpar <- parLapply(cl = cl_inla, k, par.fun)
  
  
# close parallel environment
stopCluster(cl_inla)

# merge the inla merge files
mod_CorrectionI <- inla.merge(outpar, prob = rep(1, times = 100))

# saveRDS(mod_CorrectionI, file = "results/mod_CorrectionI")

t_1 <- Sys.time()
t_1 - t_0 


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################



