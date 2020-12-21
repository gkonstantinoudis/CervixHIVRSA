

# Created 20.03.2020


# Include covariates on the spatiotemporal interaction type IV


#######################################################################################################

library(maptools)
library(INLA)
library(spdep)
library(pbapply)
library(sf)
library(dplyr)

# load the data


dat <- readRDS("data/data_DataSouthAfrica_covar")

# load the shapefiles
ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")

# Prepare the neighboring matrix
W.nb <- poly2nb(ShapeSA_mun)
nb2INLA("W.adj", W.nb) 





# set the priors
# Prior for bym2
hyper.bym <- list(theta1 = list('PCprior', c(1, 0.01)), theta2 = list('PCprior', c(0.5, 0.5)))
# Prior of iid and rws
hyper.iid <- list(theta = list(prior="pc.prec", param=c(1,0.01)))






# INLA model


finalmod_uni1 <- CCcases ~ 1 + offset(log(HIVcases)) + Urbanicity +
                               f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                 constr = TRUE, hyper = hyper.bym) + 
                               f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 

finalmod_uni2 <- CCcases ~ 1 + offset(log(HIVcases)) + SES_cat +
                               f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                 constr = TRUE, hyper = hyper.bym) + 
                               f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 

finalmod_uni3 <- CCcases ~ 1 + offset(log(HIVcases)) + NumberFacilitiesDec +
                               f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                 constr = TRUE, hyper = hyper.bym) + 
                               f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 

finalmod <- CCcases ~ 1 + offset(log(HIVcases)) + Urbanicity + SES_cat + NumberFacilitiesDec +
                               f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
                                 constr = TRUE, hyper = hyper.bym) + 
                               f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
                               f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 





list.form <- list(finalmod_uni1, finalmod_uni2, finalmod_uni3, finalmod)

INLA.RUN <- function(X){
  
  in.mod <- 
    inla(formula = as.formula(X),
         data=dat,
         family="poisson",  
         verbose = FALSE, 
         control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE), 
         control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
         control.mode=list(restart = T, theta = c(3.5525792, 1.2725765, 0.4515591, 0.2792196, 2.7877063))
    )
  
  in.mod <- inla.rerun(in.mod)
  return(in.mod)
}

list.inla <- lapply(list.form, INLA.RUN)

# and store the results

# saveRDS(list.inla[[4]], file = "results/mod_covariates")





# Get the Table of the results

CrI <- function(X){
  tmp <- exp(X$summary.fixed[-1, c("0.5quant", "0.025quant", "0.975quant")])
  datfr <- data.frame("Median" = round(tmp[,1], digits = 2), "95% CrI" = paste0("(", round(tmp[,2], digits = 2),
                                                                                ", ", round(tmp[,3], digits = 2), ")"))
  return(datfr)
}

pI <- rbind(rbind(c("1", "-"), CrI(list.inla[[1]])), 
      rbind(c("1", "-"), CrI(list.inla[[2]])), 
      rbind(c("1", "-"), CrI(list.inla[[3]])))

pII <- CrI(list.inla[[4]])

pII <- rbind(rbind(c("1", "-"), pII[1,]), 
      rbind(c("1", "-"), pII[c(2:10),]), 
      rbind(c("1", "-"), pII[c(11:19),]))

# report also exceedance probabilities

ExProb <- function(X){
  
  ranef_marginal.exp <- lapply(X= X$marginals.fixed, FUN = function(x) inla.tmarginal(exp, marginal = x))
  exceed.prob_1 <- lapply(X= ranef_marginal.exp, FUN = function(x) inla.pmarginal(marginal = x, 1))
  exceed.prob_1 <- 1 - unlist(exceed.prob_1)
  
  return(round(exceed.prob_1[-1], digits = 2))
  
}

pIex <- c(c(1,ExProb(list.inla[[1]])), 
          c(1,ExProb(list.inla[[2]])), 
          c(1,ExProb(list.inla[[3]])))

pIIex <- ExProb(list.inla[[4]])

pIIex <- c(c(1, pIIex[1]), 
           c(1, pIIex[c(2:10)]), 
           c(1, pIIex[c(11:19)]))




write.csv(cbind(pI,pIex, pII, pIIex), "resCovariates.csv")






##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# 
# 
# tmp <- dat
# tmp <- aggregate(tmp$HIVcases, by = list(ID2 = tmp$ID2), sum)
# tmp2 <- dat[,c("ID2", "NumberFacilities")]
# tmp2 <- tmp2[!duplicated(tmp2$ID2),]
# 
# 
# 
# tmp2 <- left_join(tmp2, tmp)
# tmp2$rate <- (tmp2$NumberFacilities/tmp2$x)*10000
# 
# 
# 
# ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")
# 
# tmp2 <- left_join(ShapeSA_mun, tmp2, by = c("ID" = "ID2"))
# ggplot() + geom_sf(data = tmp2, aes(fill = rate), col = NA)
# 
# head(tmp2)
# tmp2 <- tmp2[,c("ID", "rate")]
# tmp2$geometry <- NULL
# colnames(tmp2)[2] <- "FacPerHIV"
# 
# 
# dat <- left_join(dat, tmp2, by = c("ID2" = "ID"))
# 
# 
# finalmod <- CCcases ~ 1 + offset(log(HIVcases)) + Urbanicity + SES_cat + FacPerHIV +
#   f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
#   f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
#     constr = TRUE, hyper = hyper.bym) + 
#   f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
#   f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 
# 
# 
# 
# in.mod <- 
#   inla(formula = as.formula(finalmod),
#        data=dat,
#        family="poisson",  
#        verbose = FALSE, 
#        control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE), 
#        control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
#        control.mode=list(restart = T, theta = c(3.5525792, 1.2725765, 0.4515591, 0.2792196, 2.7877063))
#   )
# 
# 
# in.mod$summary.fixed
# 
# exp(-0.00760603)
# exp(0.001806357)
# exp(0.0103832)
# 
# 
# dat$FacPerHIV_cat <- cut(dat$FacPerHIV, breaks = quantile(dat$FacPerHIV, probs = seq(from = 0, to = 1, by = .1)), 
#                          include.lowest = TRUE, labels = 1:10)
# 
# 
# 
# finalmod <- CCcases ~ 1 + offset(log(HIVcases)) + Urbanicity + SES_cat + FacPerHIV_cat +
#   f(id.time, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
#   f(id.space, model='bym2', graph='W.adj', scale.model = TRUE, 
#     constr = TRUE, hyper = hyper.bym) + 
#   f(id.age, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) + 
#   f(spatim, model="iid", hyper=hyper.iid, constr = TRUE) 
# 
# 
# 
# in.mod <- 
#   inla(formula = as.formula(finalmod),
#        data=dat,
#        family="poisson",  
#        verbose = FALSE, 
#        control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE), 
#        control.inla=list(strategy="simplified.laplace", int.strategy="eb"),
#        control.mode=list(restart = T, theta = c(3.5525792, 1.2725765, 0.4515591, 0.2792196, 2.7877063))
#   )
# 
# 
# in.mod$summary.fixed
# dat
