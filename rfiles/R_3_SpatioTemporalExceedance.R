

# Calculate spatiotemporal exceedance probability


# Created 01.03.2020


#---------------------------------------------------------------------------


library(INLA)
library(patchwork)
library(ggsn)
library(dplyr)
library(sf)
library(ggplot2)


# load the model without correction
mod <- readRDS("results/mod_NoCorrection")

# load data
dat <- readRDS("data/data_DataSouthAfrica")

# load shapefiles
ShapeSA_mun <- read_sf("data/ShapeSA_mun.shp")

# Calculate the marginals for the latent spatial field and take the exp transformation to get it at the RR scale
t_0 <- Sys.time()
LF_marginal.exp <- lapply(X= mod$marginals.random$spatim, FUN = function(x) inla.tmarginal(exp, marginal = x))

# calcuate exceedance probability that RR>1

exceed.prob_1 <- lapply(X= LF_marginal.exp, FUN = function(x) inla.pmarginal(marginal = x, 1))
exceed.prob_1 <- 1 - unlist(exceed.prob_1)
boxplot(exceed.prob_1)

t_1 <- Sys.time()
t_1 - t_0




# and the plot 

lst <- list()
lst.ex <- list()
year <- 2004:2014

plotSpatTime <- list()
start <- seq(from = 1, to = 1859, by = 169)



for(i in 1:11){
  
  dat$Pr <- exceed.prob_1
  dat$year <- as.character(dat$year)
  
  ShapeSA_mun$Pr <- exceed.prob_1[start[i]:c(169*i)]
  ShapeSA_mun$Pr_cat <- NA
  ShapeSA_mun$Pr_cat[ShapeSA_mun$Pr <0.2] <- "[0, 0.2]"
  ShapeSA_mun$Pr_cat[ShapeSA_mun$Pr > 0.8] <- "(0.8, 1]"
  ShapeSA_mun$Pr_cat[is.na(ShapeSA_mun$Pr_cat)] <- "(0.2, 0.8]"
  ShapeSA_mun$Pr_cat <- factor(ShapeSA_mun$Pr_cat, levels = c("[0, 0.2]", "(0.2, 0.8]", "(0.8, 1]"))
  
  plotSpatTime[[i]] <- ggplot() + geom_sf(data = ShapeSA_mun,  aes(fill = Pr_cat)) + 
    scale_fill_viridis_d(name = "") + theme_bw() + 
    ggtitle(paste0("Pr(RR>1) in ", c(2003+i))) + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    north(ShapeSA_mun, location = "bottomright", symbol = 12) + 
    theme(text = element_text(size=9))
  
}



lst <- plotSpatTime

png("SpatTimeInteractionRR.png", width =30, height = 25, units ="cm", res = 300)
print(
  (lst[[1]]|lst[[2]]|lst[[3]])/
    (lst[[4]]|lst[[5]]|lst[[6]])/
    (lst[[7]]|lst[[8]]|lst[[9]])/
    (lst[[10]]|lst[[11]])
)
dev.off()

















# Get the overall spatial trend (Fig 4)

# Model no correction no covariates

ShapeSA_mun$ranef <- exp(mod$summary.random$id.space$`0.5quant`)[1:169]


p1 <- ggplot() + geom_sf(data = ShapeSA_mun,  aes(fill = ranef)) + 
  scale_fill_viridis_c(name = "", limits = c(0, 4.5)) + theme_bw() + 
  ggtitle("Spatial Relative Risk") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  north(ShapeSA_mun, location = "bottomright", symbol = 12) + 
  theme(text = element_text(size=9))

# and the exceedance

SP_marginal.exp <- lapply(X= mod$marginals.random$id.space, 
                           FUN = function(x) inla.tmarginal(exp, marginal = x))

# calculate exceedance probability that the RR>1

exceed.prob <- lapply(X= SP_marginal.exp, FUN = function(x) inla.pmarginal(marginal = x, 1))
exceed.prob <- 1 - unlist(exceed.prob)

ShapeSA_mun$Pr <- exceed.prob[1:169]


ShapeSA_mun$Pr_cat <- NA
ShapeSA_mun$Pr_cat[ShapeSA_mun$Pr <0.2] <- "[0, 0.2]"
ShapeSA_mun$Pr_cat[ShapeSA_mun$Pr > 0.8] <- "(0.8, 1]"
ShapeSA_mun$Pr_cat[is.na(ShapeSA_mun$Pr_cat)] <- "(0.2, 0.8]"
ShapeSA_mun$Pr_cat <- factor(ShapeSA_mun$Pr_cat, levels = c("[0, 0.2]", "(0.2, 0.8]", "(0.8, 1]"))


p2 <- ggplot() + geom_sf(data = ShapeSA_mun,  aes(fill = Pr_cat)) + 
  scale_fill_viridis_d(name = "") + theme_bw() + 
  ggtitle(paste0("Pr(RR>1)")) + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  north(ShapeSA_mun, location = "bottomright", symbol = 12) + 
  theme(text = element_text(size=9))





# Model no correction adjusted for covariates
modcov <- readRDS("results/mod_covariates")

ShapeSA_mun$ranefcov <- exp(modcov$summary.random$id.space$`0.5quant`)[1:169]


p3 <- ggplot() + geom_sf(data = ShapeSA_mun,  aes(fill = ranefcov)) + 
  scale_fill_viridis_c(name = "", limits = c(0, 4.5)) + theme_bw() + 
  ggtitle("Spatial Relative Risk") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  north(ShapeSA_mun, location = "bottomright", symbol = 12) + 
  theme(text = element_text(size=9))

# and the exceedance

SP_marginal.exp_cov <- lapply(X= modcov$marginals.random$id.space, 
                          FUN = function(x) inla.tmarginal(exp, marginal = x))

# calculate exceedance probability that the RR>1

exceed.prob_cov <- lapply(X= SP_marginal.exp_cov, FUN = function(x) inla.pmarginal(marginal = x, 1))
exceed.prob_cov <- 1 - unlist(exceed.prob_cov)

ShapeSA_mun$Prcov <- exceed.prob_cov[1:169]


ShapeSA_mun$Pr_cat_cov <- NA
ShapeSA_mun$Pr_cat_cov[ShapeSA_mun$Prcov <0.2] <- "[0, 0.2]"
ShapeSA_mun$Pr_cat_cov[ShapeSA_mun$Prcov > 0.8] <- "(0.8, 1]"
ShapeSA_mun$Pr_cat_cov[is.na(ShapeSA_mun$Pr_cat_cov)] <- "(0.2, 0.8]"
ShapeSA_mun$Pr_cat_cov <- factor(ShapeSA_mun$Pr_cat_cov, levels = c("[0, 0.2]", "(0.2, 0.8]", "(0.8, 1]"))


p4 <- ggplot() + geom_sf(data = ShapeSA_mun,  aes(fill = Pr_cat_cov)) + 
  scale_fill_viridis_d(name = "") + theme_bw() + 
  ggtitle(paste0("Pr(RR>1)")) + theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  north(ShapeSA_mun, location = "bottomright", symbol = 12) + 
  theme(text = element_text(size=9))






png("Fig4.png", width = 17, height = 13, units ="cm", res = 300)
(p1|p2)/(p3|p4)
dev.off()



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


