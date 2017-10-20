library(ape)
library(caper)
library(car)
library(nlme)
library(tidyverse)

phy <- read.tree("data/stein-et-al-single.tree")
dat <- read.csv("data/chond-data.csv", 
                stringsAsFactors = FALSE,header = TRUE)

# Creating comparative.data object for use in pgls() function
cd <- comparative.data(phy, dat, names.col = "species")

### Fitting models

# mass only
zslope_noutc_phy_2 <- pgls(log(rmax) ~ log_wt, data = cd, lambda = "ML")

# mass + temp
zslope_invtemp_phy_scaled_2 <- pgls(log(rmax) ~ log_wt + invtemp_scaled,
                                    data = cd, lambda = "ML")

# mass + depth
zslope_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt + depth_scaled, 
                                  data = cd, lambda = "ML")

# mass + temp + depth
zslope_invtemp_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt + invtemp_scaled 
                                          + depth_scaled,
                                          data = cd, lambda = "ML")

# mass + temp * depth
zslope_invtemp_x_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt + 
                                              invtemp_scaled * depth_scaled, 
                                            data = cd, lambda = "ML")

#mass * depth + temp 
zslope_invtemp_wt_x_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt * 
                                                 depth_scaled + invtemp_scaled, 
                                               data = cd, lambda = "ML")

# mass * depth
zslope_wt_x_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt * depth_scaled, 
                                       data = cd, lambda = "ML")

# mass * temp 
zslope_wt_x_invtemp_phy_scaled_2 <- pgls(log(rmax) ~ log_wt * invtemp_scaled, 
                                         data = cd, lambda = "ML")

# mass * temp + depth
zslope_wt_x_invtemp_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt * 
                                                 invtemp_scaled + depth_scaled, 
                                               data = cd, lambda = "ML")

# testing model with 2 interactions: mass * temp + mass * depth 
zslope_wt_x_invtemp_wt_x_depth_phy_scaled_2 <- pgls(log(rmax) ~ log_wt * 
                                                      invtemp_scaled + 
                                                      log_wt * depth_scaled, 
                                                    data = cd, lambda = "ML")

#  saving all models in list
capermodels <- list(zslope_noutc_phy_2, 
                    zslope_depth_phy_scaled_2,
                    zslope_invtemp_phy_scaled_2,
                    zslope_invtemp_depth_phy_scaled_2,
                    zslope_invtemp_x_depth_phy_scaled_2,
                    zslope_wt_x_depth_phy_scaled_2,
                    zslope_wt_x_invtemp_phy_scaled_2,
                    zslope_invtemp_wt_x_depth_phy_scaled_2,
                    zslope_wt_x_invtemp_depth_phy_scaled_2,
                    zslope_wt_x_invtemp_wt_x_depth_phy_scaled_2) 

# manually extracting information from each model object
aics <- sapply(capermodels, function (x) bbmle::AIC(x)) 
aiccs <- sapply(capermodels, function (x) x$aicc) 
formulas <- sapply(capermodels,   # tidying formulas for easier reading
                   function (x) deparse(formula(x), width.cutoff = 90L)) %>%
  sub("^log\\(\\w+\\)\\s\\~\\s", "", .) %>%
  gsub("_scaled", "", .) %>%
  gsub("\\_wt", "(M)", .)
r2s <- sapply(capermodels, function (x) summary(x)$r.squared) 
ar2s <- sapply(capermodels, function (x) summary(x)$adj.r.squared)
LL <- sapply(capermodels, function (x) x$model$log.lik) 
ks <- sapply(capermodels, function (x) x$k) 

models_table <- data.frame(formulas, ks, LL, aics, aiccs, r2s, ar2s) %>%
  rename(Model = formulas, n = ks, AIC = aics,
         AICc = aiccs, R_sq = r2s, adj_R_sq = ar2s) %>%
  mutate(Model = as.character(formulas), 
         LL = round(LL, 1),
         AIC = round(AIC, 1), AICc = round(AICc, 1), 
         dAIC = round(AIC - min(AIC), 2),
         dAICc = round(AICc - min(AICc), 2),
         R_sq = round(R_sq, 2), adj_R_sq = round(adj_R_sq, 2),
         Weights = round(exp(-dAICc/2)/sum(exp(-dAICc/2)), 3)) 
models_table

### Top-ranked model diagnostic plots and analyses
topmod <- which(models_table$dAICc == 0)
bestmodel <- capermodels[[topmod]]

# Re-running the model using nlme:gls for diagnostics
# parameter estimates are obtained
# corPagel() is used as it's equivalent to `lambda = "ML"` in pgls
bestmodel_gls <- nlme::gls(formula(bestmodel), data = cd$data, 
                           correlation = corPagel(1, cd$phy))

# almost identical coefficent estimates
rbind(coef(bestmodel), coef(bestmodel_gls))

# residual plots almost identical as well, homoscedastic
plot(resid(bestmodel) ~ fitted(bestmodel)) # pgls
plot(bestmodel_gls) # gls

# Variance inflation factors (on gls model)
vif(bestmodel_gls)

# qqplot (on gls model) looks normal
qqnorm(bestmodel_gls)

# profile plot of lambda
plot(pgls.profile(bestmodel))
