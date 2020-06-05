#Packages
library(tidyverse)
library(skimr)
library(pscl)
library(boot)

#load
mydat <- read.csv("Data/zeroinfl_data.csv", header = T)

#inspect
skim(mydat)

#-----
#plot response variable
#-----

#sculpin - rarest of the three
hist(mydat$Cottus_ab,
     main = "Sculpin Abundance All Sites",
     xlab = "Count",
     breaks = 10)

#var(mydat$Cottus_ab)

sc.pres <- mydat %>% 
  filter(Cottus_ab > 0) #only sites when present

histinfo<-hist(sc.pres$Cottus_ab,
     main = "Sculpin Abundance When Present",
     xlab = "Count",
     breaks = 10)
histinfo

#mean(sc.pres$Cottus_ab)
#sd.cott.pres <- sd(sc.pres$Cottus_ab)
#var(sc.pres$Cottus_ab)


#plot against predictor variables
sculpin <- mydat %>%
  select(Cottus_ab, avgT, mFlow, boulder, HAiFLS_for, med_len, BRT_100m)

plot(sculpin$avgT, sculpin$Cottus_ab,
     col="blue",
     pch=18)
plot(sculpin$mFlow, sculpin$Cottus_ab,
     col="blue",
     pch=18)
plot(sculpin$boulder, sculpin$Cottus_ab,
     col="blue",
     pch=18)
plot(sculpin$HAiFLS_for, sculpin$Cottus_ab,
     col="blue",
     pch=18)
plot(sculpin$med_len, sculpin$Cottus_ab,
     col="blue",
     pch=18)
plot(sculpin$BRT_100m, sculpin$Cottus_ab,
     col="blue",
     pch=18)


#-----
#longnose dace 
#-----
hist(mydat$LND_ab,
     main = "Longnose Dace Abundance All Sites",
     xlab = "Count",
     breaks = 10)

#var(mydat$LND_ab)

lnd.pres <- mydat %>% 
  filter(LND_ab > 0) #only sites when present

histinfo2<-hist(lnd.pres$LND_ab,
               main = "Longnose Dace Abundance When Present",
               xlab = "Count",
               breaks = 10)
histinfo2

#mean(lnd.pres$LND_ab)
#sd(lnd.pres$LND_ab)
#var(lnd.pres$LND_ab)


#plot against predictor variables
longnose <- mydat %>%
  select(LND_ab, avgT, pctcbbl, elev_m, Area_km2, med_len, BRT_100m)

plot(longnose$avgT, longnose$LND_ab,
     col="blue",
     pch=18)
plot(longnose$pctcbbl, longnose$LND_ab,
     col="blue",
     pch=18)
plot(longnose$elev_m, longnose$LND_ab,
     col="blue",
     pch=18)
plot(longnose$Area_km2, longnose$LND_ab,
     col="blue",
     pch=18)
plot(longnose$med_len, longnose$LND_ab,
     col="blue",
     pch=18)
plot(longnose$BRT_100m, longnose$LND_ab,
     col="blue",
     pch=18)




#-----
#southern redbelly dace 
#-----
hist(mydat$SRD_ab,
     main = "Southern Redbelly Dace Abundance All Sites",
     xlab = "Count",
     breaks = 10)

#var(mydat$SRD_ab)

srd.pres <- mydat %>% 
  filter(SRD_ab > 0) #only sites when present

histinfo3<-hist(srd.pres$SRD_ab,
                main = "Southern Redbelly Dace Abundance When Present",
                xlab = "Count",
                breaks = 10)
histinfo3

#mean(srd.pres$SRD_ab)
#sd(srd.pres$SRD_ab)
#var(srd.pres$SRD_ab)

#plot against predictor variables
redbelly <- mydat %>%
  select(SRD_ab, avgT, pctfines, avdep, med_len, BRT_100m)

plot(redbelly$avgT, redbelly$SRD_ab,
     col="blue",
     pch=18)
plot(redbelly$avdep, redbelly$SRD_ab,
     col="blue",
     pch=18)
plot(redbelly$pctfines, redbelly$SRD_ab,
     col="blue",
     pch=18)
plot(redbelly$med_len, redbelly$SRD_ab,
     col="blue",
     pch=18)
plot(redbelly$BRT_100m, redbelly$SRD_ab,
     col="blue",
     pch=18)

#-----
#Back to Powerpoint
#-----





#######################
#   Zero-infl models
#######################

#-----
# Just Longnose Dace to save time
#-----

#Longnose Dace Count = response, segment length = offset
#zero-inflated negative binomial model 

lnd.full.mod <- zeroinfl(LND_ab ~ avgT+pctcbbl+elev_m+Area_km2+med_len+BRT_100m | 1,
                         data = mydat,
                         dist = "negbin",
                         offset = log(SegLen))
summary(lnd.full.mod)


#-----
#extracting confidence intervals for the parameters

#top model 
dput(round(coef(lnd.full.mod, "count"), 4)) #count process
dput(round(coef(lnd.full.mod, "zero"), 4)) #extra zero process

f <- function(data, i) {
  require(pscl)
  m <- zeroinfl(LND_ab ~ avgT+pctcbbl+elev_m+Area_km2+med_len+BRT_100m | 1,
                data = data[i, ],
                dist = "negbin",
                offset = log(SegLen),
                start = list(count = c(-17.983,1.3192,0.0394,-0.0354,
                                       0.0042,-0.004,-0.044),
                             zero = c(-8.8975)))
  as.vector(t(do.call(rbind, coef(summary(m)))[, 1:2]))
}

set.seed(10)
(res <- boot(mydat, f, R = 1000))

## basic parameter estimates with percentile and bias adjusted CIs
parms <- t(sapply(c(1, 3, 5, 7, 9, 11, 13, 17), function(i) {
  out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"))
  with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
              bcaLL = bca[4], bcaUL = bca[5]))
}))

## add row names
row.names(parms) <- names(coef(lnd.full.mod))
## print results
parms

## compare with normal based approximation
confint(lnd.full.mod)

## exponentiated parameter estimates with percentile and bias adjusted CIs
expparms <- t(sapply(c(1, 3, 5, 7, 9, 11, 13, 17), function(i) {
  out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"), h = exp)
  with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
              bcaLL = bca[4], bcaUL = bca[5]))
}))

## add row names
row.names(expparms) <- names(coef(lnd.full.mod))
## print results
expparms


#########################################################################
library(emmeans)
# emmeans on continuous predictors
# See Russ Lenth's long response to this question:
# https://stackoverflow.com/questions/52381434/emmeans-continuous-independant-variable
# Also the "basics" vignette to emmeans

ref_grid(lnd.full.mod) #these are the mean values for all the covariates

#
skim(mydat)

# Plot at quantile values
emmip(lnd.full.mod, BRT_100m~avgT, at = list(avgT = c(9.68,14.4,15.58,16.78,19.83), 
                                              BRT_100m = c(0,7.33,96.14), type = "response"))  +
  theme_bw()+
  theme(panel.grid = element_blank())


# Plot effect of temp only (other covariates at their mean)
temp_rg <- ref_grid(lnd.full.mod, at = list(avgT = mydat$avgT))
temp_val <- temp_rg@grid$avgT
temp_pred <- predict(temp_rg, type = "response")

plot_df_temp <- data.frame(avgT = temp_val, effect = temp_pred)
ggplot(plot_df_temp, aes(x = avgT, y = effect)) + geom_line() +
  xlab("avwid") + ylab("Predicted effect") +
  labs(title = "Predicted effect of mean stream temp on LND Count",
       subtitle = "Other covariates held at their mean values")


# Plot predicted vs observed
lnd.df <- mydat %>%
  select(avgT, pctcbbl, Area_km2, elev_m, med_len, BRT_100m, SegLen)

all_rg <- predict(lnd.full.mod, newdata = lnd.df)# This gives you predictions at your data points

plot_df_compare <- data.frame(obs_value = mydat$LND_ab, estimate = all_rg)
ggplot(data = plot_df_compare, aes(x = estimate, y = obs_value)) + geom_point() +
  xlab("Predicted Count") + ylab("Observed Count")+
  scale_y_continuous(limits = c(0,200))+
  scale_x_continuous(limits = c(0,200))
#The line of points at zero is the zero-inflated part of the model.

ggplot(data = plot_df_compare, aes(x = estimate, y = obs_value)) + geom_point() +
  xlab("Predicted Count") + ylab("Observed Count")


df_outlier <- data.frame(obs_value = mydat$LND_ab,
                         estimate = all_rg,
                         temp = mydat$avgT,
                         cobble = mydat$pctcbbl,
                         area = mydat$Area_km2,
                         elev = mydat$elev_m,
                         length = mydat$med_len,
                         trout = mydat$BRT_100m)

outlier <- df_outlier %>%
  filter(estimate >200) #"ideal" conditions based on relationships observed with covariates

cor(mydat$LND_ab, all_rg) #0.09 -- hm. nice



#########################################################################
lnd.df$estimate <- predict(lnd.full.mod, newdata = lnd.df)# This gives you predictions at your data points
lnd.df$resid <- residuals(lnd.full.mod)

ggplot(lnd.df, aes(x = avgT, y = mydat$LND_ab)) +
  geom_segment(aes(xend = avgT, yend = estimate), alpha = .2) +  # Lines to connect points
  geom_point() +  # Points of actual values
  geom_point(aes(y = estimate), shape = 1, color = "blue") +  # Points of predicted values
  theme_bw()+
  scale_y_continuous(limits = c(0,180))
#########################################################################

