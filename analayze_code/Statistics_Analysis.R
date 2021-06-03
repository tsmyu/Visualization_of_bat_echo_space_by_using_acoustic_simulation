rm(list=ls())
setwd("D:/Publications/3_Submitted/Teshima/Statistics")

#########################################################################
## Effect of spatial learning on the echo incidence point distribution ##
#########################################################################

## libraries
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(sjPlot)
library(car)
library(effects)
library(emmeans)
library(ggplot2)

# data
d <- read.table("210304_data.txt",header=T,na.string="NA")
str(d)
d$B <- as.factor(d$bat_id)
d$P <- as.factor(d$pulse_id)
d$dist <- d$distance_to_edge_cm
d$fl <- as.factor(d$flight_id)

d$f <- recode(d$fl, "0" = "First","1"="Last", .default = levels(d$fl))

table(d$fl,d$B)
table(d$P,d$B)

# selecting data from the inner half of the wall - 0-33 cm
d_low <- d[which(d$dist<34),]

hist(d_low$dist)
summary(d_low$dist)

dotchart(d_low$dist,groups=d_low$fl)
# no outliers

d_low_1st <- d_low[(1:587),]
d_low_last <- d_low[(588:1058),]

# Raw data plot
boxplot(d_low_1st$dist~d_low_1st$fl)
boxplot(d_low_last$dist~d_low_last$fl,add=T)

# Model
m3_low <- glmmTMB(dist ~ f + (1|B/P),family="nbinom1",data=d_low) 

# Quality check
simulationOutput_m3_low <- simulateResiduals(fittedModel = m3_low, plot = T,n=1000)
testResiduals(simulationOutput_m3_low)
plotResiduals(simulationOutput_m3_low, d_low$fl)

# random effects:
plot_model(m3_low,type="re")

# overall significance:
m0_low <- glmmTMB(dist ~ 1 + (1|B/P),family="nbinom1",data=d_low)
anova(m3_low,m0_low,test="Chisq")
# Data: d_low
# Models:
#   m0_low: dist ~ 1 + (1 | B/P), zi=~0, disp=~1
# m3_low: dist ~ fl + (1 | B/P), zi=~0, disp=~1
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# m0_low  4 6226.8 6246.6 -3109.4   6218.8                             
# m3_low  5 6215.8 6240.7 -3102.9   6205.8 12.948      1  0.0003203 ***

# model outcome:
summary(m3_low)
# Family: nbinom1  ( log )
# Formula:          dist ~ f + (1 | B/P)
# Data: d_low
# 
# AIC      BIC   logLik deviance df.resid 
# 6215.8   6240.7  -3102.9   6205.8     1053 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# P:B    (Intercept) 0.36149  0.60124 
# B      (Intercept) 0.00583  0.07635 
# Number of obs: 1058, groups:  P:B, 302; B, 7
# 
# Overdispersion parameter for nbinom1 family (): 3.66 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.83942    0.05951  30.910  < 2e-16 ***
#   fLast       -0.19766    0.05484  -3.604 0.000313 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Anova(m3_low)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: dist
# Chisq Df Pr(>Chisq)    
# f 12.991  1  0.0003131 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

eff_m3_low <- allEffects(m3_low)
plot(eff_m3_low,type="response",ylab="Distance [cm]")

# Factor-level comparison
means_m3_low <- lsmeans(m3_low, pairwise~ f, adjust="bonferroni",type="response")
means_m3_low
# $lsmeans
# f     response    SE   df lower.CL upper.CL
# First     6.29 0.374 1053     5.60     7.07
# Last      5.16 0.343 1053     4.53     5.88
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast     ratio     SE   df t.ratio p.value
# First / Last  1.22 0.0668 1053 3.604   0.0003 
# 
# Tests are performed on the log scale 

# Graph for manuscript:
M <- means_m3_low$lsmeans
md <- as.data.frame(M)

t <- seq(0, 180, by = 1) * pi / 180
r <- .5
x <- r * cos(t)
y <- r*1 * sin(t)
y[20:162] <- y[20] # Flattens the arc

arc.df <- data.frame(Group = x, Value = y)

pdf("210309_DistanceToInnerEdge_Flight_Model.pdf",width = 3.5,height=3.5) 
ggplot(NULL)+
  geom_boxplot(data=d_low, aes(x=f, y=dist),color="grey",width=0.3,fill=c("palegreen","palegreen4"))+
  geom_point(data=md,aes(x=f, y=response), size=1.5, color="black") + 
  geom_errorbar(data=md, aes(ymax=upper.CL, ymin=lower.CL,x=f), width = .1)+
  theme_classic()+
  geom_line(data = arc.df, aes(Group+1.5, Value+33.4), lty = 2) +
  geom_text(data=arc.df, aes(x = 1.5, y = 33.9), label = "***") +
  labs(x="Flight",y="Distance to inner edge [cm]")
dev.off()




#########################################################
## Effects of spatial learning on flight path planning ##
#########################################################

# libraries
library(lme4)
library(DHARMa)
library(sjPlot)
library(car)
library(effects)
library(emmeans)
library(ggplot2)
library(dplyr)
library(pbkrtest)

# data
d_max <- read.table("210324_data_cor_max.txt",header=T,na.string="NA")
str(d_max)
d_max$B <- as.factor(d_max$bat)
d_max$C <- as.factor(d_max$cond)
d_max$fl <- as.factor(d_max$flight)
d_max$t <- d_max$tau/1000

hist(d_max$tau)

# Model
m1 <- lmer(t ~ fl*C + (1|B),data=d_max) 

# Quality check
simulationOutput_m1 <- simulateResiduals(fittedModel = m1, plot = T,n=1000)
# everything ok
testResiduals(simulationOutput_m1)
plotResiduals(simulationOutput_m1, d_max$fl)
plotResiduals(simulationOutput_m1, d_max$C)

# random effects:
plot_model(m1,type="re")

# overall model significance
m0 <- lmer(t ~ 1 + (1|B),data=d_max)
PB <- PBmodcomp(m1, m0, nsim = 1000)
# Bootstrap test; time: 26.64 sec; samples: 1000; extremes: 6;
# large : t ~ fl * C + (1 | B)
# t ~ 1 + (1 | B)
# stat df  p.value   
# LRT    15.992  3 0.001138 **
# PBtest 15.992    0.006993 **

# model outcome
summary(m1)
# Linear mixed model fit by REML ['lmerMod']
# Formula: t ~ fl * C + (1 | B)
# Data: d_max
# 
# REML criterion at convergence: -19.1
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.3766 -0.5955 -0.1705  0.7364  1.6228 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# B        (Intercept) 0.003861 0.06213 
# Residual             0.016182 0.12721 
# Number of obs: 28, groups:  B, 7
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    0.37714    0.05351   7.048
# flLast        -0.14429    0.06800  -2.122
# CPulse        -0.19143    0.06800  -2.815
# flLast:CPulse  0.03286    0.09616   0.342
# 
# Correlation of Fixed Effects:
#   (Intr) flLast CPulse
# flLast      -0.635              
# CPulse      -0.635  0.500       
# flLast:CPls  0.449 -0.707 -0.707

Anova(m1)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: t
# Chisq Df Pr(>Chisq)    
# fl    7.0715  1  0.0078319 ** 
#   C    13.2476  1  0.0002729 ***
#   fl:C  0.1168  1  0.7325851    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

eff_m1 <- allEffects(m1)
plot(eff_m1,type="response",ylab="Tau [sec]",x.var="C")

# facor-level comparison
means_m1 <- lsmeans(m1, pairwise~"C*fl", adjust="bonferroni",type="response")
means_m1
# $lsmeans
# C     fl    lsmean     SE   df lower.CL upper.CL
# echo  first 0.3771 0.0535 21.6   0.2661    0.488
# pulse first 0.1857 0.0535 21.6   0.0746    0.297
# echo  last  0.2329 0.0535 21.6   0.1218    0.344
# pulse last  0.0743 0.0535 21.6  -0.0368    0.185
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                 estimate    SE df t.ratio p.value
# echo first - pulse first   0.1914 0.068 18  2.815  0.0687 t
# echo first - echo last     0.1443 0.068 18  2.122  0.2879 
# echo first - pulse last    0.3029 0.068 18  4.454  0.0018 **
# pulse first - echo last   -0.0471 0.068 18 -0.693  1.0000 
# pulse first - pulse last   0.1114 0.068 18  1.639  0.7118 
# echo last - pulse last     0.1586 0.068 18  2.332  0.1891 


# Graph for manuscript
M <- means_m1$lsmeans
md <- as.data.frame(M)
md$cfl <- paste(md$fl,md$c)


library(tidyverse)
library(ggpubr)

stat.test <- tibble::tribble(
  ~fl, ~group1, ~group2,   ~p.adj,
  "First","Echo",     "Pulse", "p= 0.069",
  "Last","Echo",     "Pulse", "p = 0.19",
)

pdf("210325_TauDirection_Model.pdf",width = 4.5,height=4.5) 
p <- ggplot(md, aes(C,lsmean)) +
  geom_point(aes(colour = factor(C)),size=3.5)
p + scale_colour_manual(values = c("blue", "red"))+ 
  geom_errorbar(data=md, aes(ymax=upper.CL, ymin=lower.CL,x=C,colour=factor(C)), width = .1)+
  facet_wrap(~fl)+
  theme_classic()+
  labs(x="",y=expression(paste("Time lag ", tau, " [sec]")))+
  stat_pvalue_manual(stat.test,y.position=0.51,step.increase=0,label="p.adj")+
  theme(legend.position = "none")

dev.off()


