# load packages
library(did) # devtools::install_github("bcallaway11/did")
library(BMisc) # devtools::install_github("bcallaway11/BMisc")
library(twfeweights) # devtools::install_github("bcallaway11/twfeweights")
library(HonestDid) # devtools::install_github("asheshrambachan/HonestDID")
library(fixest) # install.packages("fixest")
library(modelsummary) # install.packages("model.summary")
library(ggplot2) # install.packages("ggplot2")

# load data
load(url("https://github.com/bcallaway11/did_chapter/raw/master/mw_data_ch2.RData"))


# setup data
# drops NE region and a couple of small groups
mw_data_ch2 <- subset(mw_data_ch2, (G %in% c(2004,2006,2007,0)) & (region != "1"))
head(mw_data_ch2[,c("id","year","G","lemp","lpop","region")])
# drop 2007 as these are right before fed. minimum wage change
data2 <- subset(mw_data_ch2, G!=2007 & year >= 2003)
# keep 2007 => larger sample size
data3 <- subset(mw_data_ch2, year >= 2003)


# twfe regression
twfe_res2 <- feols(lemp ~ post | id + year,
                  data=data2,
                  cluster="id")

modelsummary(list(twfe_res2), gof_omit=".*")



# Callaway and Sant'Anna ATT(g,t)
attgt <- did::att_gt(yname="lemp",
                     idname="id",
                     gname="G",
                     tname="year",
                     data=data2,
                     control_group="nevertreated",
                     base_period="universal")
tidy(attgt)[,1:5] # print results, drop some extra columns


# plot ATT(g,t)'s
ggdid(attgt)


# aggregate into att^O
attO <- aggte(attgt, type="group", na.rm=TRUE)
summary(attO)


# compute twfe weights
tw <- twfeweights::twfe_weights(attgt)
tw <- tw[tw$G != 0,]
tw$post <- as.factor(1*(tw$TP >= tw$G))
twfe_est <- sum(tw$wTWFEgt*tw$attgt)

# plot weights / att(g,t)
ggplot(data=tw,
       mapping=aes(x=wTWFEgt, y=attgt, color=post)) +
  geom_hline(yintercept=0, linewidth=1.5) +
  geom_vline(xintercept=0, linewidth=1.5) + 
  geom_point(size=6) +
  theme_bw() +
  ylim(c(-.15,.05)) + xlim(c(-.4,.7))


# calculate att^O weights
wO <- attO_weights(attgt)
wO <- wO[wO$G != 0,]
attO_est = sum(wO$wOgt*wO$attgt)
wO$post <- as.factor(1*(wO$TP >= wO$G))

# plot att^O weights / attgt
ggplot(data=wO,
       mapping=aes(x=wOgt, y=attgt, color=post)) +
  geom_hline(yintercept=0, linewidth=1.5) +
  geom_vline(xintercept=0, linewidth=1.5) + 
  geom_point(shape=18, size=8) +
  theme_bw() +
  ylim(c(-.15,.05)) + xlim(c(-.4,.7))


# plot a comparison of the weights
plot_df <- cbind.data.frame(tw, wOgt=wO$wOgt)
plot_df <- plot_df[plot_df$post==1,]
plot_df$g.t <- as.factor(paste0(plot_df$G,",",plot_df$TP))

ggplot(plot_df, aes(x=wTWFEgt, y=attgt, color=g.t)) +
  geom_point(size=6) +
  theme_bw() +
  ylim(c(-.15,.05)) + xlim(c(-.4,.7)) + 
  geom_point(aes(x=wOgt), shape=18, size=8) +
  geom_hline(yintercept=0, linewidth=1.5) +
  geom_vline(xintercept=0, linewidth=1.5) + 
  xlab("weight")


# comparison of weights / bias from twfe
twfe_post <- sum(tw$wTWFEgt[tw$post==1] * tw$attgt[tw$post==1])
twfe_post

# pre-treatment contamination/bias
pre_bias <- sum(tw$wTWFEgt[tw$post==0] * tw$attgt[tw$post==0])
pre_bias

twfe_bias <- twfe_est - attO_est
pre_bias/twfe_bias # bias from pre-treatment PTA violations
(twfe_post-attO_est)/twfe_bias # bias from TWFE weights instead of ATT^O weights



# twfe with covariates
twfe_x <- feols(lemp ~ post | id + region^year,
                data=data2)
modelsummary(twfe_x, gof_omit=".*")


# att(g,t) with covariates
cs_x <- att_gt(yname="lemp",
               tname="year",
               idname="id",
               gname="G",
               xformla=~region,
               control_group="nevertreated",
               base_period="universal",
               data=data2)
cs_x_res <- aggte(cs_x, type="group")
summary(cs_x_res)


# plot event study
ggdid(aggte(attgt,type="dynamic",cband=FALSE))


# honest did code
library(HonestDiD)
source("honest_did.R") # download from https://bcallaway11.github.io/presentations/NEXT-D/honest_did.R
cs_es <- aggte(attgt, type="dynamic")
hd_cs <- honest_did(es = cs_es, 
                          e = 0,
                          type="relative_magnitude")
createSensitivityPlot_relativeMagnitudes(hd_cs$robust_ci,
                                         hd_cs$orig_ci)

