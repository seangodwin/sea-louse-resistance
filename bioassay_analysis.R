## INPUTS ------------------
bioassay.datafile <- "bioassay_data.csv"
treat.datafile <- "treatment_data.csv"
  

## SET-UP ------------------
library(lme4)       # GLMs
library(tidyverse)  # dplyr etc
library(AICcmodavg) # AICc table
library(expandFunctions)# warning messages
library(MASS)       # multivariate normal for CIs

# Read in bioassay data
setwd(dir.data)
data <- read.csv(bioassay.datafile, header=T, stringsAsFactors=F)
treat <- read.csv(treat.datafile, header=T, stringsAsFactors=F)


## MODELS FIT TO BIOASSAY SURVIVAL DATA ------------------
# Empty dataframes
assay.ci <- NULL
ec50 <- NULL

# Scale concentration values
data$conc.sc <- as.numeric(scale(data$conc, center=T, scale=T))

# Fit glmm to data from each bioassay
reset.warnings()
warn <- 0
for(i in 1:length(unique(data$assay.id))) {
  
  # Model
  mod <- glmer(live ~ conc.sc*sex + (1|petri.id), data=data[data$assay.id == i,], family=binomial(link="logit"),
                 control=glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5)))
  
  
    # If there's a warning, fit model again with updated starting conditions
    # This is an ugly way to do it, but newest version of R doesn't store last.warning in workspace
    warn <- c(warn, length(warnings()))
    if(warn[i+1] > warn[i]) {
      ss <- getME(mod,c("theta","fixef"))
      mod <- update(mod,start=ss,control=glmerControl(optCtrl=list(maxfun=2e5)))
    }
  
  # Calculate conventional CIs for the model 
  # Function taken from Ben Bolker: https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html#prediction
  glmer.ci <- function(model,newdata=NULL,alpha=0.05) {
    ## baseline prediction, on the linear predictor (logit) scale:
    pred0 <- predict(model,re.form=NA,newdata=newdata)
    ## fixed-effects model matrix for new data
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
    beta <- fixef(model) ## fixed-effects coefficients
    V <- vcov(model)     ## variance-covariance matrix of beta
    pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
    ## inverse-link function
    linkinv <- family(model)$linkinv
    ## construct 95% Normal CIs on the link scale and
    ##  transform back to the response (probability) scale:
    crit <- -qnorm(alpha/2)
    linkinv(cbind(conf.low=pred0-crit*pred.se,
                  conf.high=pred0+crit*pred.se))
  }
  
  # Prediction dataframe
  conc.seq <- seq(0,1200,1)
  conc.sc.seq <- (conc.seq - mean(data$conc)) / sd(data$conc)
  newdat <- expand.grid(conc.sc = conc.sc.seq,
                        sex = factor(c("male","female")))
  
  # Generate predictions
  assay.ci <- as.data.frame(rbind(assay.ci, cbind(assay.id = i,
                                             date = data$date[data$assay.id == i][1],
                                             farm = data$farm[data$assay.id == i][1],
                                             newdat, 
                                             mean = predict(mod, newdat, type="response", re.form=NA),
                                             lwr = glmer.ci(mod, newdat)[,1],
                                             upr = glmer.ci(mod, newdat)[,2])))
  
  # Store EC50s
  sex <- c("male", "female")
  for(m in 1:length(sex)) {
    ec50 <- as.data.frame(rbind(ec50, data.frame(assay.id = i,
                                   date = data$date[data$assay.id == i][1],
                                   farm = data$farm[data$assay.id == i][1],
                                   mean = rev(conc.seq)[findInterval(0.5, sort(assay.ci$mean[assay.ci$assay.id==i & assay.ci$sex==sex[m]]))],
                                   lwr = rev(conc.seq)[findInterval(0.5, sort(assay.ci$lwr[assay.ci$assay.id==i & assay.ci$sex==sex[m]]))],
                                   upr = rev(conc.seq)[findInterval(0.5, sort(assay.ci$upr[assay.ci$assay.id==i & assay.ci$sex==sex[m]]))],
                                   sex = sex[m],
                                   prev.emb = data$prev.emb[data$assay.id == i][1])))
  }
} 
# Singular fit explained in methods section of paper
# Warnings fixed in function above (re-fit with updated starting conditions, no warnings)

# Fix date and add days since 2010-01-01
ec50$date <- as.Date(ec50$date)
ec50$day <- as.numeric(ec50$date - as.Date("2010-01-01"))

# Write EC50 data
write.csv(ec50, "ec50_data.csv", quote=F, row.names=F)


## MODELS FIT TO EC50s ------------------
# Models
mods.ec50 <- list()
mods.ec50[[1]] <- glm(mean ~ sex + prev.emb, data=ec50, family=Gamma(link="log"))
mods.ec50[[2]] <- glm(mean ~ sex + day + prev.emb, data=ec50, family=Gamma(link="log"))
mods.ec50[[3]] <- glm(mean ~ sex * day  + prev.emb, data=ec50, family=Gamma(link="log"))
mods.ec50[[4]] <- glm(mean ~ sex + day + I(day^2) + prev.emb, data=ec50, family=Gamma(link="log"))
mods.ec50[[5]] <- glm(mean ~ sex * (day + I(day^2)) + prev.emb, data=ec50, family=Gamma(link="log"))

# AICc table
aictab.ec50 <- as.data.frame(aictab(cand.set = mods.ec50, 
                               modnames = c("treated + sex", "treated + sex + day", "treated + sex * day",
                                            "treated + sex + day^2", "treated + sex * day^2"), 
                               sort = F, second.ord=T))#[,c(1,4,6)]

# Top model
top.mod.ec50 <- mods.ec50[[which(aictab.ec50$Delta_AICc == min(aictab.ec50$Delta_AICc))]]

# Prediction dataframe
day.seq <- seq(0,4500,1)
newdat.ec50 <- expand.grid(day = day.seq,
                      sex = factor(c("male","female")),
                      prev.emb = mean(ec50$prev.emb))

# Link function and alpha
linkinv <- family(top.mod.ec50)$linkinv
alpha <- 0.05

# Generate predictions
pred0 <- predict(top.mod.ec50, newdat.ec50, type="link", se.fit=T)
pred.ec50 <- cbind(newdat.ec50,
                   date = as.Date("2010-01-01") + newdat.ec50$day,
                   linkinv(cbind(mean = pred0$fit,
                                 lwr = pred0$fit - -qnorm(alpha/2) * pred0$se.fit,
                                 upr = pred0$fit + -qnorm(alpha/2) * pred0$se.fit)))


## MODELS FIT TO RELATIVE COUNTS ------------------
# Remove treatments that are not slice, that have NA for rel.count, or that were treated again within 3 months
treat.slice <- treat[treat$type == "slice" & is.na(treat$rel.count)==F & treat$treated.again==0,]
treat.slice$non.zero <- ifelse(treat.slice$rel.count > 0, 1, 0)
treat.slice$date <- as.Date(treat.slice$date)
  
# Binomial
mods.treat.binom <- list()
mods.treat.binom[[1]] <- glm(non.zero ~ prev.emb, data=treat.slice, family=binomial(link="logit"))
mods.treat.binom[[2]] <- glm(non.zero ~ prev.emb + day, data=treat.slice, family=binomial(link="logit"))
mods.treat.binom[[3]] <- glm(non.zero ~ prev.emb + day + I(day^2), data=treat.slice, family=binomial(link="logit"))

# Gamma models
mods.treat.gamma <- list()
mods.treat.gamma[[1]] <- glm(rel.count ~ prev.emb, data=treat.slice[treat.slice$rel.count != 0,], family=Gamma(link="log"))
mods.treat.gamma[[2]] <- glm(rel.count ~ prev.emb + day, data=treat.slice[treat.slice$rel.count != 0,], family=Gamma(link="log"))
mods.treat.gamma[[3]] <- glm(rel.count ~ prev.emb + day + I(day^2), data=treat.slice[treat.slice$rel.count != 0,], family=Gamma(link="log"))

# Calculate AICc values for each hurdle model
aicc.fun <- function(mbinom, mgamma) {
  k <- length(coef(mbinom)) + length(coef(mgamma))
  L <- as.numeric(logLik(mbinom) + logLik(mgamma))
  n <- length(mbinom$residuals)
  aic <- -2*L + 2*k
  aicc <- aic + (2*k^2 + 2*k)/(n - k - 1)
  return(aicc)
}
aicc.hurdle <- mapply(aicc.fun, mods.treat.binom, mods.treat.gamma)
delta.aicc.hurdle <- aicc.hurdle - min(aicc.hurdle)
weights.hurdle <- exp(-0.5 * delta.aicc.hurdle) / sum(exp(-0.5 * delta.aicc.hurdle))
aictab.hurdle <- data.frame(modnames = c("treated", "treated + day", "treated + day^2"),
                            delta = delta.aicc.hurdle,
                            weights = weights.hurdle)

# Top models
top.mod.treat.binom <- mods.treat.binom[[which(aictab.hurdle$delta == min(aictab.hurdle$delta))]]
top.mod.treat.gamma <- mods.treat.gamma[[which(aictab.hurdle$delta == min(aictab.hurdle$delta))]]



## HURDLE MODEL CIS ---------------------
# Prediction dataframe (kind of silly but for consistency)
day.seq.treat <- seq(0,4500,5)
newdat.treat <- expand.grid(prev.emb = mean(treat$prev.emb, na.rm=T),
                            day = day.seq)

# Simulate data from the models
nsim <- 10000
sims.gamma <- simulate(top.mod.treat.gamma, nsim)
sims.binom <- simulate(top.mod.treat.binom, nsim)

# Parametric bootstrap
refit <- function(x1, x2){
  data.binom <- top.mod.treat.binom$data
  data.binom$non.zero <- x1
  newmod.binom <- update(top.mod.treat.binom, data = data.binom)
  
  data.gamma <- top.mod.treat.gamma$data
  data.gamma$rel.count <- x2
  newmod.gamma <- update(top.mod.treat.gamma, data = data.gamma)
  
  predict(newmod.binom, newdat.treat, type="response") * predict(newmod.gamma, newdat.treat, type="response")
}
boot <- mapply(refit, sims.binom, sims.gamma) # warnings checked
# it's not perfect but best we can do, and just used for visual

# Generate predictions
pred.treat <- data.frame(date = as.Date("2010-01-01") + newdat.treat$day,
                         mean = predict(top.mod.treat.binom, newdat.treat, type="response") *
                                predict(top.mod.treat.gamma, newdat.treat, type="response"),
                         lwr=apply(boot, 1, function(x) {quantile(x, probs=c(0.025,0.975))})[1,], 
                         upr=apply(boot, 1, function(x) {quantile(x, probs=c(0.025,0.975))})[2,])

