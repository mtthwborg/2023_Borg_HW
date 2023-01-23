
# OI(I) refers to occupational injury or illness

# Percentile values. Will be updated in S2
predper <- c(seq(0,1,0.1),2,2.5,3:33,100/3,34:66,200/3,67:97,97.5,98,seq(99,100,0.1)) # percentiles (%). 100/3 and 200/3 work even without rounding
predper.short <- c(1,2.5,10,25,50,75,90,97.5,99) # c(1,10,90,99) # percentiles to report RRs for



##################################################
### MODEL, INDOOR/OUTDOOR AND DISTRIBUTION
##################################################

# Type of outcome
if (str_detect(outcome.var, paste(c('claims','injur','disease','llness','OIIs'), collapse="|"))) {
  type.outcome <- 'oi' # injury and/or illness (disease or condition)
} else if (any(str_detect(outcome.var, c('non-zero','non-0','no zero','no 0'))) | isTRUE(plus1)) {
  type.outcome <- 'non-0' # costs without 0 values
} else {
  type.outcome <- 'cost' # costs with 0 values
}

# Outcome distribution
if (distribution.choice %in% c('auto','auto2','automatic','default',NULL)) {
  if (type.outcome=='oi') { # if number of claims/injuries/diseases
    if (distribution.choice=='auto2') {distribution <- 'poisson'}
    else {distribution <- 'quasipoisson'}
  } else if (type.outcome=='cost') { # if cost
    distribution <- tweedie(var.power=1.7, link.power=0) # parameters restimated later; initial choice makes no difference
    tweedie.profile. <- var.power <- list() # if tweedie distribution, create tweedie.profile
  } else {
    distribution <- NULL
    print('NO DISTRIBUTION ASSIGNED')
  }
} else {
  distribution <- distribution.choice
}

# Model type (GLM, GAM, GNM)
if(model.type %in% c('a','gam','GAM')) { # if gam, assign gam
  modeltype <- 'gam'
} else if(model.type=='auto') { # if auto
  if(type.outcome=='oi') { # if number of claims/injuries/diseases
    modeltype <- 'glm' # use glm
  } else {
    modeltype <- 'gam' # use gam
  }
} else {
  modeltype <- model.type # use specified model.type
}

# Select gam or non-gam trend as appropriate
if(modeltype %in% c('a','gam','GAM') ) {  
  trnd <- trend.gam
} else { # glm (Date or date)
  trnd <- trend
}

# Add trnd to equations
mformula <- paste0(m.formula, trnd)



##################################################
### Yaxis for graphs as percentages. Will be adjusted to RR in code as required
##################################################

if (rr.yaxis %in% c('auto','automatic','default',NULL)) { # affect individual model plots
  if(type.outcome=='oi') { # if number of claims/injuries/diseases
    rryaxis <- c(-20,60)
    rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
  } else {
    rryaxis <- c(-20,60)
    rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
  }
} else {
  rryaxis <- rr.yaxis
  rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
}

if (bigvar.yaxis[1] %in% c('auto','automatic','default',NULL)) { # affect bigvar plots
  bigvaryaxis <- rryaxis
} else {
  bigvaryaxis <- bigvar.yaxis
}
# c(-20,30) appropriate for OIs, c(-10,60) for costs. Want to keep axes identical

# overall plots determined per plot



##################################################
### Create folders to store results
##################################################

base <- paste0('Results/',outcome.var, ' ', exposure.var)
dir.create(base) # Warning if exists (doesn't replace)

# folder <- paste0(base, '/')
s1results <- paste0(base,'/Stage 1/')
dir.create(paste0(base,'/Stage 1')) # Warning if exists (doesn't replace)
s2results <- paste0(base,'/Stage 2/')
dir.create(paste0(base,'/Stage 2')) # Warning if exists (doesn't replace)
projresults <- paste0(base,'/Projections/')
dir.create(paste0(base,'/Projections')) # Warning if exists (doesn't replace)


##################################################
### Restrict to study period
##################################################

# Restrict to study period, retaining lagged climate data prior to study period to use for DLNMs
daily.ds <- daily.ds[Date %in% seq(as.Date('2005-07-01')-lmax, as.Date('2018-06-30'), by="day")]
daily.ds[Date < '2005-07-01', (outcome.var):=NA] # ignore outcome in lagged days perior to stud period

# Limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
daily.ds[Month %in% 4:9, (outcome.var):=NA]

# shol
day.before.melbourne.cup <- as.Date(c("2004-11-01","2005-10-31","2006-11-06","2007-11-05","2008-11-03","2009-11-02","2010-11-01","2011-10-31","2012-11-05","2013-11-04","2014-11-03","2015-11-02","2016-10-31","2017-11-06","2018-11-05")) # do manually, as shift seems to miss 31st October
daily.ds[City=='Melbourne' & Date %in%  day.before.melbourne.cup, shol:='Day before Melbourne Cup'] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
daily.ds[Month==12 & Day %in% c(24,25), shol:='Xmas Day or Eve'] # 23rd Dec seems to results in high res for Melbourne and Sydney contrasting with generally negative res in Xmas break


##################################################
### Create objects to store results
##################################################

# Matrices with 1 column: Number of outcomes, formula, exposure at 0, AIC and dispersion (if poisson used)
no.outcome <- no.years <- m.aic <- m.dispersion <- m.r2 <- m.devexp <- matrix(NA, length.ds.city, 1, dimnames=list(ds.city)) # matrix of NAs, with a column to denote total number of outcomes
colnames(m.aic) <- 'AIC'
colnames(m.dispersion) <- 'Dispersion parameter'
colnames(m.r2) <- 'R^2'
colnames(m.devexp) <- 'Deviance explained'

# Exposure. Matrix of NAs, with a column to denote total number of outcomes
exposure <- matrix(NA, length.ds.city, 4 + length(predper.short), dimnames=list(ds.city, c('mean','range','min','max',paste0(predper.short,'%'))))

# Percentile exposure per by variables
exposure.by <- matrix(NA, length.ds.city, length(predper), dimnames=list(ds.city, predper))

# EHF. Heatwave if >0, though as a continuous variable, they are essentially identical
exposure_ehf <- matrix(NA, length.ds.city, 4, dimnames=list(ds.city, c('Heatwave','Severe','Extreme','Exorbitant'))) # matrix of NAs, with a column to denote total number of outcomes
ehf_threshold <- matrix(NA, length.ds.city, 4, dimnames=list(ds.city, c('Heatwave','Severe','Extreme','Exorbitant'))) # matrix of NAs, with a column to denote total number of outcomes

# Coefficients for overall cumulative summary
if(espline=='bs') {
  espline.length <- length(eknots) + 3 # not sure if correct, based on testing. If wrong, error: Error in coef[i, ] <- coef(red[[i]]) : number of items to replace is not a multiple of replacement length
  coef <- matrix(NA, length.ds.city, espline.length, dimnames=list(ds.city, paste0('b',rep(1:espline.length)))) # matrix of NAs, with rows per by.vars
} else {
  coef <- matrix(NA, length.ds.city, edf, dimnames=list(ds.city, paste0('b',rep(1:edf)))) # matrix of NAs, with rows per by.vars
}

# Covariance for overall cumulative summary
vcov <- vector("list", length.ds.city)
names(vcov) <- ds.city

# Tweedie max shape parameter values (if Tweedie distribution used)
m.tweedie.shape <- mean_rhs <- matrix(NA, length.ds.city, 1, dimnames=list(ds.city)) # matrix of NAs, with a column to denote optimised shape parameter

# List objects
rhs <- temps <- outcomes <- model <- model.omit <- model2 <- m.coef <- res <- model.checks <- check <- model.checks.month <- model.checks.week <- check.res.coef <- check.res.coef.sig <- exposure.rr.s1 <- red <- list() #  Model, residuals, residual length + plots and earity



################################################################################
# FIRST-STAGE ANALYSIS: MODEL FOR EACH BY-VARIABLE COMBINATION, REDUCE AND SAVE, NO POOLING, include dummy dlnm rows
################################################################################

gc()
time <- proc.time()[3]
for(i in ds.city) {
  print(paste('Stage 1 model:',i))
  .ds <- daily.ds[City==i] # dataset for each unique value
  .d <- .ds[-(1:lmax)] # dataset for each unique value without lagged values prior to study period
  .name <- unique(.ds$stratum) # names with all of a, b and c together

  # Dependent and independent variables using .ds
  .outcome <- outcomes[[i]] <- .d[,get(outcome.var)]
  no.outcome[i,] <- sum(outcomes[[i]] != 0) # number of non-zero (and non-missing) outcomes
  temps[[i]] <- .d[,get(exposure.var)] # exposure. Include all temperature values used for dlnm including lagged values (early ones are used less, but so are day 0 values towards end of study)
  # if (no.outcome[i,] < 100) {next} # skip iteration if number of outcomes is less than 100, which can lead to non-convergence. Still results in model output
  
  # Centre (reference value) on mean for crosspred
  .cen <- mean(temps[[i]], na.rm=T) 
  
  # Calculate crossbasis. Add group if restrict analysis to summer
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag)
  
  # Mean exposure and range, for meta-predictors. Note, this is only dependent on location and indoor/outdoor
  exposure[i,1] <- mean(temps[[i]], na.rm=T) # mean
  exposure[i,2] <- diff(range(temps[[i]], na.rm=T)) # range
  exposure[i,3] <- min(temps[[i]], na.rm=T) # min
  exposure[i,4] <- max(temps[[i]], na.rm=T) # max
  exposure[i,5:(length(predper.short)+4)] <- quantile(temps[[i]], predper.short/100, na.rm=T)
  exposure.by[i,] <- quantile(temps[[i]], predper/100,na.rm=T) # all percentiles and related average exposure
  exposure_ehf[i,1] <- 0
  exposure_ehf[i,2] <- quantile(temps[[i]][temps[[i]]>=0],0.85)
  exposure_ehf[i,3] <- quantile(temps[[i]][temps[[i]]>=0],0.85)*2
  exposure_ehf[i,4] <- quantile(temps[[i]][temps[[i]]>=0],0.85)*3
  ehf_threshold[i,1] <- ecdf(temps[[i]])(0) # percentile corresponding with 0. Needed for EHF. If does not exist, it becomes 1
  ehf_threshold[i,2] <- ecdf(temps[[i]])(exposure_ehf[i,2]) # percentile corresponding with severe ehf heatwave
  ehf_threshold[i,3] <- ecdf(temps[[i]])(exposure_ehf[i,3]) # percentile corresponding with severe ehf heatwave * 2
  ehf_threshold[i,4] <- ecdf(temps[[i]])(exposure_ehf[i,4]) # percentile corresponding with severe ehf heatwave * 3
  rhs[[i]] <- .d[,get('Average relative humidity')] # relative humidity
  mean_rhs[i,] <- mean(.d[,get('Average relative humidity')], na.rm=T) # relative humidity
  
  # Number of years in model
  no.years[i,] <- max(.d[,FYear], na.rm=T) - min(.d[,FYear], na.rm=T) + 1
  .no.years <- no.years[i,]
  
  # Generate model
  if(modeltype %in% c('a','gam','GAM')) {
    if (distribution[1]=='Tweedie') {
      model[[i]] <- gam(as.formula(mformula), family=tw(), .ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
      m.tweedie.shape[i,] <- as.numeric(str_remove_all(model[[i]]$family[["family"]], "[Twediep=()]")) # save shape parameter
    } else {
      model[[i]] <- gam(as.formula(mformula), family=distribution, data=.ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
    }
   
    # GAM SE, t-values, P-values and R^2
    .m.se <- summary(model[[i]])[["se"]]
    .m.t <- summary(model[[i]])[["p.t"]] # values not produced for gam smoothers
    .m.p <- summary(model[[i]])[["p.pv"]] # values not produced for gam smoothers, that is summary(model[[i]])$s.table[,4]
    m.r2[i,] <- summary(model[[i]])[["r.sq"]] # R2
  } else {
    if (distribution[1]=='Tweedie') { 
      var.power[[i]] <- tweedie.profile(formula=mformula, data=.ds, p.vec=var.power.values, method='series', do.ci=F, do.smooth=tp.smooth) # estimate optimal shape parameter. Turn off CI for speed
      m.tweedie.shape[i,] <- var.power[[i]]$p.max # save shape parameter
      model[[i]] <- glm(as.formula(mformula), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # "na.exclude". na.omit works better for Tweedie e.g. AIC calculation
      # model.omit[[i]] <- glm(as.formula(mformula), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # NAs will cause tweedie AIC calculation to fail 
    } else {
      model[[i]] <- glm(as.formula(mformula), family=distribution, .d, na.action="na.exclude", control=list(maxit=max.iter))
    }
    
    # GLM SE, t-values, P-values and R^2
    .glm.coef <- summary(model[[i]])$coefficients
    .m.se <- .glm.coef[,2]
    .m.t <- .glm.coef[,3]
    .m.p <- .glm.coef[,4]
    m.r2[i,] <- 1 - model[[i]]$deviance/model[[i]]$null.deviance # ?Nagelekre R&2
  }
  
  # Save coefficients
  .m.coef <- coef(model[[i]])
  .m.row <- names(.m.coef)
  .m.lci <- .m.coef - qnorm(0.975)*.m.se
  .m.uci <- .m.coef + qnorm(0.975)*.m.se
  suppressWarnings(m.coef[[i]] <- data.table(location=i, Parameter=.m.row, Coefficient=.m.coef, SE=.m.se, LCI=.m.lci, UCI=.m.uci, RR=exp(.m.coef), 'RR LCI'=exp(.m.lci), 'RR UCI'=exp(.m.uci), 't-value'=.m.t, 'p-value'=round(.m.p,8))) # combine coefficients into a data.table
  m.coef[[i]][!1, location:=NA] # only keep location for 1st row
  m.coef[[i]][`p-value` <0.1, sig:='.']
  m.coef[[i]][`p-value` <0.05, sig:='*']
  m.coef[[i]][`p-value` <0.01, sig:='**']
  m.coef[[i]][`p-value` <0.001, sig:='***']
  if(modeltype %in% c('a','gam','GAM')) {
    m.coef[[i]][str_detect(Parameter,'\\('), ':='(`t-value`=NA, `p-value`=NA)] # t-values and p-values not calculated for gam smoothers, so set NA instead of repeat for them
  }
  
  if(distribution[1]=='quasipoisson') { #[1] is to avoid warnings from using multiple elements
    if(modeltype %in% c('a','gam','GAM')) {
      m.aic[i,] <- fqaic.fn(model[[i]], gam=T) # if gam, use gam option
    } else {
      m.aic[i,] <- fqaic.fn(model[[i]])
    } 
  } else if(distribution[1]=='Tweedie') {
    if(modeltype %in% c('a','gam','GAM')) {
      m.aic[i,] <- AIC(model[[i]]) # conditional AIC, which is better for GAM # model[[i]]$aic is uncorrected AIC
    } else {
      m.aic[i,] <- AICtweedie(model[[i]])
    }
  } else {
    m.aic[i,] <- AIC(model[[i]]) # model[[i]]$aic is identical
  }
  m.devexp[i,] <- (model[[i]]$null.deviance - model[[i]]$deviance) / model[[i]]$null.deviance # deviance explained
  m.dispersion[i,] <- sum(residuals(model[[i]], type="pearson")^2, na.rm=T)/df.residual(model[[i]]) # summary(model[[i]])[["dispersion"]]
  
  # Sum (accumulate) effects of all lags in order to eliminate one dimension of the association
  # Predicted effects: extract parameters from model corresponding to .cb variables through functions coef and vcov
  .pred1 <- crosspred(basis=.cb, model=model[[i]], cen=.cen, model.link='log') # must set log for tweedie package 
 
  # Overall cumulative summary for main model
  red[[i]] <- crossreduce(.cb, model[[i]], cen=.cen, model.link='log') # reduce exposure-lag-response association to overall cumulative exposure-response association, using mean
  coef[i,] <- coef(red[[i]])
  vcov[[i]] <- vcov(red[[i]])
  
  exposure.rr.s1[[i]] <- data.frame(matrix(nrow=length(red[[i]]$predvar), ncol=4))
  colnames(exposure.rr.s1[[i]]) <- c("temp", "RRfit", "RRlow", "RRhigh")
  exposure.rr.s1[[i]]$temp <- red[[i]]$predvar # = .pred1$predvar
  
  if(any(str_detect(names(red[[i]]), 'RRfit'))) { # default dlnm behaviour for log-link
    exposure.rr.s1[[i]]$RRfit  <- red[[i]]$RRfit # already %, do not need to convert like in supp. Replaced RRfit/low/high with fit/low/high. Crosspred does output se, ?can manully calculate robust SE
    exposure.rr.s1[[i]]$RRlow <- red[[i]]$RRlow
    exposure.rr.s1[[i]]$RRhigh <- red[[i]]$RRhigh
  }  else if(any(str_detect(names(red[[i]]), 'fit'))) { # if Tweedie glm, which outputs a log-link value not recognised by dlnm
    exposure.rr.s1[[i]]$RRfit  <- exp(red[[i]]$fit) # already %, do not need to convert like in supp. Replaced RRfit/low/high with fit/low/high. Crosspred does output se, ?can manully calculate robust SE
    exposure.rr.s1[[i]]$RRlow <- exp(red[[i]]$low)
    exposure.rr.s1[[i]]$RRhigh <- exp(red[[i]]$high)
  } 
  
  # Overall cumulative exposure-response relationship plot. We're not using these plots for results, but for checking
  # Convert RR into %
  .exposure.rr.s1.rrfit <- (exposure.rr.s1[[i]]$RRfit-1)*100 
  .exposure.rr.s1.rrlow <- (exposure.rr.s1[[i]]$RRlow-1)*100
  .exposure.rr.s1.rrhigh <- (exposure.rr.s1[[i]]$RRhigh-1)*100
  
  .ind1 <- .pred1$predvar<=.cen # values below centering value (blue)
  .ind2 <- .pred1$predvar>=.cen # values above centering value (red)
  
  # Overall exposure response
  plot(exposure.rr.s1[[i]]$temp, .exposure.rr.s1.rrfit, type="n", ylim=c(rryaxis[1], rryaxis[2]), lwd=2, col="white",
       main=i, ylab="Percent change (%)", xlab=paste(exposure.var,'(Â°C)'), 
       cex.main=cexmain, cex.lab=cexlab, cex.axis=cexaxis, lab=c(6,5,7))
  .erplot %<a-% { # save main plot additions.
    polygon(c(exposure.rr.s1[[i]]$temp,rev(exposure.rr.s1[[i]]$temp)),c(.exposure.rr.s1.rrlow,rev(.exposure.rr.s1.rrhigh)), col="grey89", border=F) # 95% CI envelope
    lines(exposure.rr.s1[[i]]$temp[.ind1],.exposure.rr.s1.rrfit[.ind1],col=4,lwd=2); # cold, left of cen
    lines(exposure.rr.s1[[i]]$temp[.ind2],.exposure.rr.s1.rrfit[.ind2],col=2,lwd=2); # hot, right of cen
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure[i,c("1%","99%")]),lty=3) # 1st and 99th percentiles
  }
  .erplot # insert main plot additions
  .lines # insert lines
} 
proc.time()[3]-time
# traceback()



################################################################################
# Combine and review results of Stage 1 Analysis
################################################################################

# Combine number of outcomes, DLNM coefficients, AIC and dispersion parameter, including totals for outcomes and AIC
m1.results <- rbind(cbind(no.outcome, coef, m.aic,  m.dispersion), c(sum(no.outcome),rep(NA, ncol(coef)),sum(m.aic), NA))
colnames(m1.results) <- c('Number of outcomes',colnames(coef),'AIC','Dispersion parameter')
write.csv(m1.results, file=paste0(s1results,'Reduced coef and model fit.csv'), na='', row.names=T) # output AIC per model and its sum

# Output model coefficients and RRs
write.csv(do.call(rbind, m.coef), file=paste0(s1results,'Coefficients.csv'), na='') # Model coefficients, SE, 95% CI and t-tests
# do.call(rbind, m.coef)[str_detect(parameter,'day1')] 

# Tweedie shape parameters
if (distribution[1]=='Tweedie') { 
  write.csv(m.tweedie.shape, file=paste0(s1results,'Tweedie shape parameters.csv'), na='', row.names=T)
} 



######################### END ############################
######################### END ############################
