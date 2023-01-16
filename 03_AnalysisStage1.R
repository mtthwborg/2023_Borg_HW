
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
### Create objects to store results
##################################################

# Matrices with 1 column: Number of outcomes, formula, exposure at 0, AIC and dispersion (if poisson used)
no.outcome <- no.years <- m.aic <- m.dispersion <- m.r2 <- m.devexp <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote total number of outcomes
colnames(m.aic) <- 'AIC'
colnames(m.dispersion) <- 'Dispersion parameter'
colnames(m.r2) <- 'R^2'
colnames(m.devexp) <- 'Deviance explained'

# Exposure. Matrix of NAs, with a column to denote total number of outcomes
exposure <- matrix(NA, length.ds.stratum, 4 + length(predper.short), dimnames=list(ds.stratum, c('mean','range','min','max',paste0(predper.short,'%'))))

# Percentile exposure per by variables
exposure.by <- matrix(NA, length.ds.stratum, length(predper), dimnames=list(ds.stratum, predper))

# Coefficients for overall cumulative summary
if(espline=='bs') {
  espline.length <- length(eknots) + 3 # not sure if correct, based on testing. If wrong, error: Error in coef[i, ] <- coef(red[[i]]) : number of items to replace is not a multiple of replacement length
  coef <- matrix(NA, length.ds.stratum, espline.length, dimnames=list(ds.stratum, paste0('b',rep(1:espline.length)))) # matrix of NAs, with rows per by.vars
} else {
  coef <- matrix(NA, length.ds.stratum, edf, dimnames=list(ds.stratum, paste0('b',rep(1:edf)))) # matrix of NAs, with rows per by.vars
}

# Covariance for overall cumulative summary
vcov <- vector("list", length.ds.stratum)
names(vcov) <- ds.stratum

# Tweedie max shape parameter values (if Tweedie distribution used)
m.tweedie.shape <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote optimised shape parameter

# List objects
temps <- outcomes <- model <- model.omit <- model2 <- m.coef <- res <- model.checks <- check <- model.checks.month <- model.checks.week <- check.res.coef <- check.res.coef.sig <- exposure.rr.s1 <- red <- list() #  Model, residuals, residual length + plots and earity



################################################################################
# FIRST-STAGE ANALYSIS: MODEL FOR EACH BY-VARIABLE COMBINATION, REDUCE AND SAVE, NO POOLING, include dummy dlnm rows
################################################################################

gc()
time <- proc.time()[3]
for(i in ds.stratum) {
  print(paste('Stage 1 model:',i))
  .ds <- daily.ds[stratum==i] # dataset for each unique value #
  if(any(is.na(.ds[, get(outcome.var)]))) {
    .d <- .ds[!is.na(get(outcome.var))] # no missing outcome. SHOULD NOT, BUT BE CAUTIOUS, AS IT COULD UPSET WBGT ANALYSIS
    # .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
  }
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  adf <- asdf[stratum==i, get(s.sdf)] # select appropriate sdf, if it's to be used
  
  # Dependent and independent variables using .ds
  .outcome <- .ds[,get(outcome.var)]
  outcomes[[i]] <- .d[,get(outcome.var)]
  no.outcome[i,] <- sum(outcomes[[i]] != 0) # number of days with non-zero (and non-missing) outcomes
  no.oii[i,] <- sum(.ds[,get("Number of OIIs")], na.rm=T) # number of claims (OIIs)
  temps[[i]] <- .ds[,get(exposure.var)] # exposure. Include all temperature values used for dlnm including lagged values (early ones are used less, but so are day 0 values towards end of study)
  hums[[i]] <- .ds[,get('Average specific humidity')] # exposure. Include all temperature values used for dlnm including lagged values (early ones are used less, but so are day 0 values towards end of study)
  
  # Define exposure
  .cen <- 0 # centred on 0
  
  # Calculate crossbasis. Add group if restrict analysis to summer
  li_argvar[[i]] <- list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T), Bound=range(temps[[i]], na.rm=T)) # need Bound for projections, specifcally do.call(onebasis,...)
  .cb <- crossbasis(temps[[i]], argvar=li_argvar[[i]], lag=lmax, arglag=li_arglag)
  .hcb <- crossbasis(hums[[i]], argvar=list(fun='lin'), lag=2, arglag=list(fun='integer'))
  
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
  rhs[[i]] <- .d[,ave_rh] # relative humidity
  mean_rhs[i,] <- mean(.d[,ave_rh], na.rm=T) # relative humidity
  
  # Number of years in model: 13, 12 for TAS
  no.years[i,] <- max(.d[,FYear], na.rm=T) - min(.d[,FYear], na.rm=T) + 1
  .no.years <- no.years[i,]
  
  # Adjust formula as needed
  form[i,] <- mformula # use indoor formula
  if(!is.null(perth.patch) & isTRUE(fyearp) & str_detect(i, 'Perth')) { # if perth.patch is not NULL, insert formula patchin formula. In this case, only if relevant (costs and Perth)
    form[i,] <- paste0(form[i,], perth.patch)
  } 
  
  if(str_detect(i,'Adelaide')) { # no PH in Adelaide on Sun. Remove else perfect multicollinearity results. Does improve fit for other cities, even Hobart which only has 5 days
    form[i,] <- str_remove_all(form[i,], ' \\+ Sun:public.hol') 
    form[i,] <- str_remove_all(form[i,], ' \\+ Sun\\*public.hol')
  }
  if(str_detect(i,'Hobart')) {  # for Oct-Nov, perfect collinearity results
    form[i,] <- str_remove_all(form[i,], ' \\+ Sat:public.hol')
    form[i,] <- str_remove_all(form[i,], ' \\+ Sat\\*public.hol')
  }
  # Generate model
  if(modeltype %in% c('a','gam','GAM')) {
    if (distribution[1]=='Tweedie') {
      model[[i]] <- gam(as.formula(form[i,]), family=tw(), .ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
      m.tweedie.shape[i,] <- as.numeric(str_remove_all(model[[i]]$family[["family"]], "[Twediep=()]")) # save shape parameter
    } else {
      model[[i]] <- gam(as.formula(form[i,]), family=distribution, data=.ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
    }
    if(no.smooth==1) {
      tryCatch({gam.k[i,] <- k.check.MAB(model[[i]], subsample=nrow(.ds)-lmax)}, warning=function(w) print(i)) # check k' and edf. MAB version removes NA residuals, which is required with DLNM. Default subsample is 5000, which exceeds data limit. Lower data limit to be no more than the number of rows in dataset excluding the dlnm lag NA rows
    } else if(no.smooth==2) {
      .no <- which(ds.stratum==i) 
      gam.k[c(2*.no-1,2*.no),] <- k.check.MAB(model[[i]], subsample=nrow(.ds)-lmax) # use if have >1 smooth
    } 
    # GAM SE, t-values, P-values and R^2
    .m.se <- summary(model[[i]])[["se"]]
    .m.t <- summary(model[[i]])[["p.t"]] # values not produced for gam smoothers
    .m.p <- summary(model[[i]])[["p.pv"]] # values not produced for gam smoothers, that is summary(model[[i]])$s.table[,4]
    m.r2[i,] <- summary(model[[i]])[["r.sq"]] # R2
    
    # GAM specific model parameters
    gam.concurvity[[i]] <- concurvity(model[[i]], full=T) # check concurvity
    
    print(paste0(i,'. Convergence: ',model[[i]]$converged,'. Smoothing parameter: ',model[[i]]$sp,'. Total df: ',sum(model[[i]]$edf)))
    png(file = paste0(outcome.exposure.loc.s1.smoother, .name,'- Date smoother', '.png')) # save plot for smoother
    plot(model[[i]], pages=0, main=paste0("Smoother for date - ",.name))
    dev.off()
  } else {
    if (distribution[1]=='Tweedie') { 
      var.power[[i]] <- tweedie.profile(formula=form[i,], data=.ds, p.vec=var.power.values, method='series', do.ci=F, do.smooth=tp.smooth) # estimate optimal shape parameter. Turn off CI for speed
      m.tweedie.shape[i,] <- var.power[[i]]$p.max # save shape parameter
      model[[i]] <- glm(as.formula(form[i,]), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # "na.exclude". na.omit works better for Tweedie e.g. AIC calculation
      # model.omit[[i]] <- glm(as.formula(form[i,]), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # NAs will cause tweedie AIC calculation to fail 
    } else {
      model[[i]] <- glm(as.formula(form[i,]), family=distribution, .ds, na.action="na.exclude", control=list(maxit=max.iter))
    }
    .collin <- suppressMessages(vif(model[[i]], terms='high-order')) # car::vif(model[[i]], terms='high-order') # suppress "there are higher-order terms (interactions) in this model consider setting type = 'predictor'; see ?vif"
    collin[[i]] <- data.table(Location=i, Coefficient=rownames(.collin), .collin)
    collin[[i]][!1, Location:=NA] # only keep location for 1st row
    
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
  
  # Residuals. Default residuals are deviance, but changed to quantile if Tweedie
  if (distribution[1]=='Tweedie' | isTRUE(res.quantile)) {
    res[[i]] <- qresiduals(model[[i]]) # quantile residuals for glm, recommended when deviance and Pearson residuals are grossly non-normal. Doesn't work with nb
  } else {
    res[[i]] <- residuals(model[[i]]) # default is deviance for gam and glm
  }
  
  # Variables for model checking
  # if see Warning in as.data.table.list(x, keep.rownames = keep.rownames, check.names = check.names,  : Item 1 has <4748 rows but longest item has 4748, it ooccurs in .model checks
  suppressWarnings(.cooks <- cooks.distance(model[[i]])[!is.na(res[[i]])]) # as dlnm with missing residuals, will return warning In res/(1 - hat) : longer object length is not a multiple of shorter object length
  .model.checks <- model.checks[[i]] <- data.table(res=res[[i]][!is.na(res[[i]])], fitted=fitted.values(model[[i]])[!is.na(res[[i]])], cooksd=.cooks,
                                                   Date=.d[,Date], DoW=.d[,DoW], Week=as.factor(.d[,week]), month=.d[,month], Year=.d[,Year], FYear=.d[,FYear], public.hol=.d[,public.hol],
                                                   shol=.d[,shol], PH=.d[,get('Public holiday')], day1=.d[,day1], school.hols=.d[,school.hols],
                                                   outcome=outcomes[[i]], temp=.d[,get(exposure.var)])
  
  # Aggregated model checks
  model.checks.month[[i]] <- .model.checks[, lapply(.SD, mean, na.rm=T), by=.(month,Year), .SDcols='res']
  model.checks.month[[i]][,Date:=as.Date(paste(1, month, Year, sep="."), format = "%e.%b.%Y")] # create date for 1st month
  
  model.checks.week[[i]] <- merge(model.checks[[i]][, lapply(.SD, mean, na.rm=T), by=Week, .SDcols='res'], .ds[DoW=='Mon',.(Week,Date)], all.x=T) # add Date for easier checking later
  
  # Residual plots
  # TO UPDATE: From left to right, Sun=yellow, Sat=purple, Mon=black, Tues=red, Wed=green, Thurs=blue, Fri=teal
  plot.h %<a-% {hist(.model.checks[,res], main='Residual histogram', xlab='Residuals', breaks=50); abline(v=0, lty=2, lwd=2, col='blue')} # residual distribution
  plot.qq %<a-% {qqnorm(.model.checks[,res], main="QQ-plot", las=1, col=.model.checks[,DoW]); qqline(.model.checks[,res])} # QQplot
  plot.of %<a-% {plot(x=.model.checks[,outcome], y=.model.checks[,fitted], pch=16, cex=0.5, col=.model.checks[,DoW], main='Outcome vs fitted values', ylab="Outcome", xlab="Fitted values"); abline(a=0, b=1, lty=2, lwd=2, col='blue')} # Outcome vs fitted values, with straight ilne
  plot.rf %<a-% {plot(x=.model.checks[,fitted], y=.model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs fitted values', ylab="Residuals", xlab="Fitted values"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  plot.rd %<a-% {plot(x=.model.checks[,Date], .model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs date', ylab="Residuals", xlab="Date"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  plot.re %<a-% {plot(x=.model.checks[,temp], .model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs exposure', ylab="Residuals", xlab=exposure.var); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  bplot.rd %<a-% boxplot(.model.checks[,res] ~ .model.checks[,DoW], main='Residuals vs DoW', ylab="Residuals", xlab="Day of the week")
  bplot.rm %<a-% boxplot(.model.checks[,res] ~ .model.checks[,month], main='Residuals vs month', ylab="Residuals", xlab="Day of the week")
  bplot.ry %<a-% boxplot(.model.checks[,res] ~ .model.checks[,FYear], main='Residuals vs Fyear', ylab="Residuals", xlab="Financial year")
  bplot.rp %<a-% boxplot(.model.checks[,res] ~ .model.checks[,public.hol], main='Residuals vs PH', ylab="Residuals", xlab="Public holiday?")
  bplot.ro %<a-% boxplot(.model.checks[,res] ~ .model.checks[,shol], main='Residuals vs specific PH', ylab="Residuals", xlab="Special PH")
  bplot.r1 %<a-% boxplot(.model.checks[,res] ~ .model.checks[,day1], main='Residuals vs Day 1 of month', ylab="Residuals", xlab="Day 1 of month")
  bplot.rs %<a-% boxplot(.model.checks[,res] ~ .model.checks[,school.hols], main='Residuals vs school holidays', ylab="Residuals", xlab="School holidays")
  plot.acf %<a-% {Acf(.model.checks[,res], na.action=na.omit, ylab='ACF', xlab="Lag (days)", main=''); title(main='Autocorrelation function', line=1)} # autocorrelation, without day 0. title sep to move closer to plot, as default is further away than normal
  plot.pacf %<a-% {Pacf(.model.checks[,res], na.action=na.omit, ylab='PACF', xlab="Lag (days)", main=''); title(main='Partial autocorrelation function', line=1)} # partial autocorrelation
  # plot(model[[i]], xlab='Date', ylab='Trend', se=T, main='Trend vs date') # se=T is default, plotting +/-2 SE 
  # plot(x=.model.checks[,Date], y=.model.checks[,cooksd], pch="*", cex=2, main="Cook's distance", xlab='Date', ylab="Cook's distance'"); abline(h = 4/nrow(.model.checks), col="red")  # Cook's distance
  # text(x=1:.model.checks[,.N], y=.model.checks[,cooksd], labels=.model.checks[cooksd > 4/.N, Date], col="blue")  # add labels
  plot.mrd %<a-% {plot(x=model.checks.month[[i]][,Date], model.checks.month[[i]][,res], pch=16, cex=0.5, main='Monthly residuals vs date', ylab="Residuals", xlab="Date"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  plot.wrd %<a-% {plot(x=model.checks.week[[i]][,Week], model.checks.week[[i]][,res], pch=16, cex=0.5, main='Weekly residuals vs Week', ylab="Residuals", xlab="Week"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  
  # Combine plots (large)
  png(file = paste0(outcome.exposure.loc.s1.mcheck,'zL ',.name,' mc.png')) # plot location
  par(mfrow=c(4,4), mar=rep(2.8,4), mgp=c(1.5,0.5,0)) # lower margin sizes
  plot.qq; plot.h; plot.of; plot.rf
  plot.rd; plot.re; bplot.rd; bplot.rm
  bplot.ry; bplot.rp ; bplot.ro; bplot.r1
  bplot.rs; plot.acf; plot.pacf; plot.wrd 
  dev.off()
  
  # Combine plots (small)
  png(file = paste0(outcome.exposure.loc.s1.mcheck, .name,' mc.png')) # plot location
  par(mfrow=c(3,3), mar=rep(3,4), mgp=c(1.5,0.5,0)) # lower margin sizes
  plot.qq; plot.h; plot.of
  plot.rf; plot.rd; plot.re
  bplot.rd; plot.acf; plot.pacf
  dev.off() 
  
  # Plotting residuals against covariates, to determine if modelling patterns missed 
  .check.fourier <- harmonic(.d[,Date], nfreq=6, period = nrow(.d)/.no.years)
  .check.res.coef <- coef(summary(lm(model.checks[[i]][,res] ~ log(.d[,n]) + .d[,Tue]+.d[,Wed]+.d[,Thu]+.d[,Fri]+.d[,Sat]+.d[,Sun] +
                                       .d[,public.hol] *.d[,DoW]*+ .d[,month]*.d[,DoW] + .d[,fyear]*.d[,DoW] + .d[,shol]*.d[,DoW] + .d[,day1]*.d[,DoW] + .d[,school.hols]*.d[,DoW] + 
                                       .d[,month]:.d[,fyear] + .d[,month]:.d[,day1] + .d[,month]:.d[,school.hols] +
                                       .d[,fyear]:.d[,day1] + .d[,fyear]:.d[,school.hols] + .d[,day1]:.d[,school.hols] + .check.fourier)))
  .check.res.coef.names <- rownames(.check.res.coef)
  .check.res.coef.sig <-  case_when(
    .check.res.coef[,4] < 0.001 ~ '***',
    .check.res.coef[,4] < 0.01 ~ '**',
    .check.res.coef[,4] < 0.05 ~ '*',
    .check.res.coef[,4] < 0.1 ~ '.',
    TRUE ~ as.character(NA))
  check.res.coef[[i]] <- data.table(Location=i, Coefficient=.check.res.coef.names, .check.res.coef, Significance=.check.res.coef.sig)
  check.res.coef[[i]][!1, Location:=NA] # only keep location for 1st row
  check.res.coef.sig[[i]] <- check.res.coef[[i]][,.(Coefficient, Significance)] # condensed for side-to-side comparison
  setnames(check.res.coef.sig[[i]], old=c('Coefficient','Significance'), new=c('Coefficient',paste(i,'Sig'))) # rename with i for merging
  
  if(distribution[1]=='quasipoisson') { #[1] is to avoid warnings from using multiple elements
    if(modeltype %in% c('a','gam','GAM')) {
      m.aic[i,] <- fqaic.fn(model[[i]], gam=T) # if gam, use gam option. CURRENTLY REULSTS IN INF. Not sure (I suspect not) if it applies corrected AIC, but also not sure if GAM AIC works on it
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
  # for(i in ds.stratum) { # chi-square test for over-dispersion
  #   print(c(i, round(sum(residuals(model[[i]], type="pearson")^2, na.rm=T)/df.residual(model[[i]]),3),
  #           round(pchisq(sum(residuals(model[[i]], type="pearson")^2, na.rm=T), df.residual(model[[i]]), lower.tail=F),3)))
  # }
  
  # Sum (accumulate) effects of all lags in order to eliminate one dimension of the association
  # Predicted effects: extract parameters from model corresponding to .cb variables through functions coef and vcov
  .pred1 <- crosspred(basis=.cb, model=model[[i]], cen=.cen, model.link='log') # must set log for tweedie package 
  if(str_detect(form[i,],'.hcb')) { # if above When used with .cb2, Error in crosspred(basis = .cb, model = model[[i]], cen = .cen) : coef/vcov not consistent with basis matrix. See help(crosspred). Does work with pred2 though
    .pred2 <- crosspred(basis=.hcb, model=model[[i]], cen=.cen2, model.link='log')
  }
  # .pred1 <- crosspred(basis=.cb, cen=.cen, coef=model[[i]]$coefficients, vcov=model[[i]]$coefficients, model.link='log')
  
  # Plot exposure-lag-response relationship
  png(file = paste0(outcome.exposure.loc.s1.r, .name,', e-l-r relationship.png')) # plot location
  plot(.pred1, "3d", xlab=exposure.var, ylab="Lag", zlab="Relative risk",
       main=paste0("e-l-r relationship: ",.name))
  dev.off()
  
  # Plot exposure-lag-response relationship as contour plot
  png(file = paste0(outcome.exposure.loc.s1.r, .name,', e-l-r contour plot.png'))
  plot(.pred1, "contour", xlab=exposure.var, ylab="Lag", main=paste0("e-l-r relationship: ",.name))
  dev.off()
  
  # plot(.pred1, "3d", xlab=exposure.var, ylab="Lag", zlab="Relative risk", main=paste0("e-l-r relationship: ",.name))
  
  # Plot exposure-response relationship at day 0 (no lag)
  # : allow in title names but not file names
  # add ylim for claims, but not for costs
  png(file = paste0(outcome.exposure.loc.s1.r,.name, ', e-r relationship day 0', '.png'))
  plot(.pred1, "slices", lag=0, lwd=1, col=2,  xlab=paste(exposure.var,'?(C)'), # ylim=c(0.85,1.10),
       ylab="Relative risk", main= paste0("Exposure-response relationship at day 0: ", .name)) 
  dev.off() # Save image + clear settings
  
  # Plot exposure-response relationship at all lag days
  png(file = paste0(outcome.exposure.loc.s1.r,.name, ', e-r relationships', '.png'))
  plot(.pred1, "slices", lag=0, lwd=1, col=2, xlab=paste(exposure.var,'?(C)'), # ylim=c(0.85,1.10), 
       ylab="Relative risk", main= paste0("Exposure-response relationships: ", .name))
  for(l in seq(1:lmax)) {
    lines(.pred1, "slices", lag=l, col=l+2) # for each lag day excluding lag 0 (main plot)
  }
  legend("bottomright", paste("Lag",seq(0:lmax)-1), col=(seq(0:lmax)+1), lwd=1) # ?topleft
  dev.off()
  
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
  
  png(file = paste0(outcome.exposure.loc.s1.oer, .name,', Overall e-r.png')) # plot location. File name based on heat metric
  plot(exposure.rr.s1[[i]]$temp, .exposure.rr.s1.rrfit, type="n", ylim=c(rryaxis[1], rryaxis[2]), lwd=2, col="white",
       main=i, ylab="Percent change (%)", xlab=paste(exposure.var,'(°C)'), 
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
  dev.off() # Save image + clear settings
} 
proc.time()[3]-time
# traceback()



################################################################################
# Combine and review results of Stage 1 Analysis
################################################################################

# Combine number of outcomes (OIIs and days with data), DLNM coefficients, AIC and dispersion parameter, including totals for outcomes and AIC
m1.results <- rbind(cbind(no.oii, no.outcome, coef, m.aic, m.r2, m.devexp, m.dispersion), c(sum(no.oii),sum(no.outcome),rep(NA, ncol(coef)),sum(m.aic), rep(NA,3)))
colnames(m1.results) <- c('n','Number of days with outcome',colnames(coef),'AIC','R^2','Deviance explained','Dispersion parameter')
write.csv(m1.results, file=paste0(outcome.exposure.loc.s1,'Outcome, coefficients and model fit.csv'), na='', row.names=T) # create csv file
save(m1.results, file=paste0(outcome.exposure.loc.s1,'Outcome, coefficients and model fit.rda')) # save dataset for easier access
write.csv(sum(m.aic), file=paste0(outcome.exposure.loc.s1,'Total AIC.csv'), na='', row.names=T)
write.csv(sum(no.oii), file=paste0(outcome.exposure.loc.s1,'Total OII.csv'), na='', row.names=T)

# Output model coefficients and RRs
write.csv(do.call(rbind, m.coef), file=paste0(outcome.exposure.loc.s1,'Coefficients.csv'), na='') # Model coefficients, SE, 95% CI and t-tests

# Modelling against residuals. Not output varies
check.res.coef.s <- check.res.coef.sig %>% purrr::reduce(full_join, by = "Coefficient") # p-values only, next to each other
check.res.coef.f <- do.call(rbind, check.res.coef) # all measurements, arranged vertically
write.csv(check.res.coef.s, file=paste0(outcome.exposure.loc.s1.mcheck,'zLM residuals with covariates.csv'), na='', row.names=F)
write.csv(check.res.coef.f, file=paste0(outcome.exposure.loc.s1.mcheck,'zLM residuals with covariatesf.csv'), na='', row.names=F) # f short for full

# Collinearity or GAM specific checks
if(modeltype %in% c('a','gam','GAM')) {
  write.csv(gam.k, file=paste0(outcome.exposure.loc.s1,'k-index.csv'), na='', row.names=T) # create csv file
  write.csv(do.call(rbind, gam.concurvity), file=paste0(outcome.exposure.loc.s1,'Concurvity.csv'), na='', row.names=F)
} else {
  write.csv(do.call(rbind, collin), file=paste0(outcome.exposure.loc.s1,'Collinearity.csv'), na='', row.names=F)
}

# Tweedie shape parameters
if (distribution[1]=='Tweedie') { 
  write.csv(m.tweedie.shape, file=paste0(outcome.exposure.loc.s1,'Tweedie shape parameters.csv'), na='', row.names=T)
} 


######################### END ############################
######################### END ############################
