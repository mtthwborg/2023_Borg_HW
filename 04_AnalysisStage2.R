
##################################################
### Meta-analysis: results of reduced coef and compute BLUPS (best linear unbiased prediction)
##################################################

exposure.country.mean <- colMeans(exposure.by) # could arguably weight by area number of years (study frame), but Gasaparrini didn't do this either and impact likely negligible

## Meta-analysis. Default methods is REML
exposure.mean <- exposure[,1] # mean exposure
mv <- mixmeta(coef~1, vcov)
mv.results <- mv.results.fn(mv) 
write.csv(mv.results, file=paste0(s2results, 'MM tests.csv'), na='', row.names=T) # create csv file


## Pool estimated overall cumulative exposure-response associations to obtain overall estimation at each location level for main model
blups <- blup(mv, vcov=T) # BLUPs. Element for each City



##################################################
### Define minimum occupational OII values (not used in this analysis, but still calculated should you prefer this approach)
##################################################

# Generate matrix for storing re-centred results
minpercity <- mintempcity <- rep(NA,length.ds.city)
names(mintempcity) <- names(minpercity) <- ds.city

# Define minimum injuries values, excluding very low and hot Ts (<1st and >99th percentiles)
for(i in ds.city) {
  .ds <- daily.ds[City==i] # dataset for each unique value # daily.ds[City=='Hobart: Indoor']
  
  # Use same exposure-response relationship, but for exposure percentiles 1-99
  .predvar <- quantile(temps[[i]], 1:99/100,na.rm=T)
  .argvar.bound <- list(x=.predvar, fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T), Bound=range(temps[[i]], na.rm=T)) # list of arguments
  
  # Redefine function using all arguments, boundary knots included
  .blups <- blups[[which(ds.city==i)]] 
  .bvar.bound <- do.call(onebasis, .argvar.bound) # basis for exposure-response
  .bvar1 <- .bvar.bound%*%.blups$blup # 
  minpercity[i] <- (1:99)[which.min(.bvar1)] # if stop using 1:99, may need to change code so that position links to actual percentile
  mintempcity[i] <- quantile(temps[[i]], minpercity[i]/100, na.rm=T)
}



################################################################################
# Predict pooled overall cumulative associations
################################################################################

## Mean values of mixmeta predictors as a data frame to combine. In this case no meta predictors
datanew <- as.data.frame(t(rep(1,length.ds.city)))
mvpred <- predict.mixmeta(mv, datanew, vcov=T, format="list")

## Define exposure percentile with lowest incidence of injuries
bvar <- crossbasis(exposure.country.mean, argvar=list(fun=espline, knots=exposure.country.mean[paste0(eknots*100)]), arglag=li_arglag) # crossbasis
bvar1 <- bvar%*%mvpred$fit # multiply coefficients for combined effect
bvar1 <- bvar1[order(bvar1[,1]),, drop=F] # sort in ascending order (not required, but easier for visualisation
bvar1 <- bvar1[between(as.numeric(rownames(bvar1)), cenpen.minmax[1], cenpen.minmax[2]),] # exclude percentiles outside range for selecting optimal percentile, converts to vector
cenindcountry <- as.numeric(names(bvar1[which.min(bvar1)])) # exposure percentile with lowest incidence of outcome (optimal percentile)

# Define centering percentile for country, using either median or above exposure percentile 
if (is.null(cenpen)) {
  cenpercountry <- pmin(pmax(predper[cenindcountry],cenpen.minmax[1]),cenpen.minmax[2]) # min/max is 10th/90th. chooses 10, as incidence lowest at 10
} else if(is.numeric(cenpen)) {
  cenpercountry <- cenpen # use defined centering value
} else if(cenpen %in% c('mean','Mean','average','Average')) {
  cenpercountry <- 'mean' # to be used as code for mean
} else if(cenpen=='ehf') {
  cenpercountry <- 'ehf' # to be used as code for mean
}

if(is.numeric(cenpercountry)) {
  centre <- exposure.country.mean[paste0(cenpercountry)] # centered prediction (from corresponding percentile)
} else if(cenpercountry=='ehf') {
  centre <- 0 # centred on 0
} else {
  centre <- sum(sapply(temps, sum)) / sum(sapply(temps, length)) # mean exposure across all strata
}

# Predict pooled overall cumulative associations (RR) for main model
exposure_ehf_overall <- colMeans(exposure_ehf)
ehfmean <- sum(sapply(temps, sum)) / sum(sapply(temps, length))
cp <- crosspred(bvar, coef=mvpred$fit, vcov=mvpred$vcov, model.link="log", at=c(exposure.country.mean,exposure_ehf_overall,ehfmean), cen=centre)
save(cp, file=paste0(s2results, 'cp.rda')) # save for fwald1



################################################################################
# LAG-RESPONSE ASSOCIATIONS, USING SELECTED CENTERING VALUE
################################################################################

# Objects to store overall cumulative exposure results (lag-response association at set percentiles vs centering value)
cvlag.length <- length(crossreduce(.cb, model=model[[i]], "var", value=exposure.by[1,'1'], cen=.cen, model.link='log')$coefficients) # only interested in length, uses last cb from S1 but that is irrelevant

coeflag_ehf1 <- coeflag_ehf2 <- matrix(NA, length.ds.city, cvlag.length, dimnames=list(ds.city)) # not really sure why columns is 3. I though 8 was for columns for lag 0 + each lag day. CHANGED WHEN I INCREASED PARAMETES FROM edf TO edf + 1
cpmodel <- vcovlag_ehf1 <- vcovlag_ehf2 <- vector("list", length.ds.city)
names(cpmodel) <- names(vcovlag_ehf1) <- names(vcovlag_ehf2) <- ds.city

# Run model for each by.vars with re-centered values to obtain lag-reponse relationship at the 1st and 99th percentile

for(i in ds.city) {
  # print(paste('Stage 2 individual models:',i))
  .ds <- daily.ds[City==i] # dataset for each unique value
  .name <- unique(.ds$City) # names with all of a, b and c together
  
  .outcome <- .ds[,get(outcome.var)] # outcome
  
  # Calculate crossbasis (comments in Gasparrini and Martinez-Solanas code state it is centred, although this isn't the case). Model appears same, with changes made for crosspred and crossreduce with cen 
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag) # lag=lag2
  
  # Model using crossbasis
  .no.years <- no.years[i,]
  
  # Predictions and reduction to lag-response at 1st (extreme cold) and 99th (extreme hot) percentiles with new centering (changes coef and vcov)
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T) # quantile
  }  else {
    .cen <- exposure.mean[i] # mean
  }
  
  # Predictions and reduction to lag-response at set percentiles. Centering required as it changes coef-vcov
  .perc <-  exposure.by[i,] 
  
  # Severe heatwave threshold
  redlag_ehf1 <- crossreduce(.cb, model=model[[i]], "var", value=exposure_ehf[i,2], cen=.cen, model.link='log')
  coeflag_ehf1[i,] <- coef(redlag_ehf1)
  vcovlag_ehf1[[i]] <- vcov(redlag_ehf1)
  # Extreme (>=2*EHF85) heatwave threshold
  redlag_ehf2 <- crossreduce(.cb, model=model[[i]], "var", value=exposure_ehf[i,3], cen=.cen, model.link='log')
  coeflag_ehf2[i,] <- coef(redlag_ehf2)
  vcovlag_ehf2[[i]] <- vcov(redlag_ehf2)
}

# Run meta-analysis with lag models 
mvlag_ehf1 <- mixmeta(coeflag_ehf1~1, vcovlag_ehf1, data=list(ds.city))
mvlag_ehf2 <- mixmeta(coeflag_ehf2~1, vcovlag_ehf2, data=list(ds.city)) 

# Predict pooled coefficients
mvpredlag_ehf1 <- predict(mvlag_ehf1,datanew,vcov=T,format="list")
mvpredlag_ehf2 <- predict(mvlag_ehf2,datanew,vcov=T,format="list")

# Obtain predictions for lag 0 to "lmax"
blag <- do.call(onebasis,c(list(x=seq(0,lmax)),attr(.cb,"arglag"))) # uses CB attributes (same for all by.vars)

# Predict pooled lag-response associations for heat and cold
cplag_ehf1 <- crosspred(blag,coef=mvpredlag_ehf1$fit,vcov=mvpredlag_ehf1$vcov, model.link="log", at=0:lmax)
cplag_ehf2 <- crosspred(blag,coef=mvpredlag_ehf2$fit,vcov=mvpredlag_ehf2$vcov, model.link="log", at=0:lmax)



################################################################################
# Main results: RRs (95% CI) of set exposure percentiles compared to minimum occupational-injuries percentile (MOIP)
################################################################################

# Short percentile
predprer.short.rep <- 1:length(predper.short) # number of percentiles
results.short <- matrix(NA, nrow=length(predper.short), ncol=4) # create empty matrix of results
results.short[,1] <- predper.short # 1st column is percentile values
results.short[predprer.short.rep,2:4] <- c(cp$allRRfit[predper %in% predper.short], # column 2 is RR corresponding to each percentile (row), 
                                           cp$allRRlow[predper %in% predper.short], # column 3 is RRlow
                                           cp$allRRhigh[predper %in% predper.short]) # column 4 is RRhigh
results.short <- round(results.short, digits=3)



###############################################################################
# COMPUTE ATTRIBUTABLE INJURIES FOR EACH CITY, WITH EMPIRICAL CI ESTIMATED USING RE-CENTERED BASES
################################################################################

# Create vectors to store total injuries, accounting for missing
totclaims <- rep(NA,length.ds.city)
names(totclaims) <- ds.city

# Objects to store attributable injuries (simulations)
sim.names <- c("Non-heatwave days","Heatwaves","Low-intensity heatwaves","Severe heatwaves","Extreme heatwaves")
length.sim <- length(sim.names)
matsim <- matrix(NA, length.ds.city, length.sim, dimnames=list(ds.city, sim.names)) # matrix: attributable injuries
arraysim <- array(NA, dim=c(length.ds.city, length.sim, nsim), dimnames=list(ds.city, sim.names)) # array: attributable injuries CI

# Run loop
gc()
for(i in ds.city){
  print(paste('Attributable risk:',i))
  
  .ds <- daily.ds[City==i] # dataset for each unique value
  .outcome <- .ds[,get(outcome.var)] # .define outcome
  
  # Predictions and reduction to lag-response at 1st (extreme cold) and 99th (extreme hot) percentiles with new centering (changes coef and vcov)
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else if(cenpercountry=='ehf') {
    .cen <- 0 # centred on 0
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  # Derive cross-basis. NB: centering point different than original choice of 75th (i see no difference)
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag)
  
  # Compute attributable numbers with reduced coefficients, based on centering value
  .perc <-  exposure.by[i,] # all percentiles and related average exposure
  .blups <- blups[[which(ds.city==i)]]
  
  matsim[i,sim.names[1]] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                   range=c(-100,.cen))
  matsim[i,sim.names[2]] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                   range=c(.cen,100))
  matsim[i,sim.names[3]] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                   range=c(.cen,exposure_ehf[i,2]))
  matsim[i,sim.names[4]] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                   range=c(exposure_ehf[i,2],exposure_ehf[i,3]))
  matsim[i,sim.names[5]] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                   range=c(exposure_ehf[i,3],100))
  
  arraysim[i,sim.names[1],] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                      range=c(-100, .cen), sim=T, nsim=nsim)
  arraysim[i,sim.names[2],] <-  attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                       range=c(.cen,100), sim=T, nsim=nsim)
  arraysim[i,sim.names[3],] <-  attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                       range=c(.cen,exposure_ehf[i,2]), sim=T, nsim=nsim)
  arraysim[i,sim.names[4],] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                      range=c(exposure_ehf[i,2],exposure_ehf[i,3]), sim=T, nsim=nsim)
  arraysim[i,sim.names[5],] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir=attrdl.dir, cen=.cen,
                                      range=c(exposure_ehf[i,3],100), sim=T, nsim=nsim)

  totclaims[i] <- sum(.outcome, na.rm=T) # store total injuries (account for missing)
}


### Attributable numbers ###

if (type.outcome=='oi') { # if number of claims/injuries/diseases
  afyr.r <- 2 # round yearly results to 2 df to show Darwin stats clearly
} else {
  afyr.r <- 1 # round yearly results to 1 df to show Darwin stats clearly
}

## City-specific attributable numbers. Already have row names
ancitylow <- apply(arraysim,c(1,2),quantile,0.025) # quantile across dimensions 1 and 2
ancityhigh <- apply(arraysim,c(1,2),quantile,0.975)
matsim.yr <- sweep(matsim,MARGIN=1,no.years,'/') # AN per year, factoring in different years per city (Hobart)
ancitylow.yr <- sweep(ancitylow,MARGIN=1,no.years,'/')
ancityhigh.yr <- sweep(ancityhigh,MARGIN=1,no.years,'/')

ancity <- matrix(paste0(round.fn(matsim),' (',round.fn(ancitylow),' to ',round.fn(ancityhigh),')'), ncol=length.sim) # insert comma every 3 digits, no rounding
ancity.yr <- matrix(paste0(round.fn(matsim.yr,afyr.r),' (',round.fn(ancitylow.yr,afyr.r),' to ',round.fn(ancityhigh.yr,afyr.r),')'), ncol=length.sim) # insert comma every 3 digits, no rounding
rownames(ancity) <- rownames(ancity.yr) <- ds.city

## Total attributable numbers
antot <- colSums(matsim) # sum through strata
antotlow <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.025) # sum by dimensions 2 and 3 (across all other dimension i.e. 1st), then quantile by new dimension 1 (old 2nd dimension) (quantile across new dimension 2 which is old 3rd dimension)
antothigh <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.975)
antot.yr <- colSums(matsim.yr) # per year results
antotlow.yr <- colSums(ancitylow.yr) # per year results. quantile calculation slightly differs from above, but necessary as above approach doesn't consider different years per city (Hobart)
antothigh.yr <- colSums(ancityhigh.yr) # per year results

antotal <- matrix(paste0(round.fn(antot),' (',round.fn(antotlow),' to ',round.fn(antothigh),')'), ncol=length.sim)
antotal.yr <- matrix(paste0(round.fn(antot.yr,afyr.r),' (',round.fn(antotlow.yr,afyr.r),' to ',round.fn(antothigh.yr,afyr.r),')'), ncol=length.sim)
rownames(antotal) <- rownames(antotal.yr) <- 'Total'
colnames(ancity) <- colnames(ancity.yr) <- sim.names # colnames

## Export attributable numbers
ancountry <- rbind(antotal, ancity) # ANs without extra stratification
ancountry.yr <- rbind(antotal.yr, ancity.yr) # ANs without extra stratification
write.csv(ancountry, file=paste0(s2results, 'Numbers.csv'), na='', row.names=T) # create csv file
write.csv(ancountry.yr, file=paste0(s2results, 'Numbers per year.csv'), na='', row.names=T) # create csv file


### Attributable fractions ###
totclaimtot <- sum(totclaims) # total injuries

# City-specific
afcit <- matsim/totclaims*100
afcitylow <- ancitylow/totclaims*100
afcityhigh <- ancityhigh/totclaims*100
afcity <- matrix(paste0(round.fn(afcit, af.round),' (',round.fn(afcitylow, af.round),' to ',round.fn(afcityhigh, af.round),')'), ncol=length.sim)
rownames(afcity) <- ds.city

# Total
aftot <- antot/totclaimtot*100
aftotlow <- antotlow/totclaimtot*100
aftothigh <- antothigh/totclaimtot*100
aftotal <- matrix(paste0(round.fn(aftot, af.round),' (',round.fn(aftotlow, af.round),' to ',round.fn(aftothigh, af.round),')'), ncol=length.sim)
rownames(aftotal) <- 'Total'
colnames(afcity) <- colnames(aftotal) <- sim.names # colnames

# Export attributable fractions
afcountry <- rbind(aftotal, afcity) # AFs without extra stratification
write.csv(afcountry, file=paste0(s2results, 'Fractions.csv'), na='', row.names=T) # create csv file

# Export AF and AN together
save(ancountry, ancountry.yr, afcountry, file=paste0(s2results, 'acountry.rda')) # save dataset after all changes for easier access



################################################################################
# Plots: overall relationship
################################################################################

indlab <- predper %in% c(0,2.5,10,25,50,75,90,75,97.5,100) # requires predper, so must be in S2

# Plot: Overall cumulative exposure-response association (all exposure values, lag reduced)
oer.yaxis <- seq(0.9,1.5,by=0.1)

png(file = paste0(s2results, 'Overall e-r.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
par(mar=c(4.1,3,1.6,0)) # inner graph margins, as much whitespace removed as possible
plot(cp,lwd=2,col="white",yaxt='n', ylim=c(min(oer.yaxis), max(oer.yaxis)), xlab='', ylab='') # str_remove_all(paste0('Percent change in ',tolower(outcome.var)," (%)"),'\\(000s\\)') # ylim=c(floor(min(cp[["allRRfit"]])*10)/10-0.05, ceiling(max(cp[["allRRfit"]])*10)/10+0.05)
ind1 <- cp$predvar<=exposure_ehf_overall[1]
ind2 <- cp$predvar>=exposure_ehf_overall[1]
ind85 <- cp$predvar>=exposure_ehf_overall[2]
ind852 <- cp$predvar>=exposure_ehf_overall[3]
lines(cp$predvar[ind1],cp$allRRfit[ind1],col=hwcolours[1],lwd=2)
lines(cp$predvar[ind2],cp$allRRfit[ind2],col=hwcolours[3],lwd=2) # low-intensity heatwave
lines(cp$predvar[ind85],cp$allRRfit[ind85],col=hwcolours[4],lwd=2) # severe heatwave
lines(cp$predvar[ind852],cp$allRRfit[ind852],col=hwcolours[5],lwd=2) # extreme heatwave
title(xlab=paste(exposure.var, "(°K)"), line=2) 
abline(v=c(exposure_ehf_overall[2:3]), lty=c(3,3)) # severe and extreme heatwaves
axis(2, at=oer.yaxis, labels=(oer.yaxis-1)*100, las=1, mgp=c(2.5,0.8,0))
title(main=str_replace(str_remove_all(outcome.var,'\\(000s\\)'),'Number of illnesses','Number of injuries and illnesses'), line=nline) # title, move closer to graph
title(ylab='Percent change (%)', line=2.1) # title, move closer to graph
abline(v=cp$cen,lty=2) # centre line
dev.off() # Save image + clear settings


# Plot: overall cumulative lag-response associations
.cexaxis <- 1
lag.yaxis <- seq(0.97,1.05,by=0.01)
.axis %<a-% {axis(2,at=lag.yaxis, labels=(lag.yaxis-1)*100, cex.axis=.cexaxis, las=1)} # shared axis as percentage change, tick labels are horizontal
lag.ylim <- c(min(lag.yaxis),max(lag.yaxis))

png(file = paste0(s2results, 'Overall lag.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
par(mfrow=c(1,2), mar=c(3,2.5,1,1), oma=c(0,0,0,0), mgp=c(1.5,0.5,0)) # 3*3, space between borders and text

plot(cplag_ehf1,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col=hwcolours[4], ci.arg=list(density=20,col=hwcolours[4])) # practically identical to mean, thus no relationship apparent, but included for comparison and c(3,3) instead of c(2,4)
title(main='Severe', line=tline)
.axis
plot(cplag_ehf2,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col=hwcolours[5], ci.arg=list(density=20,col=hwcolours[5])) # practically identical to mean, thus no relationship apparent, but included for comparison and c(3,3) instead of c(2,4)
title(main='Extreme', line=tline)
.axis
dev.off() # Save image + clear settings




################################################################################
#### Plots: overall cumulative exposure-response associations by by-variables
################################################################################

.exposure.rr <- list()
s2city <- paste0(s2results,'City level/')
dir.create(paste0(s2results,'City level')) # Warning if exists (doesn't replace)

# Loop over each strata. Plots are individualised
for(i in ds.city) {
  .ds <- daily.ds[City==i] # dataset for each unique value
  .name <- unique(.ds$City) # names with all of a, b and c together
  
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else if(cenpercountry=='ehf') {
    .cen <- 0 # centred on 0
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  .argvar <- list(x=temps[[i]], fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T))
  
  .bvar <- do.call(onebasis, .argvar)
  .blups <- blups[[which(ds.city==i)]] # get correct blups
  .pred2 <- crosspred(.bvar, coef=.blups$blup, vcov=.blups$vcov, model.link="log", by=0.1, cen=.cen) #  overall cumulative exposure-response relationship, for each percentile. Specific humidity fails at this step because of Error in seq.default(from = min(pretty), to = to, by = by) : 'from' must be a finite number
  
  .exposure.rr[[i]] <- data.frame(matrix(nrow=length(.pred2$predvar), ncol=4))
  colnames(.exposure.rr[[i]]) <- c("temp", "RRfit", "RRlow", "RRhigh")
  .exposure.rr[[i]]$temp <- .pred2$predvar
  .exposure.rr[[i]]$RRfit <- (.pred2$allRRfit-1)*100
  .exposure.rr[[i]]$RRlow <- (.pred2$allRRlow-1)*100
  .exposure.rr[[i]]$RRhigh <- (.pred2$allRRhigh-1)*100
  
  # Values above or below centering value (for red/blue colour)
  .ind1 <- .pred2$predvar<=.cen
  .ind2 <- .pred2$predvar>=.cen
  .ind85 <- .pred2$predvar>=exposure_ehf[i,2] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  .ind852 <- .pred2$predvar>=exposure_ehf[i,3] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  
  # Overall cumulative exposure-response relationship plot
  # Could technically extract specific RRs values from RRfit/RRlow/RRhigh, but AF more useful as it covers a range of RRs
  png(file = paste0(s2city, .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=rryaxis, lwd=2, col="white",
       main=i, ylab="Percent change (%)", xlab=paste(exposure.var,'(°C)'),
       cex.main=cexmain, cex.lab=cexlab, cex.axis=cexaxis, lab=c(6,5,7)) # main plot
  .erplot %<a-% { # save main plot additions
    polygon(c(.exposure.rr[[i]]$temp,rev(.exposure.rr[[i]]$temp)),c(.exposure.rr[[i]]$RRlow,rev(.exposure.rr[[i]]$RRhigh)), col="grey89", border=F) # 95% CI envelope
    lines(.exposure.rr[[i]]$temp[.ind1],.exposure.rr[[i]]$RRfit[.ind1],col=hwcolours[1],lwd=2); # cold, left of cen
    lines(.exposure.rr[[i]]$temp[.ind2],.exposure.rr[[i]]$RRfit[.ind2],col=hwcolours[3],lwd=2); # heatwave
    lines(.exposure.rr[[i]]$temp[.ind85],.exposure.rr[[i]]$RRfit[.ind85],col=hwcolours[4],lwd=2); # severe heatwave
    lines(.exposure.rr[[i]]$temp[.ind852],.exposure.rr[[i]]$RRfit[.ind852],col=hwcolours[5],lwd=2); # extreme heatwave
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure_ehf[i,c(2,3)]),lty=3)
  }
  .erplot # insert main plot additions
  .lines # insert lines
  dev.off() # Save image + clear settings
  
  # Plot + histogram (slice plot)
  png(file = paste0(s2city, 'z', .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=c(rryaxis[1]-50,rryaxis[2]), lwd=2, col="white", ylab="Percent change (%)", xlab='',
       mgp=c(1.3,0.4,0), cex.lab=cexlab*0.7, cex.axis=cexaxis*0.7, yaxt='n', axes=F) # main plot, but no axes # cexlab*.75, cexaxis*.75
  .erplot
  title(main=i, line=tline*2, cex.main=cexmain*0.7) # closer to plot
  title(xlab=paste(exposure.var,'(°C)'), line=1.1, cex.lab=cexlab*0.7) # xtitle, moved away from graph
  axis(1, col.axis="black", cex.axis=cexaxis*0.7, tck=tck.length, mgp=c(2.5,0.3,0)) # x-axis
  axis(2, at=seq(rryaxis[1],rryaxis[2],by=20), col.axis="black", cex.axis=cexaxis*0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # y-axis on L
  
  .breaks <- c(min(temps[[i]],na.rm=T)-1, seq(.pred2$predvar[1], .pred2$predvar[length(.pred2$predvar)],length=30), max(temps[[i]],na.rm=T)+1)
  hist <- hist(temps[[i]],breaks=.breaks,plot=F) # histogram
  hist$density <- hist$density/max(hist$density)*0.7 # density value
  prop <- max(hist$density)/max(hist$counts) # convert proportion so that it fits on same scale as other plot
  counts <- pretty(hist$count,3) # counts per bin
  par(new=TRUE) # add to graph
  plot(hist,ylim=c(0,max(hist$density)*2.0),axes=F,ann=F,col=grey(0.95),freq=F) # plot histogram
  axis(4, at=counts*prop, labels=counts, cex.axis=0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # axis on R indicating bin counts
  .lines
  dev.off() # Save image + clear settings
}

# Loop over each strata. Plots are combined
png(file = paste0(s2city, 'zOverall e-rs.png'), res=350, width=1700, height=2400) # plot location # paste0('/Users/MatthewBorg/', 'zOverall e-rs.png')
layout(matrix(c(1:7,0),ncol=2,byrow=T))
par(mfrow=c(4,2), mar=c(2.1,2.2,1,1.5), oma=c(0,0,0,0), mgp=c(1.5,0.4,0), las=1) # 4*4, space between borders and text

for(i in ds.city) {
  .ds <- daily.ds[City==i] # dataset for each unique value
  .name <- unique(.ds$City) # names with all of a, b and c together
  # .name <- str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(.name, ',',''), 'Adelaide','Ade'), 'Brisbane','Bri'), 'Canberra','Can'), 'Darwin','Dar'), 'Hobart','Hob'), 'Melbourne','Mel'), 'Perth','Per'), 'Sydney','Syd')
  
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else if(cenpercountry=='ehf') {
    .cen <- 0 # centred on 0
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  .argvar <- list(x=temps[[i]], fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T))
  
  .bvar <- do.call(onebasis, .argvar)
  .blups <- blups[[which(ds.city==i)]] # get correct blups
  .pred2 <- crosspred(.bvar, coef=.blups$blup, vcov=.blups$vcov, model.link="log", by=0.1, cen=.cen) #  overall cumulative exposure-response relationship, for each percentile. Specific humidity fails at this step because of Error in seq.default(from = min(pretty), to = to, by = by) : 'from' must be a finite number
  
  # Values above or below centering value (for red/blue colour)
  .ind1 <- .pred2$predvar<=.cen
  .ind2 <- .pred2$predvar>=.cen
  .ind85 <- .pred2$predvar>=exposure_ehf[i,2] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  .ind852 <- .pred2$predvar>=exposure_ehf[i,3] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  
  
  # Overall cumulative exposure-response relationship plot
  .erplot %<a-% { # save main plot additions
    polygon(c(.exposure.rr[[i]]$temp,rev(.exposure.rr[[i]]$temp)),c(.exposure.rr[[i]]$RRlow,rev(.exposure.rr[[i]]$RRhigh)), col="grey89", border=F) # 95% CI envelope
    lines(.exposure.rr[[i]]$temp[.ind1],.exposure.rr[[i]]$RRfit[.ind1],col=hwcolours[1],lwd=2); # cold, left of cen
    lines(.exposure.rr[[i]]$temp[.ind2],.exposure.rr[[i]]$RRfit[.ind2],col=hwcolours[3],lwd=2); # heatwave
    lines(.exposure.rr[[i]]$temp[.ind85],.exposure.rr[[i]]$RRfit[.ind85],col=hwcolours[4],lwd=2); # severe heatwave
    lines(.exposure.rr[[i]]$temp[.ind852],.exposure.rr[[i]]$RRfit[.ind852],col=hwcolours[5],lwd=2); # extreme heatwave
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure_ehf[i,c(2,3)]),lty=3)
  }
  
  # Plot + histogram (slice plot)
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=c(rryaxis[1]-50,rryaxis[2]), lwd=2, col="white", ylab="Percent change (%)", xlab='',
       mgp=c(1.3,0.4,0), cex.lab=cexlab*0.7, cex.axis=cexaxis*0.7, yaxt='n', axes=F) # main plot, but no axes # cexlab*.75, cexaxis*.75
  .erplot
  title(main=i, line=tline*2, cex.main=cexmain*0.7) # closer to plot
  title(xlab=paste(exposure.var,'(°C)'), line=1.1, cex.lab=cexlab*0.7) # xtitle, moved away from graph
  axis(1, col.axis="black", cex.axis=cexaxis*0.7, tck=tck.length, mgp=c(2.5,0.3,0)) # x-axis
  axis(2, at=seq(rryaxis[1],rryaxis[2],by=20), col.axis="black", cex.axis=cexaxis*0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # y-axis on L
  
  .breaks <- c(min(temps[[i]],na.rm=T)-1, seq(.pred2$predvar[1], .pred2$predvar[length(.pred2$predvar)],length=30), max(temps[[i]],na.rm=T)+1)
  hist <- hist(temps[[i]],breaks=.breaks,plot=F) # histogram
  hist$density <- hist$density/max(hist$density)*0.7 # density value
  prop <- max(hist$density)/max(hist$counts) # convert proportion so that it fits on same scale as other plot
  counts <- pretty(hist$count,3) # counts per bin
  par(new=TRUE) # add to graph
  plot(hist,ylim=c(0,max(hist$density)*2.0),axes=F,ann=F,col=grey(0.95),freq=F) # plot histogram
  axis(4, at=counts*prop, labels=counts, cex.axis=0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # axis on R indicating bin counts
  .lines
  # dev.off() # Save image + clear settings
}
dev.off()

# Export exposure RRs
exposure.rr <- do.call(rbind,.exposure.rr) # convert to df
exposure.rr <- cbind(gsub("\\..*","",rownames(exposure.rr)), exposure.rr) # combine with rownames formatted to remove .
rownames(exposure.rr) <- NULL
colnames(exposure.rr) <- c('City',exposure.var,'RRfit','RRlow','RRhigh')
write.csv(exposure.rr, file=paste0(s2results, 'Exposure values RR.csv'), na='', row.names=F) # create csv file



################################################################################
# CITY OVERALL E-R CURVES
################################################################################

# Plot: Overall cumulative exposure-response association (all exposure values, lag reduced)

png(file = paste0(s2city, 'City vs BLUP.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
par(mar=c(4.1,3,1.6,0)) # inner graph margins, as much whitespace removed as possible
layout(matrix(1:2,1,2))

plot(cp,lwd=2,yaxt='n', ylim=c(min(oer.yaxis.blup),max(oer.yaxis.blup)), xlab='', ylab='', ci="n") # str_remove_all(paste0('Percent change in ',tolower(outcome.var)," (%)"),'\\(000s\\)') # ylim=c(floor(min(cp[["allRRfit"]])*10)/10-0.05, ceiling(max(cp[["allRRfit"]])*10)/10+0.05)
title(xlab=paste(exposure.var,'(°K)'), line=2) 
abline(v=c(cp$cen), lty=2, lwd=0.8) # severe and extreme heatwaves
# abline(v=c(cp$cen, exposure_ehf_overall[2:3]), lty=c(2,3,3), lwd=c(0.8,0.8,0.8)) # severe and extreme heatwaves
axis(2, at=oer.yaxis.blup, labels=(oer.yaxis.blup-1)*100, las=1, mgp=c(2.5,0.8,0))
for(i in ds.city) {
  lines(crosspred(bvar,coef=coef[i,],vcov=vcov[[i]], model.link="log",cen=centre),lty=6,lwd=1.5,col=oer.colours[which(ds.city==i)])
} 
title(main='Study-specific', line=nline) # title, move closer to graph
title(ylab='Percent change (%)', line=2.1) # title, move closer to graph
legend(legend=ds.city,'bottomright',col=oer.colours, lty=1, cex=0.45)

plot(cp,lwd=2,yaxt='n', ylim=c(min(oer.yaxis.blup),max(oer.yaxis.blup)), xlab='', ylab='', ci="n") # str_remove_all(paste0('Percent change in ',tolower(outcome.var)," (%)"),'\\(000s\\)') # ylim=c(floor(min(cp[["allRRfit"]])*10)/10-0.05, ceiling(max(cp[["allRRfit"]])*10)/10+0.05)
title(xlab=paste(exposure.var,'(°K)'), line=2) 
abline(v=c(cp$cen), lty=2, lwd=0.8) # severe and extreme heatwaves
# abline(v=c(cp$cen, exposure_ehf_overall[2:3]), lty=c(2,3,3), lwd=c(0.8,0.8,0.8)) # severe and extreme heatwaves
axis(2, at=oer.yaxis.blup, labels=(oer.yaxis.blup-1)*100, las=1, mgp=c(2.5,0.8,0))
for(i in ds.city) {
  lines(crosspred(bvar,coef=blups[[which(ds.city==i)]]$blup, vcov=blups[[which(ds.city==i)]]$vcov, model.link="log",cen=centre),lty=6,lwd=1.5,col=oer.colours[which(ds.city==i)])
} 
title(main='BLUPs', line=nline) # title, move closer to graph
title(ylab='Percent change (%)', line=2.1) # title, move closer to graph
legend(legend=ds.city,'bottomright',col=oer.colours, lty=1, cex=0.45)

dev.off() # Save image + clear settings



######################### END ############################
######################### END ############################

  