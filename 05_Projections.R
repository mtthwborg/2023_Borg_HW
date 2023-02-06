
########################################################################################################
# Prepare projected data. No leap days, like in Vicedo-Cabrera example
########################################################################################################

load(paste0(cciaf.loc,'ccia_future2.rda')) # load data

## Define projected outcome series as historical average per day of year, then repeated along same projection period of modelled temperatures series
ccia_future2[,Month:=month(Date)] # define Month to limit to warm season. much faster to do outside of loop
ccia_future2[RCP %in% c('rcp26','RCP26'),RCP:='RCP2.6'] # format RCP
ccia_future2[RCP %in% c('rcp45','RCP45'),RCP:='RCP4.5']
ccia_future2[RCP %in% c('rcp60','RCP60'),RCP:='RCP6.0']
ccia_future2[RCP %in% c('rcp85','RCP85'),RCP:='RCP8.5']

## EHF projected data, assuming no adaptation. Set to use EHF automatically
fdmt <- ccia_future2[Month %in% c(10,11,12,1,2,3) & City!='Canberra' & RCP %in% c('RCP4.5','RCP8.5'), # no restrictions on ehf
                     .(`Average DMT (SD)`=paste0(round.fn(mean(tas_ehf, na.rm=T),2),' (',round.fn(sd(tas_ehf, na.rm=T),2),')')),
                     by=c('City','Model','RCP','Period')]
fehf <- unique(ccia_future2[Month %in% c(10,11,12,1,2,3) & City!='Canberra' & RCP %in% c('RCP4.5','RCP8.5') & ehf>0, # automatically only include EHF>0
                            .(`DMT~95~`=round.fn(`dmt95`,2),
                              `HW days per year`=round.fn(sum(ehf>0)/30,1), #`Mean EHF~p~`=round.fn(mean(ehf),2),
                              `EHF~50p~`=round.fn(median(ehf),2),
                              `EHF~85p~`=round.fn(quantile(ehf, 0.85),2)),
                            by=c('City','Model','RCP','Period')])
fehfa <- unique(ccia_future2[Month %in% c(10,11,12,1,2,3) & City!='Canberra' & RCP %in% c('RCP4.5','RCP8.5') & ehfa>0, # automatically only include EHFa>0
                             .(`DMT~95~ a`=round.fn(`dmt95a`,2),
                               `HW days per year a`=round.fn(sum(ehfa>0)/30,1),
                               `EHF~50p~ a`=round.fn(median(ehfa),2),
                               `EHF~85p~ a`=round.fn(quantile(ehfa, 0.85),2)),
                             by=c('City','Model','RCP','Period')])
dmt <- merge(fdmt, merge(fehf, fehfa, by=c('City','Model','RCP','Period')), by=c('City','Model','RCP','Period'))
write.csv(dmt, file=paste0(outcome.exposure.loc.des,'Projected DMT.csv'), na='', row.names=F) # create csv file


## Rename exposures to match daily.ds
ccia.newnames <- c('ave_sh','Excess heat factor','Excess heat factor (forward)','Excess heat index factor','Excess heat factor index (forward)')
if(adapt==1) { # use EHF metrics with adaptation
  setnames(ccia_future2, old=c('sh','ehfa','ehffa','ehfa.hi','ehffa.hi'), new=ccia.newnames, skip_absent=T) 
} else if(adapt==0.5) { # use EHF metrics with partial adaptation
  setnames(ccia_future2, old=c('sh','ehfpa','ehffpa','ehfpa.hi','ehffpa.hi'), new=ccia.newnames, skip_absent=T) 
} else { # use EHF metrics without adaptation
  setnames(ccia_future2, old=c('sh','ehf','ehff','ehf.hi','ehfa.hi'), new=ccia.newnames, skip_absent=T) 
}


# Stage 3 folders
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'S3 Cen',cenpen,' ',mmpred,'Adapt',adapt))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s3 <- paste0(outcome.exposure.loc,'S3 Cen',cenpen,' ',mmpred,'Adapt',adapt,'/') # folder destination (add /)



########################################################################################################
# Loop components / dimensions
########################################################################################################

gcm <- c('ACCESS1-0','CanESM2','CESM1-CAM5','CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','MIROC5','NorESM1-M') # General circulation models 
# rcp <-  c('rcp26','rcp45','rcp60','rcp85') # not all have rcp26 or rcp60, and it does change their precision (improve, as less models, but uncertainty not captured)
rcp <-  c('RCP4.5','RCP8.5') # not all have rcp26 or rcp60
# if(adapt==0.5) { # Projection periods
#   prd <- '2036-2065' # Partial adaptation is identical to adaptation for 2016-2045
# } else {
prd <- c('2016-2045','2036-2065') 
# } # having only one period will cause apply with an array to lose that dimension, causing error when creating anabs
absrel <- c("abs","rel") # Absolute and relative (difference from future compared to baseline) ANs



########################################################################################################
# Projected outcome data
########################################################################################################

oiproj30 <- oiproj50 <- oidoy <- oidom <- list() # list objects

tot.city <- c(ds.city, 'Total') # all strata plus total
oiperiod <- rep(NA, length(tot.city)) # Average number of outcomes across period
names(oiperiod) <- tot.city # names of locations

for(a in ds.city) {
  .ds <- daily.ds[City==a & !(Month==2 & Day==29)] # dataset per stratum. 29th Feb is not in projected data, so remove this day
  .ds[,doy := leap.ever.year.fn(Date)] # day per year. in case you wish to include 29 Feb, this function will
  oidom[[a]] <- .ds[Month %in% c(10,11,12,1,2,3), lapply(.SD, mean, na.rm=T), by=Day, .SDcols=outcome.var][,get(outcome.var)] # averaged estimate per day of the month
  oidoy[[a]] <- .ds[Month %in% c(10,11,12,1,2,3), lapply(.SD, mean, na.rm=T), by=doy, .SDcols=outcome.var][,get(outcome.var)] # averaged estimate per day of the year
  oiperiod[a] <- sum(oidoy[[a]]) # sum of outcome across projection period averaged annually (minus day 366) per city
  oiproj30[[a]] <- rep(oidoy[[a]], length=length(unique(ccia_future2[City=='Adelaide' & Model==gcm[1] & RCP==rcp[1] & Period==prd[1] & Month %in% c(10,11,12,1,2,3),Date]))) # repeat along 2016-45. Choice of city, gcm, rcp and prd doesn't matter
  oiproj50[[a]] <- rep(oidoy[[a]], length=length(unique(ccia_future2[City=='Adelaide' & Model==gcm[1] & RCP==rcp[1] & Period==prd[2] & Month %in% c(10,11,12,1,2,3),Date]))) # repeat along 2036-65
}

oiperiod[tot.city[length(tot.city)]] <- sum(oiperiod[-length(tot.city)]) # add total to oiperiod
# sum(oiperiod[1:7]) # adds correctly



########################################################################################################
# (1) Extrpaolate exposure-response curve, (2) project and quantify ANNUAL impact and (3) ensemble estimates & quantification of uncertainty
########################################################################################################

# Define objects
ansim <- past_ansim <- array(NA, dim=c(length(prd), length(sim.names), length(absrel), length(gcm), length(rcp), length(tot.city), nsim2+1), 
                             dimnames=list(prd, sim.names, absrel, gcm, rcp, tot.city, c("est",paste0("sim",seq(nsim2))))) # Array to store AN + a dummy for past results

proj.exp.size <- nrow(ccia_future2[City=='Adelaide' & RCP==rcp[1] & Model==gcm[1] & Period==prd[1] & Month %in% c(10,11,12,1,2,3)]) # Choice of city, gcm, rcp and prd doesn't matte
projexpmark <- matrix(NA, proj.exp.size, length(sim.names))
colnames(projexpmark) <- sim.names


gc()
time <- proc.time()[3]
for(a in ds.city) { # loop per city
  for (h in prd) { # loop per period
    for (i in rcp) { # loop per RCP
      for(j in gcm) { # loop per GCM
        .projexposure <- ccia_future2[City==a & RCP==i & Model==j & Period==h & Month %in% c(10,11,12,1,2,3), get(exposure.var)] # projected exposure
        
        if(length(.projexposure)==0) { # skip when no RCP/Model combination and print the combination
          print(paste('Skip:', a, i, j)); next
        }  
        
        # (1) Extrapolation of curve: derive centred basis using projected exposure series and extract parameters
        .li_argvarproj <- list(fun=espline, knots=quantile(.projexposure, eknots, na.rm=T), Bound=range(.projexposure, na.rm=T)) # need Bound for projections, specifcally do.call(onebasis,...)
        
        .cenvec <- do.call(onebasis,c(list(x=centre),.li_argvarproj))
        .projbvar <- do.call(onebasis,c(list(x=.projexposure),.li_argvarproj))
        .bvarcen <- scale(.projbvar,center=.cenvec,scale=F)
        
        # Indicators to meet heatwave thresholds
        projexpmark[,1] <- .projexposure<=exposure_ehf[a,1] # not heatwave
        projexpmark[,2] <- .projexposure>exposure_ehf[a,1] # heatwave
        projexpmark[,3] <- .projexposure>exposure_ehf[a,1] & .projexposure<exposure_ehf[a,2] # low-intensity
        projexpmark[,4] <- .projexposure>exposure_ehf[a,2] & .projexposure<exposure_ehf[a,3] # severe
        projexpmark[,5] <- .projexposure>exposure_ehf[a,3] # extreme (>=2 * EHF85)
        
        
        # (2) Impact projections
        if(h == prd[1]) {.oiproj <- oiproj30[[a]]} # set oiproj based on period
        else if(h == prd[2]) {.oiproj <- oiproj50[[a]]}
        
        .blups <- blups[[which(ds.city==a)]] # store BLUP
        # an <- (1-exp(-.bvarcen%*%coef[a,]))*.oiproj # compute daily contributions of attributable outcomes, pre-BLUP
        an <- (1-exp(-.bvarcen%*%.blups$blup))*.oiproj # compute daily contributions of attributable outcomes, BLUP, for each day in projected period
        
        # Store sum of AN (attributable numbers) by exposure range and period. Store in array before iterations. Do 2030 and 2050 separately
        for(t in 1:5) { # attributable categories: below centre, above centre, low-intensity heatwave, severe heatwave, extreme heatwave
          ansim[h,sim.names[t],"abs",j,i,a,1] <- sum(an[projexpmark[,t]], na.rm=T)/30 # divide by 30 for annual AN
        } 
        
        # (3) Estimate uncertainty of projected AN
        set.seed(1 + which(gcm==j) + which(rcp==i)*10 + which(prd==h)*100)
        # .coefsim <- mvrnorm(nsim2,coef[a,],vcov[[a]]) # sample coef assuming multivariate normal distribution, pre-BLUP
        .coefsim <- mvrnorm(nsim2,.blups$blup,.blups$vcov) # sample coef assuming multivariate normal distribution, BLUP
        
        for(s in seq(nsim2)) { # loop across iterations
          .an <- (1-exp(-.bvarcen%*%.coefsim[s,]))*.oiproj # compute daily contributions of attributable outcomes for each day in projected period
          for(t in 1:5) { # Store sum of AN as above, but for each simulation
            ansim[h,sim.names[t],"abs",j,i,a,s+1] <- sum(.an[projexpmark[,t]], na.rm=T)/30
          }
        }
      }
    }
  }
}
proc.time()[3]-time # about 1 min on Mac


## Add national results to array, filling in the Total/National category
ansim[,,,,,tot.city[length(tot.city)],] <- apply(ansim, c(1:5,7), sum, na.rm=T) # sum by all dimensions but city (sums city results)


## Estimate AN relative "rel" to baseline (absolute/projected - baseline) for all cities and total, assuming no change in pop
for(a in tot.city) { # fill in dummy ansim with past results, including total/national
  # past_ansim[,,,,,a,] <- matrix(rep(rbind(matsim, 'Total'=antot)[a,],2),nrow=2, byrow=T) # duplicate for both rows (periods)
  past_ansim[,,,,,a,] <- matrix(rep(rbind(matsim.yr, 'Total'=antot.yr)[a,],2),nrow=2, byrow=T) # duplicate for both rows (periods), annual results
}
ansim[,,"rel",,,,] <- ansim[,,"abs",,,,] - past_ansim[,,"abs",,,,] # calculate rel. Limit to one dimension for rel/ans



########################################################################################################
# Summarize results: compute AN/AF with 95% CI in ensemble by range, period and RCP
########################################################################################################

estci <- c("est","ci.l","ci.u") # store estimate, lci and uci

# Create 4 new arrays to store ensemble estimates (average impacts across GCMs) with empirical 95% CIs
# 2 arrays store abs and rel AN, another 2 to store estimated rel and abs AF
# 5 array dimensions: period, point estimate and CI, attributable category, scenario, and city
anabs <- afabs <- anrel <- afrel <- array(
  NA, dim=c(length(prd), length(estci), length(sim.names), length(rcp), length(tot.city)),
  dimnames=list(prd, estci, sim.names, rcp, tot.city))

# Store attributable numbers. GCM dimension in ansim not looped, as results calculated across all GCMs
for(a in tot.city) { # by city including total
  for (i in rcp) { # loop per RCP
    anabs[,"est",,i,a] <- apply(ansim[,,"abs",,i,a,1], 1:2, mean, na.rm=T) # Absolute AN
    anabs[,"ci.l",,i,a] <- apply(ansim[,,"abs",,i,a,-1], 1:2, quantile, 0.025, na.rm=T)
    anabs[,"ci.u",,i,a] <- apply(ansim[,,"abs",,i,a,-1], 1:2, quantile, 0.975, na.rm=T)
    
    anrel[,"est",,i,a] <- apply(ansim[,,"rel",,i,a,1], 1:2, mean, na.rm=T) # Relative AN
    anrel[,"ci.l",,i,a] <- apply(ansim[,,"rel",,i,a,-1], 1:2, quantile, 0.025, na.rm=T)
    anrel[,"ci.u",,i,a] <- apply(ansim[,,"rel",,i,a,-1], 1:2, quantile, 0.975, na.rm=T)
  }
  afabs[,,,,a] <- anabs[,,,,a]/oiperiod[a]*100 # Absolute AF
  afrel[,,,,a] <- anrel[,,,,a]/oiperiod[a]*100 # Relative AF
} # na.rm=T only makes a difference if using RCP2.6 and 6.0, which do not have data for all GCMs



########################################################################################################
# Format AF
########################################################################################################

## Combine AF results 
projafa <- as.data.table(afabs)
projafr <- as.data.table(afrel)
projafa[,type:='Absolute']
projafr[,type:='Relative']
projaf <- rbindlist(list(projafa, projafr))
colnames(projaf) <- c('Period','Estimate','Category','RCP','City','Result','Value')

## Reshape AF results to include estimates, LCI and UCI on one row
projafby <- c('Period','Category','RCP','City','Value')
projaf1 <- merge(merge(projaf[Estimate=='est',-"Estimate"], projaf[Estimate=='ci.l',-"Estimate"], by=projafby), projaf[Estimate=='ci.u',-"Estimate"], by=projafby)
colnames(projaf1) <- c('Period','Category','RCP','City','Value','Estimate','LCI','UCI')
projaf1[,AF:=paste0(round.fn(Estimate, af.round),' (',round.fn(LCI, af.round),' to ',round.fn(UCI, af.round),')')]
projaf1[,':='(Estimate=NULL, LCI=NULL, UCI=NULL)] # do not want these

## Reshape AF results to include outcomes on one row
projaf1by <- c('Period','RCP','City','Value')
projaf2 <- merge(merge(merge(merge(projaf1[Category==sim.names[1],-"Category"],
                                   projaf1[Category==sim.names[2],-"Category"], by=projaf1by, suffixes=c('nhw','hw')),
                             projaf1[Category==sim.names[3],-"Category"], by=projaf1by),
                       projaf1[Category==sim.names[4],-"Category"], by=projaf1by, suffixes=c('l','s')),
                 projaf1[Category==sim.names[5],-"Category"], by=projaf1by)
colnames(projaf2) <- c('Period','RCP','City','Value',sim.names)

## Save AF
save(projaf,projaf1,projaf2, file=paste0(outcome.exposure.loc.s3, 'projaf.rda')) # save AF results
write.csv(projaf2, file=paste0(outcome.exposure.loc.s3, 'Fractions.csv'), na='', row.names=F) # create csv file



########################################################################################################
# Format AN
########################################################################################################

## Combine AN results 
projana <- as.data.table(anabs)
projanr <- as.data.table(anrel)
projana[,type:='Absolute']
projanr[,type:='Relative']
projan <- rbindlist(list(projana, projanr)) # combine absolute, relative, and projections based on population
colnames(projan) <- c('Period','Estimate','Category','RCP','City','Result','Value')

## Attributable numbers considering change in population
## "RELATIVE" VALUE (CHANGE RELATIVE TO BASELINE) WILL NOT LONGER BE ACCURATE WHEN APPLYING SEIRES (NO CHANGE STILL OK)
projans <- merge(projan,projpop2, by=c('City','Period')) # add population projections to AN
projans[,':='(VA=Result*popratioA, VB=Result*popratioB, VC=Result*popratioC)] # create values for each series
setnames(projans, 'Result','V0',skip_absent=T) # to more clearily separate from other values
projans1 <- melt(projans, id.vars=c('City','Period','RCP','Value','Estimate','Category'), measure.vars=c('V0','VA','VB','VC'), variable.name='Series', value.name='Result') # melt by series
projans1[str_detect(Series,'0'), Series:='No change'] 
projans1[str_detect(Series,'A'), Series:='A']
projans1[str_detect(Series,'B'), Series:='B']
projans1[str_detect(Series,'C'), Series:='C']

## Reshape AN results to include estimates, LCI and UCI on one row
projanby <- c(projafby,'Series')
projan1 <- merge(merge(projans1[Estimate=='est',-"Estimate"], projans1[Estimate=='ci.l',-"Estimate"], by=projanby), projans1[Estimate=='ci.u',-"Estimate"], by=projanby)
colnames(projan1) <- c(projanby,'Estimate','LCI','UCI')
projan1[,AN:=paste0(round.fn(Estimate,afyr.r),' (',round.fn(LCI,afyr.r),' to ',round.fn(UCI,afyr.r),')')]
projan1[,':='(Estimate=NULL, LCI=NULL, UCI=NULL)] # do not want these

## Reshape AN results to include outcomes on one row
projan1by <- c(projaf1by,'Series')
projan2 <- merge(merge(merge(merge(projan1[Category==sim.names[1],-"Category"],
                                   projan1[Category==sim.names[2],-"Category"], by=projan1by, suffixes=c('nhw','hw')),
                             projan1[Category==sim.names[3],-"Category"], by=projan1by),
                       projan1[Category==sim.names[4],-"Category"], by=projan1by, suffixes=c('l','s')),
                 projan1[Category==sim.names[5],-"Category"], by=projan1by)
colnames(projan2) <- c(projan1by, sim.names)

## Save AN
save(projan,projan1,projan2, file=paste0(outcome.exposure.loc.s3, 'projan.rda')) # save AN result
write.csv(projan2, file=paste0(outcome.exposure.loc.s3, 'Numbers.csv'), na='', row.names=F) # create csv file



################################################################################
# Checks
################################################################################

# View(projan2[Series=='No change'])
# View(projan2[Series=='B'])
# projaf2

# results definitely a little odd from rcp2.6 and particularly rcp6.0, likely from lack of GCMs. These should be treated with caution or removed
# ccia_future2[RCP %in% c('rcp45','rcp85'), lapply(.SD, mean, na.rm=T), by=.(City,RCP,Period), .SDcols=c('tas','tasmax','Excess heat factor','ehisig','ehiaccl')] # do increase slightly from rcp45 to rp85
# despite this, EHF AF is slightly lower for rcp85 compared to rcp45, despite the generally higher temperatures and ?more heatwave conditions. Perhaps acclimatisation was improved, with dmt3 being more comparable more often than before to last 30 days



######################### END ############################
######################### END ############################
  