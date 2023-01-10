##################################################
### Load packages
##################################################

# The versions for each package used in testing are included as comments next to their name

# Data manipulation
library(readxl) # Read Excel files (read_excel)
library(data.table) # data.table manipulation
library(lubridate) # Date manipulation
library(stringr) # String commands
library(zoo) # rollapply()

# Statistical analysis
library(statmod) # Tweedie distribution
library(tweedie) # Estimate Tweedie shape parameter
library(mgcv) # Generalized additive models
library(dlnm) # Distributed lag non-linear models
library(splines) # Splines
library(mixmeta) # Multivariate meta-analysis
library(FluMoDL) # Compute attributable risk
library(weathermetrics) # Compute attributable risk

# Graphing
library(pryr) # %<a-%



##################################################
### Prepare meteorological data and save as .rda file
###   This step is already done for you and is commented out
###   However, the code is included to demonstrate how the original data was extracted and edited
##################################################

# # Source: Brambilla et al. 2022, Data in Brief, https://doi.org/10.1016/j.dib.2022.108291
# 
# ## Extract data

brambilla.loc <- paste0("your_directory_with_the_files") # Location of Brambilla et al. 2022 Appendix 1 Climate_files
brambilla.loc <- paste0("/Users/MatthewBorg/Library/CloudStorage/Box-Box/Data original/Brambilla 2022/Appendix 1/Climate_files/") # Location of Brambilla et al. 2022 Appendix 1 Climate_files
brambilla.data <- list() # Save data in list
for(a in c('Melbourne','Sydney')) {
  for(b in 21:25) { # 2021-2015
    if (a=='Sydney') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'5_Sydney_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
    else if (a=='Melbourne') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'6_Melbourne_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
  }
}

brambilla <- do.call(rbind.data.frame, brambilla.data) # Combine list into a data frame
brambilla <- as.data.table(brambilla) # Save as data.table. This will remove row names
colnames(brambilla) <- c('City','Year','Month','Day','Hour','temp','rh','ws','ap','dnr','dhr')
brambilla[,':='(ws=NULL,ap=NULL,dnr=NULL,dhr=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT

## Calculate het index
source('heat.index2.r') # weathermetrics::heat.index modified with a 79 threshold changed to 80 and no rounding
brambilla[,hi := heat.index2(t=temp, rh=rh, temperature.metric='celsius', round=9999999),] # Heat index. Use self-code that allows rounding to >2 digits

## Obtain daily metrics at time of maximum temperature using BoM classification
brambilla[,predate:=ISOdate(Year, Month, Day, Hour)] # Date
brambilla[,Date:=as.Date(predate-9*3600+1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on same day to 9am next day
brambilla.max <- setDT(brambilla)[Date!='2010-12-31', .SD[which.max(temp)], by=c('Date','City')] # Metrics at time of maximum temperature
brambilla[,Date:=as.Date(predate+15*3600-1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on previous day to 9am same day
brambilla.min <- setDT(brambilla)[Date!='2010-12-31', .SD[which.max(temp)], by=c('Date','City')] # Metrics at time of minimum temperature
brambilla.max[,':='(Year=NULL,Month=NULL,Day=NULL,Hour=NULL,predate=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT
brambilla.min[,':='(Year=NULL,Month=NULL,Day=NULL,Hour=NULL,predate=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT

### Combine max and min
brambilla.all <- merge(brambilla.max, brambilla.min, by=c('Date','City'), all=T, suffixes = c(".max",".min"))


### Save meteorolgical dataset
save(brambilla.all, file='brambilla.all.rda')
load('brambilla.all.rda')



##################################################
### EHF
################################################## 

### EHF variables
## Use City == shift(City, -x) to ensure that same city is used, but this misses x rows based on date and lag (start, Adelaide) /lead (end, Sydney)
brambilla.all[City == shift(City, -1) ,  ':='(temp.min_lead=shift(temp.min,-1), min_hi_lead=shift(hi.min,-1))] # EHF uses min T from same day
brambilla.all[City == 'Sydney' & Date>='2019-02-26', ':='(temp.min_lead=shift(temp.min,-1), min_hi_lead=shift(hi.min,-1))] # above misses Sydney 27/2/19, the 2nd to last day, which this addresses
brambilla.all <- brambilla.all[Date < '2019-02-28',] # now that EHF is calculated, remove '2019-02-28' as it lacks enough data to be calculated for max T
brambilla.all[, temp_ehf := rowMeans(.SD), .SDcols = c("temp.max","temp.min_lead")]
brambilla.all[, hi_ehf := rowMeans(.SD), .SDcols = c("hi.max","min_hi_lead")]

for(i in unique(brambilla.all$City)) {
  brambilla.all[City==i, dmt3 := rollmean(temp_ehf, 3, fill=NA, align='right')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
  brambilla.all[City==i, dmtp30 := (shift(temp_ehf,1+2) + shift(temp_ehf,2+2) + shift(temp_ehf,3+2) + shift(temp_ehf,4+2) + shift(temp_ehf,5+2)
                                + shift(temp_ehf,6+2) + shift(temp_ehf,7+2) + shift(temp_ehf,8+2) + shift(temp_ehf,9+2) + shift(temp_ehf,10+2)
                                + shift(temp_ehf,11+2) + shift(temp_ehf,12+2) + shift(temp_ehf,13+2) + shift(temp_ehf,14+2) + shift(temp_ehf,15+2)
                                + shift(temp_ehf,16+2) + shift(temp_ehf,17+2) + shift(temp_ehf,18+2) + shift(temp_ehf,19+2) + shift(temp_ehf,20+2)
                                + shift(temp_ehf,21+2) + shift(temp_ehf,22+2) + shift(temp_ehf,23+2) + shift(temp_ehf,24+2) + shift(temp_ehf,25+2)
                                + shift(temp_ehf,26+2) + shift(temp_ehf,27+2) + shift(temp_ehf,28+2) + shift(temp_ehf,29+2) + shift(temp_ehf,30+2))/30] # retrospective (with +2)
  brambilla.all[City==i, dmt3f := rollmean(temp_ehf, 3, fill=NA, align='left')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
  brambilla.all[City==i, dmtp30f := (shift(temp_ehf,1) + shift(temp_ehf,2) + shift(temp_ehf,3) + shift(temp_ehf,4) + shift(temp_ehf,5)
                                 + shift(temp_ehf,6) + shift(temp_ehf,7) + shift(temp_ehf,8) + shift(temp_ehf,9) + shift(temp_ehf,10)
                                 + shift(temp_ehf,11) + shift(temp_ehf,12) + shift(temp_ehf,13) + shift(temp_ehf,14) + shift(temp_ehf,15)
                                 + shift(temp_ehf,16) + shift(temp_ehf,17) + shift(temp_ehf,18) + shift(temp_ehf,19) + shift(temp_ehf,20)
                                 + shift(temp_ehf,21) + shift(temp_ehf,22) + shift(temp_ehf,23) + shift(temp_ehf,24) + shift(temp_ehf,25)
                                 + shift(temp_ehf,26) + shift(temp_ehf,27) + shift(temp_ehf,28) + shift(temp_ehf,29) + shift(temp_ehf,30))/30] # retrospective (with +2)
  brambilla.all[City==i, dmhi3 := rollmean(hi_ehf, 3, fill=NA, align='right')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
  brambilla.all[City==i, dmhip30 := (shift(hi_ehf,1+2) + shift(hi_ehf,2+2) + shift(hi_ehf,3+2) + shift(hi_ehf,4+2) + shift(hi_ehf,5+2)
                                 + shift(hi_ehf,6+2) + shift(hi_ehf,7+2) + shift(hi_ehf,8+2) + shift(hi_ehf,9+2) + shift(hi_ehf,10+2)
                                 + shift(hi_ehf,11+2) + shift(hi_ehf,12+2) + shift(hi_ehf,13+2) + shift(hi_ehf,14+2) + shift(hi_ehf,15+2)
                                 + shift(hi_ehf,16+2) + shift(hi_ehf,17+2) + shift(hi_ehf,18+2) + shift(hi_ehf,19+2) + shift(hi_ehf,20+2)
                                 + shift(hi_ehf,21+2) + shift(hi_ehf,22+2) + shift(hi_ehf,23+2) + shift(hi_ehf,24+2) + shift(hi_ehf,25+2)
                                 + shift(hi_ehf,26+2) + shift(hi_ehf,27+2) + shift(hi_ehf,28+2) + shift(hi_ehf,29+2) + shift(hi_ehf,30+2))/30] # retrospective (with +2)
  brambilla.all[City==i, dmhi3f := rollmean(hi_ehf, 3, fill=NA, align='left')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
  brambilla.all[City==i, dmhip30f := (shift(hi_ehf,1) + shift(hi_ehf,2) + shift(hi_ehf,3) + shift(hi_ehf,4) + shift(hi_ehf,5)
                                  + shift(hi_ehf,6) + shift(hi_ehf,7) + shift(hi_ehf,8) + shift(hi_ehf,9) + shift(hi_ehf,10)
                                  + shift(hi_ehf,11) + shift(hi_ehf,12) + shift(hi_ehf,13) + shift(hi_ehf,14) + shift(hi_ehf,15)
                                  + shift(hi_ehf,16) + shift(hi_ehf,17) + shift(hi_ehf,18) + shift(hi_ehf,19) + shift(hi_ehf,20)
                                  + shift(hi_ehf,21) + shift(hi_ehf,22) + shift(hi_ehf,23) + shift(hi_ehf,24) + shift(hi_ehf,25)
                                  + shift(hi_ehf,26) + shift(hi_ehf,27) + shift(hi_ehf,28) + shift(hi_ehf,29) + shift(hi_ehf,30))/30] # retrospective (with +2)
}
brambilla.all[, dmt95 := quantile(temp_ehf, 0.95, na.rm=T), by=City] # 95th quantile of dmt (average air t), by city
brambilla.all[, dmhi95 := quantile(hi_ehf, 0.95, na.rm=T), by=City] # 95th quantile of dmhi (average hi), by city
brambilla.all[, ehisig := dmt3 - dmt95] # EHIsig
brambilla.all[, ehisigf := dmt3f - dmt95] # EHIsig forward
brambilla.all[, ehiaccl := dmt3 - dmtp30] # EHIaccl
brambilla.all[, ehiacclf := dmt3f - dmtp30f] # EHIaccl forward
brambilla.all[, ehisig.hi := dmhi3 - dmhi95] # EHIsig HI
brambilla.all[, ehisigf.hi := dmhi3f - dmhi95] # EHIsig forward HI
brambilla.all[, ehiaccl.hi := dmhi3 - dmhip30] # EHIaccl HI
brambilla.all[, ehiacclf.hi := dmhi3f - dmhip30f] # EHIaccl forward HI
brambilla.all[, ehf := ehisig * pmax(1, ehiaccl)] # EHF
brambilla.all[, ehff := ehisigf * pmax(1, ehiacclf)] # EHF forward
brambilla.all[, ehf.hi := ehisig.hi * pmax(1, ehiaccl.hi)] # EHF HI
brambilla.all[, ehff.hi := ehisigf.hi * pmax(1, ehiacclf.hi)] # EHF forward EHF
# table(is.na(brambilla.all$temp.min_lead), useNA='ifany') # 0 as expected after date correction
# table(is.na(brambilla.all$temp_ehf), useNA='ifany') # 0 as expected
# table(is.na(brambilla.all$dmt3), useNA='ifany') # 16 as expected after date correction
# table(is.na(brambilla.all$dmtp30), useNA='ifany') # 256 as expected
# table(is.na(brambilla.all$dmt95), useNA='ifany') # none as expected
# table(is.na(brambilla.all$ehisig), useNA='ifany') # 16 as expected
# table(is.na(brambilla.all$ehiaccl), useNA='ifany') # 256 as expected
# table(is.na(brambilla.all$ehf), useNA='ifany') # 256 as expected
# table(brambilla.all$ehf >= 0, useNA = 'ifany') # 2980 / (nrow(brambilla.all)-256) = 0.03508235, looks right
# table(is.na(brambilla.all$ehff.hi), useNA='ifany') # 256 as expected
# table(brambilla.all$ehff.hi >= 0, useNA = 'ifany') # 2980 / (nrow(brambilla.all)-256) = 0.03508235, looks right

# Save EHF DMT95 reference values
brambilla_ehfr <- unique(brambilla.all[, .(City,dmt95,dmhi95)])
save(brambilla_ehfr, file='Results/brambilla_ehfr.rda')


### HW variables using EHF.  If any day in 3-day period meet criteria, than all days are hw

# HW 
brambilla.all[, hw := fifelse(ehf > 0 | shift(ehf > 0) | shift(ehf > 0, 2), 1, 0)]
brambilla.all[is.na(hw), hw := fifelse(ehf > 0 | shift(ehf > 0), 1, 0)] # cover missing NA for first lag
brambilla.all[is.na(hw), hw := fifelse(ehf > 0, 1, 0)]# cover missing NA for same day

# Severe heatwave
# for(i in unique(brambilla.all$City)) {
#   print(i)
#   print(summary(brambilla.all[City==i, ehf]))
#   print(summary(brambilla.all[City==i & ehf > 0, ehf]))
# }
brambilla.all[ehf > 0, ehf85 := quantile(ehf, 0.85, na.rm=F), by=City] # 85th quantile of all positive EHF values, by city
brambilla.all[ehff > 0, ehff85 := quantile(ehff, 0.85, na.rm=F), by=City] # 85th quantile of all positive EHF values, by city
brambilla.all[ehf.hi > 0, ehf.hi85 := quantile(ehf.hi, 0.85, na.rm=F), by=City] # 85th quantile of all positive EHF values, by city
brambilla.all[ehff.hi > 0, ehff.hi85 := quantile(ehff.hi, 0.85, na.rm=F), by=City] # 85th quantile of all positive EHF values, by city

brambilla.all[, hw1 := fifelse(ehf > ehf85 | shift(ehf > ehf85) | shift(ehf > ehf85, 2), 1, 0)]
brambilla.all[is.na(hw1), hw1 := fifelse(ehf > ehf85 | shift(ehf > ehf85), 1, 0)] # cover missing NA for first lag
brambilla.all[is.na(hw1), hw1 := fifelse(ehf > ehf85, 1, 0)] # cover missing NA for same day

# Extreme heatwave 
brambilla.all[, hw2 := fifelse(ehf > ehf85*2 | shift(ehf > ehf85*2) | shift(ehf > ehf85*2, 2), 1, 0)]
brambilla.all[is.na(hw2), hw2 := fifelse(ehf > ehf85*2 | shift(ehf > ehf85*2), 1, 0)] # cover missing NA for first lag
brambilla.all[is.na(hw2), hw2 := fifelse(ehf > ehf85*2, 1, 0)] # cover missing NA for same day
# table(is.na(brambilla.all$ehf85), useNA = 'ifany') # 0 as expected
# table(brambilla.all$hw, useNA = 'ifany') # 4962
# table(brambilla.all$hw1, useNA = 'ifany') # 907
# table(brambilla.all$hw2, useNA = 'ifany') # 205
# head(brambilla.all[, .(date, City, dmt3, dmt95, ehisig, ehiaccl, ehf, hw, ehf85, hw1, hw2)],50) # looks fine
# View(brambilla.all[is.na(hw1), .(date, City, dmt3, ehisig, ehiaccl, ehf, hw, ehf85, hw1, hw2)])
# brambilla.all[City=='Darwin' & date <= '1990-01-05', .(date, City, dmt3, ehisig, ehiaccl, ehf, hw, ehf85, hw1, hw2)]


### Save
save(brambilla.all, file=paste0(barra.loc,'brambilla.all.rda'))


### EHF apparent T
# ehf.odappt <- with(daily.bom1, ehf.fn(ave.od.app.t, Date, station.no, sdate.="2005-07-01", edate.="2018-07-01")) # Outdoor apparent T
# ehf.shappt <- with(daily.bom1, ehf.fn(ave.sh.app.t, Date, station.no, sdate.="2005-07-01", edate.="2018-07-01")) # Shade apparent T
# ehf.idappt <- with(daily.bom1, ehf.fn(ave.id.app.t, Date, station.no, sdate.="2005-07-01", edate.="2018-07-01")) # Indoor apparent T
# ehf.vars.keep <- c('station.no', 'Date', 'ehf', 'hw', 'hw1', 'hw2', 'hw3')
# ehf.odappt <- rename(ehf.odappt[,ehf.vars.keep], odappt.ehf=ehf, odappt.hw=hw, odappt.hw1=hw1, odappt.hw2=hw2, odappt.hw3=hw3)
# ehf.shappt <- rename(ehf.shappt[,ehf.vars.keep], shappt.ehf=ehf, shappt.hw=hw, shappt.hw1=hw1, shappt.hw2=hw2, shappt.hw3=hw3)
# ehf.idappt <- rename(ehf.idappt[,ehf.vars.keep], idappt.ehf=ehf, idappt.hw=hw, idappt.hw1=hw1, idappt.hw2=hw2, idappt.hw3=hw3)
# 
# daily.bom2 <- left_join(daily.bom2, ehf.meant[,ehf.vars.keep], by=c('station.no','Date'))
# daily.bom2 <- left_join(daily.bom2, ehf.odappt, by=c('station.no','Date'))
# daily.bom2 <- left_join(daily.bom2, ehf.shappt, by=c('station.no','Date'))
# daily.bom2 <- left_join(daily.bom2, ehf.idappt, by=c('station.no','Date'))
# rm(ehf.meant,ehf.odappt,ehf.shappt,ehf.idappt) # save room



##################################################
### Create simulated claims data
##################################################

## Create randomly-generated claims data with Poisson and Tweedie assumed distributions for number of OIIs (claims) and total costs, respectively
## This simulated data has no relationship to any predictor variables
no.models <- 4 # 2 cities (Melbourne and Sydney) and their indoor/outdoor combinations
set.seed(7)
dates <- rep(seq(as_date("2011-01-01"), as_date("2015-12-31"), 1),each=no.models)
ldates <- length(dates)/no.models
claims <- data.table(
  Date = dates,
  City = rep(c('Melbourne','Melbourne','Sydney','Sydney'),ldates),
  outin = rep(c('Indoors','Outdoors'),ldates*2),
  'Number of OIIs' = c(rpois(ldates,60), rpois(ldates,22), rpois(ldates,145), rpois(ldates,50)), # Number of OIIs
  'Total costs' = c(rtweedie(ldates, xi=1.5, mu=240, phi=8), rtweedie(ldates, xi=1.7, mu=85, phi=9),
                    rtweedie(ldates, xi=1.4, mu=515, phi=11), rtweedie(ldates, xi=1.6, mu=200, phi=7)) # Total costs
)
claims[,'Costs per OII':=get('Total costs')/get('Number of OIIs')] # Costs per OII, should you wish to model this

# By variables for commands
by.vars <- c("Date","City","outin")
by.vars2 <-by.vars[!by.vars == 'outin']
by.vars3 <- by.vars[!by.vars == 'Date']



##################################################
### Make public holiday data long
##################################################

load(file='Public holidays.rda') 
public.holidays <- melt(public.hols, id.vars=c('Date','Public holiday'),
                        measure.vars=c('Adelaide','Brisbane','Canberra','Darwin','Hobart','Melbourne','Perth','Sydney'),
                        variable.name='City', value.name='phol')
public.holidays[, Date:=as.Date(Date)]



##################################################
### Merge data sets together
##################################################

load(file='brambilla.max.rda') # Climate data
load(file='school.holidays.rda')
load(file='pop.rda') # Worker's population using ABS data. Due to estimating indoor/outdoor proportions, n can have decimal places


daily <- merge(merge(merge(claims, brambilla.max, by=c('Date','City')),
               public.holidays[!is.na(phol)], by=c('Date','City'), all.x=T),
                school.holidays, by=c('Date','City'), all.x=T)
daily[is.na(phol), phol := 0] # Replace NA in public holidays with 0

# Merge worker's month.y population
daily.ds <- merge(daily, pop, by=c("Year","Month","City","outin"), all.x=T)

# Create WBGT based on indoors or outdoors
daily.ds[, 'Maximum WBGT':=fifelse(outin=='Outdoors',wbgtl,wbgtb)]

# Stratum variables
daily.ds[,stratum := do.call(paste, c(mget(by.vars3), sep=' '))] # combine by.vars3 into a string as stratum
ds.stratum <- sort(unique(daily.ds[,stratum]))
length.ds.stratum <- length(ds.stratum)



##################################################
### Time variables
##################################################

## Financial year (July to June in Australia)
daily.ds[, FYear := fifelse(Month <= 6, Year-1, Year)]

## Month as a factor variable
daily.ds[,month := factor(Month, labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))]

## Day of the week variables (one with all categories, rest are binary variables per day
daily.ds[, dow:=weekdays(Date)]
daily.ds[, Mon:=fifelse(str_detect(dow,'Mon'), 1, 0)]
daily.ds[, Tue:=fifelse(str_detect(dow,'Tue'), 1, 0)]
daily.ds[, Wed:=fifelse(str_detect(dow,'Wed'), 1, 0)]
daily.ds[, Thu:=fifelse(str_detect(dow,'Thu'), 1, 0)]
daily.ds[, Fri:=fifelse(str_detect(dow,'Fri'), 1, 0)]
daily.ds[, Sat:=fifelse(str_detect(dow,'Sat'), 1, 0)]
daily.ds[, Sun:=fifelse(str_detect(dow,'Sun'), 1, 0)]

## Format public holidays into a factor
daily.ds[, public.hol := factor(phol, levels=c(0,1), labels=c('No','Yes'))]

## Special holidays
daily.ds[, shol:=0] # Default value for special holidays
daily.ds[Month==12 & Day %in% c(23:30), shol:=1] # Christmas break
daily.ds[Month==12 & Day==31, shol:=2] # NYE
daily.ds[Month==1 & Day==1, shol:=3] # NYD
daily.ds[Month==1 & Day %in% c(2:4), shol:=4] # 2nd to 4th January
daily.ds[City=='Sydney' & str_detect(`Public holiday`, 'Australia Day'), shol:=5] # Australia Day has a public celebration at the Sydney Opera House
day.before.melbourne.cup <- as.Date(c("2004-11-01","2005-10-31","2006-11-06","2007-11-05","2008-11-03","2009-11-02","2010-11-01","2011-10-31","2012-11-05","2013-11-04","2014-11-03","2015-11-02","2016-10-31","2017-11-06","2018-11-05","2019-11-04","2020-11-02","2021-11-01")) # List of days before Melboure Cup from 2004 to 2021 
daily.ds[City=='Melbourne' & Date %in%  day.before.melbourne.cup, shol:=6] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
daily.ds[, shol := factor(shol, labels=c('No','Xmas period','NYE','NYD','2-4 Jan','Australia Day','Day before Melbourne Cup'))] # Make sphol an unordered factor

## 1st day of the month, excluding New Year's Day
daily.ds[, day1 := as.factor(fifelse(Day==1 & Month!=1, 'Yes','No'))]

## Numeric date for gam()
daily.ds[, date := as.numeric(Date)]



######################### END ############################
######################### END ############################
