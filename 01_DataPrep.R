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
library(weathermetrics) # Modified version of heat index
library(MASS) # mnvnorm(), only used for projections

# Graphing
library(pryr) # %<a-%

rm(list=ls())   # remove existing variables


##################################################
### Prepare retrospective meteorological data and save as .rda file
###   This step is already done for you and is commented out
###   However, the code is included to demonstrate how the original data was extracted and edited
##################################################

# Source: Brambilla et al. 2022, Data in Brief, https://doi.org/10.1016/j.dib.2022.108291

# ## Extract data
#
# brambilla.loc <- paste0("your_directory_with_the_files") # Location of Brambilla et al. 2022 Appendix 1 Climate_files
# brambilla.loc <- paste0("/Users/MatthewBorg/Library/CloudStorage/Box-Box/Data original/Brambilla 2022/Appendix 1/Climate_files/") # Location of Brambilla et al. 2022 Appendix 1 Climate_files
# brambilla.loc <- paste0("C:/Users/a1210385/Box/Data original/Brambilla 2022/Appendix 1/Climate_files/") # Location of Brambilla et al. 2022 Appendix 1 Climate_files
#
# brambilla.data <- list() # Save data in list
# for(a in c('Melbourne','Sydney')) {
#   for(b in 1:30) { # 1991-2020, one b per year. This forms a 30-year period for the EHF reference threshold
#     if (a=='Sydney') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'5_Sydney_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
#     else if (a=='Melbourne') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'6_Melbourne_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
#   }
# }
#
# brambilla <- do.call(rbind.data.frame, brambilla.data) # Combine list into a data frame
# brambilla <- as.data.table(brambilla) # Save as data.table. This will remove row names
# colnames(brambilla) <- c('City','Year','Month','Day','Hour','temp','rh','ws','ap','dnr','dhr')
# brambilla[,':='(ws=NULL,ap=NULL,dnr=NULL,dhr=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT
#
# ## Calculate het index
# source('heat.index2.r') # weathermetrics::heat.index modified with a 79 threshold changed to 80 and no rounding
# brambilla[,hi := heat.index2(t=temp, rh=rh, temperature.metric='celsius', round=9999999),] # Heat index. Use self-code that allows rounding to >2 digits
#
# ## Obtain daily metrics at time of maximum temperature using BoM classification
# brambilla[,predate:=ISOdate(Year, Month, Day, Hour)] # Date
# brambilla[,Date:=as.Date(predate-9*3600+1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on same day to 9am next day
# brambilla.max <- setDT(brambilla)[Date!='1990-12-31', .SD[which.max(temp)], by=c('Date','City')] # Metrics at time of maximum temperature
# brambilla[,Date:=as.Date(predate+15*3600-1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on previous day to 9am same day
# brambilla.min <- setDT(brambilla)[Date!='1990-12-31', .SD[which.min(temp)], by=c('Date','City')] # Metrics at time of minimum temperature
# brambilla.max[,':='(Year=NULL,Month=NULL,Day=NULL,Hour=NULL,predate=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT
# brambilla.min[,':='(Year=NULL,Month=NULL,Day=NULL,Hour=NULL,predate=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT
#
# ### Combine max and min
# brambilla.all <- merge(brambilla.max, brambilla.min, by=c('Date','City'), all=T, suffixes = c(".max",".min"))
# brambilla.all[,rh.ave := (rh.max+rh.min)/2] # Average relative humidity for supplementary analysis
#
# ### EHF variables
# ## Use City == shift(City, -x) to ensure that same city is used, but this misses x rows based on date and lag (start, Adelaide) /lead (end, Sydney)
# setorder(brambilla.all, City, Date)
# brambilla.all[City == shift(City, -1) , ':='(temp.min_lead=shift(temp.min,-1), min_hi_lead=shift(hi.min,-1))] # EHF uses min T from same day
# brambilla.all <- brambilla.all[Date < '2020-12-31',] # now that min_lead is calculated, remove the last day of the study period, as it lacks enough data to be calculated for max T
# brambilla.all[, temp_ehf := rowMeans(.SD), .SDcols = c("temp.max","temp.min_lead")]
# brambilla.all[, hi_ehf := rowMeans(.SD), .SDcols = c("hi.max","min_hi_lead")]
#
# for(i in unique(brambilla.all$City)) {
#   brambilla.all[City==i, dmt3 := rollmean(temp_ehf, 3, fill=NA, align='right')]
#   brambilla.all[City==i, dmtp30 := (shift(temp_ehf,1+2) + shift(temp_ehf,2+2) + shift(temp_ehf,3+2) + shift(temp_ehf,4+2) + shift(temp_ehf,5+2)
#                                 + shift(temp_ehf,6+2) + shift(temp_ehf,7+2) + shift(temp_ehf,8+2) + shift(temp_ehf,9+2) + shift(temp_ehf,10+2)
#                                 + shift(temp_ehf,11+2) + shift(temp_ehf,12+2) + shift(temp_ehf,13+2) + shift(temp_ehf,14+2) + shift(temp_ehf,15+2)
#                                 + shift(temp_ehf,16+2) + shift(temp_ehf,17+2) + shift(temp_ehf,18+2) + shift(temp_ehf,19+2) + shift(temp_ehf,20+2)
#                                 + shift(temp_ehf,21+2) + shift(temp_ehf,22+2) + shift(temp_ehf,23+2) + shift(temp_ehf,24+2) + shift(temp_ehf,25+2)
#                                 + shift(temp_ehf,26+2) + shift(temp_ehf,27+2) + shift(temp_ehf,28+2) + shift(temp_ehf,29+2) + shift(temp_ehf,30+2))/30] # retrospective (with +2)
#   brambilla.all[City==i, dmt3f := rollmean(temp_ehf, 3, fill=NA, align='left')]
#   brambilla.all[City==i, dmtp30f := (shift(temp_ehf,1) + shift(temp_ehf,2) + shift(temp_ehf,3) + shift(temp_ehf,4) + shift(temp_ehf,5)
#                                  + shift(temp_ehf,6) + shift(temp_ehf,7) + shift(temp_ehf,8) + shift(temp_ehf,9) + shift(temp_ehf,10)
#                                  + shift(temp_ehf,11) + shift(temp_ehf,12) + shift(temp_ehf,13) + shift(temp_ehf,14) + shift(temp_ehf,15)
#                                  + shift(temp_ehf,16) + shift(temp_ehf,17) + shift(temp_ehf,18) + shift(temp_ehf,19) + shift(temp_ehf,20)
#                                  + shift(temp_ehf,21) + shift(temp_ehf,22) + shift(temp_ehf,23) + shift(temp_ehf,24) + shift(temp_ehf,25)
#                                  + shift(temp_ehf,26) + shift(temp_ehf,27) + shift(temp_ehf,28) + shift(temp_ehf,29) + shift(temp_ehf,30))/30] # retrospective (with +2)
#   brambilla.all[City==i, dmhi3 := rollmean(hi_ehf, 3, fill=NA, align='right')]
#   brambilla.all[City==i, dmhip30 := (shift(hi_ehf,1+2) + shift(hi_ehf,2+2) + shift(hi_ehf,3+2) + shift(hi_ehf,4+2) + shift(hi_ehf,5+2)
#                                  + shift(hi_ehf,6+2) + shift(hi_ehf,7+2) + shift(hi_ehf,8+2) + shift(hi_ehf,9+2) + shift(hi_ehf,10+2)
#                                  + shift(hi_ehf,11+2) + shift(hi_ehf,12+2) + shift(hi_ehf,13+2) + shift(hi_ehf,14+2) + shift(hi_ehf,15+2)
#                                  + shift(hi_ehf,16+2) + shift(hi_ehf,17+2) + shift(hi_ehf,18+2) + shift(hi_ehf,19+2) + shift(hi_ehf,20+2)
#                                  + shift(hi_ehf,21+2) + shift(hi_ehf,22+2) + shift(hi_ehf,23+2) + shift(hi_ehf,24+2) + shift(hi_ehf,25+2)
#                                  + shift(hi_ehf,26+2) + shift(hi_ehf,27+2) + shift(hi_ehf,28+2) + shift(hi_ehf,29+2) + shift(hi_ehf,30+2))/30] # retrospective (with +2)
#   brambilla.all[City==i, dmhi3f := rollmean(hi_ehf, 3, fill=NA, align='left')]
#   brambilla.all[City==i, dmhip30f := (shift(hi_ehf,1) + shift(hi_ehf,2) + shift(hi_ehf,3) + shift(hi_ehf,4) + shift(hi_ehf,5)
#                                   + shift(hi_ehf,6) + shift(hi_ehf,7) + shift(hi_ehf,8) + shift(hi_ehf,9) + shift(hi_ehf,10)
#                                   + shift(hi_ehf,11) + shift(hi_ehf,12) + shift(hi_ehf,13) + shift(hi_ehf,14) + shift(hi_ehf,15)
#                                   + shift(hi_ehf,16) + shift(hi_ehf,17) + shift(hi_ehf,18) + shift(hi_ehf,19) + shift(hi_ehf,20)
#                                   + shift(hi_ehf,21) + shift(hi_ehf,22) + shift(hi_ehf,23) + shift(hi_ehf,24) + shift(hi_ehf,25)
#                                   + shift(hi_ehf,26) + shift(hi_ehf,27) + shift(hi_ehf,28) + shift(hi_ehf,29) + shift(hi_ehf,30))/30] # retrospective (with +2)
# }
# brambilla.all[, dmt95 := quantile(temp_ehf, 0.95, na.rm=T), by=City] # 95th quantile of dmt (average air t), by city
# brambilla.all[, dmhi95 := quantile(hi_ehf, 0.95, na.rm=T), by=City] # 95th quantile of dmhi (average hi), by city
# brambilla.all[, ehisig := dmt3 - dmt95] # EHIsig
# brambilla.all[, ehisigf := dmt3f - dmt95] # EHIsig forward
# brambilla.all[, ehiaccl := dmt3 - dmtp30] # EHIaccl
# brambilla.all[, ehiacclf := dmt3f - dmtp30f] # EHIaccl forward
# brambilla.all[, ehisig.hi := dmhi3 - dmhi95] # EHIsig HI
# brambilla.all[, ehisigf.hi := dmhi3f - dmhi95] # EHIsig forward HI
# brambilla.all[, ehiaccl.hi := dmhi3 - dmhip30] # EHIaccl HI
# brambilla.all[, ehiacclf.hi := dmhi3f - dmhip30f] # EHIaccl forward HI
# brambilla.all[, ehf := ehisig * pmax(1, ehiaccl)] # EHF
# brambilla.all[, ehff := ehisigf * pmax(1, ehiacclf)] # EHF forward
# brambilla.all[, ehf.hi := ehisig.hi * pmax(1, ehiaccl.hi)] # EHF HI
# brambilla.all[, ehff.hi := ehisigf.hi * pmax(1, ehiacclf.hi)] # EHF forward EHF
#
# ## EHF DMT95 reference values
# brambilla_ehfr <- unique(brambilla.all[, .(City,dmt95,dmhi95)])
#
# ### Save
# save(brambilla.all, brambilla_ehfr, file='brambilla.all.rda')

load(file='brambilla.all.rda') # Climate data



##################################################
### Create simulated claims data
##################################################

# Stratum variables
ds.city <- sort(unique(c('Melbourne','Sydney'))) # 2 cities (Melbourne and Sydney)
length.ds.city <- length(ds.city) # Number of cities

## Create randomly-generated claims data with Poisson and Tweedie assumed distributions for number of OIIs (claims) and total costs, respectively
## This simulated data has no relationship to any predictor variables
dates <- rep(seq(min(brambilla.all[,Date]), max(brambilla.all[,Date]), 1),each=length.ds.city)
ldates <- length(dates)/length.ds.city
set.seed(3) # for random distributions
claims <- data.table(
  Date = dates,
  City = rep(c('Melbourne','Sydney'),ldates),
  'Number of OIIs' = c(rpois(ldates,100), rpois(ldates,180)), # Number of OIIs
  'Total costs' = c(rtweedie(ldates, xi=1.5, mu=300, phi=8), rtweedie(ldates, xi=1.4, mu=600, phi=11)) # Total costs
)
claims[,'Costs per OII':=get('Total costs')/get('Number of OIIs')] # Costs per OII, should you wish to model this

# By variables for commands
by.vars <- c("Date","City")
by.vars1 <- c('Year', replace(by.vars, by.vars=="Date", 'Month'))


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

load(file='school.holidays.rda')
load(file='pop.rda') # Worker's population using ABS data. Due to estimating indoor/outdoor proportions, n can have decimal places


daily <- merge(merge(merge(claims, brambilla.all, by=c('Date','City')),
               public.holidays[!is.na(phol)], by=c('Date','City'), all.x=T),
                school.holidays, by=c('Date','City'), all.x=T)
daily[is.na(phol), phol := 0] # Replace NA in public holidays with 0

# Merge worker's monthly population
daily[, ':='(Year=year(Date), Month=month(Date), Day=day(Date), dow=weekdays(Date))]
pop <- pop[, lapply(.SD, sum, na.rm=T), by=mget(by.vars1), .SDcols = c("n","fulltime","parttime")] # combine Indoors and Outdoos
daily.ds <- merge(daily, pop, by=c("Year","Month","City"), all.x=T)



##################################################
### Rename climate variables
##################################################

exposure.old <- c('temp.max','temp.min',
                  'rh.max','rh.min','rh.ave',
                  'hi.max','hi.min',
                  'ehf','ehff',
                  'ehf.hi','ehff.hi')
exposure.new <- c("Maximum temperature","Minimum temperature",
                  'Maximum relative humidity','Minimum relative humidity','Average relative humidity',
                  'Maximum heat index','Minimum heat index',
                  'Excess heat factor','Excess heat factor (forward)',
                  'Excess heat index factor','Excess heat index factor (forward)')
setnames(daily.ds, old=exposure.old, new=exposure.new)  # rename exposures



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
