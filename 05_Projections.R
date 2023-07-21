
##################################################
### Prepare future meteorological data and save as .rda file
###   This step is already done for you and is commented out
###   However, the code is included to demonstrate how the original data was extracted and edited
###   This dataset and code, along with the original datasets used to derive it, are available on FigShare
##################################################
# 
# # Libraries
# library(raster); library(ncdf4)
# library(data.table)
# library(dplyr)
# library(HeatStress) # apparent temperature calculations by Casanueva. Need to run the following commands to install: # install.packages("devtools") # devtools::install_github("anacv/HeatStress")
# library(lubridate) 
# library(readr) # import csv
# library(readxl) # import Excel files 
# library(stringr) # string commands
# library(zoo)

# # Loop components
# cities <- c('Adelaide','Brisbane','Canberra','Darwin','Hobart','Melbourne','Perth','Sydney')
# vrbl <- c('hurs','rsds','tas','tasmax','tasmin')
# gcm <- c('ACCESS1-0','CanESM2','CESM1-CAM5','CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','MIROC5','NorESM1-M')
# rcp <-  c('rcp45','rcp85')
# prd <- c('2016-2045','2036-2065')
# 
# byvar <- c('City','Model','RCP','Period','Date')
# vara <- c('City','Model','RCP','Period','Date','value')
# 
# ## Centroids
# coords <- dput(structure(list(city = c("Adelaide", "Brisbane", "Darwin", "Hobart", "Melbourne", "Perth", "Sydney", "Canberra"),
#                lon = c(138.6244, 153.02, 130.84, 147.3, 144.95, 115.9, 151.2, 148.99),
#                lat = c(-34.92661, -27.47, -12.46, -42.85, -37.8, -32, -33.85, -27.51010751)),
#           row.names = c(NA, -8L),
#           class = c("data.table", "data.frame"), 
#           index = structure(integer(0),"`__city`" = c(1L, 2L, 8L, 3L, 4L, 5L, 6L, 7L))))
# 
# ## Extract files and save as csv.s
# 
# # cciaf.loc <- paste0("your_directory_with_the_files") # Location of Climate Change in Australia (CCIA) files
# # 
# time <- proc.time()[3]
# for(i in cities) {
#   for(a in vrbl) {
#     for(b in gcm) {
#       for(c in rcp) {
#         for(d in prd) {
#           .all <- paste(i,a,b,c,d, sep='_') # unique name combination that is easily referenced in R
#           .string <- paste0(cciaf.loc, i, '/', a, '_aus_', b, '_', c, '_r1i1p1_CSIRO-MnCh-wrt-1986-2005-Scl_v1_day_', d, '_g-latSSlonEE.nc') # filename
#           skip_to_next <- F # reset object result from tryCatch
#           tryCatch(.brick <- brick(.string), error = function(e) { skip_to_next <<- TRUE}) # attempt to extract file as brick. Error will occur if file doesn't exist, save as T
#           if(skip_to_next) { next } # if error, skip to next loop
#           .obs <- raster::extract(.brick, coords[city==i,.(lon,lat)], buffer=10^20, fun=mean) # buffer of 100000 correctly gets mean of 7*7 grid. Increasing to be safe as only coordinates are 7*7 grid. na.rm=T not required, though assumingly would be if NAs
#           .obs1 <- t(.obs) # convert to a 1 column matrix of the values
#           colnames(.obs1) <- str_remove_all(.all,'-') # unique column
#           write.csv(.obs1, paste0(cciaf.loc,'_csv/',i,'_',a,'_',b,'_',c,'_',d,'.csv'))
#         }
#       }
#     }
#   }
# } 
# proc.time()[3]-time
# 
# 
# ## Load files and combine
# time <- proc.time()[3]
# ccia_data <- list()
# for(i in cities) { # combine into a list, and assign loop variables into data
#   for(b in gcm) {
#     for(c in rcp) {
#       for(d in prd) {
#         for(a in vrbl) {
#           .all <- paste(i,a,b,c,d, sep='_') # unique name combination that is easily referenced in R
#           tryCatch(ccia_data[[.all]] <- cbind(i,b,c,d,a, read_csv(paste0(cciaf.loc,'_csv/',i,'_',a,'_',b,'_',c,'_',d,'.csv'))),
#                    error = function(e) { skip_to_next <<- TRUE}) # attempt to save to list. Error will occur if file doesn't exist, save as T
#           if(skip_to_next) { next } # if error, skip to next loop
#         }
#       }
#     }
#   }
# }
# for(z in seq(length(ccia_data))) { # combine into a list, and assign loop variables into data
#   names(ccia_data[[z]]) <- c('City','Model','RCP','Period','variable','Date','value')
# }
# proc.time()[3]-time
# 
# time <- proc.time()[3]
# ccia_future <- as.data.table(do.call(rbind.data.frame, ccia_data)) # combine lists into one data frame
# proc.time()[3]-time
# detach("package:raster", unload = TRUE) # remove raster, which loads terra and its shift that interferes with data.table shift
# 
# 
# ## Make wide by climate variable
# ccia_future[,Date:=as_date(str_remove(Date,'X'))] # convert Date to an actual Date
# time <- proc.time()[3]
# ccia_future1 <- merge(merge(merge(merge(rename(ccia_future[variable=='tas',mget(vara)], c('tas'='value')),
#                                         rename(ccia_future[variable=='tasmax',mget(vara)], c('tasmax'='value')), by=byvar),
#                                   rename(ccia_future[variable=='tasmin',mget(vara)], c('tasmin'='value')), by=byvar),
#                             rename(ccia_future[variable=='hurs',mget(vara)], c('hurs'='value')), by=byvar),
#                       rename(ccia_future[variable=='rsds',mget(vara)], c('rsds'='value')), by=byvar)
# proc.time()[3]-time
# 
# 
# ## Create apparent T metrics
# source(paste(code.loc,'heat.index2.r',sep='')) # slightly modified Heat Index for proper rounding
# ccia_future1[hurs>100, hurs:=100] # clouds should not be forming at 2m above surface, so cap RH at 100%
# ccia_future1[,ef := ef.fn(tas)] # Enhancement factor based on Buck's 1996 and air P of 1013.25 (default). EF is very small regardless of air P
# ccia_future1[,svp := vp.fn(tas, ef.=ef)] # SVP based on Buck's 1996 and air P of 1013.25 (default). May be less accurate without air P variable
# ccia_future1[,vp:=hurs*svp/100] # VP, by rearranging calc used to calculate RH with BARRA data. May be less accurate without air P variable
# ccia_future1[,dewpt:=vp.to.tdp.fn(vp.=vp, t.=tas)] # dew point T. Reverses effect of EF, so not affected by lack of air P variable
# ccia_future1[,sh:=sh.vp.fn(vp, ef.=ef)] # specific humidity. Affected by lack of air P and EF
# ccia_future1[,heat.index := heat.index2(t=tas, rh=hurs, temperature.metric='celsius', round=9999999),] # Heat index. Use self-code that allows rounding to >2 digits and matches website exactly
# ccia_future1[,heat.index.max := heat.index2(t=tasmax, rh=hurs, temperature.metric='celsius', round=9999999),]
# ccia_future1[,humidex := tas + 5/9 * (vp - 10)] # Humidex
# ccia_future1[,humidex.max := tasmax + 5/9 * (vp - 10)] # Humidex
# ccia_future1[,wbt := wbt.Stull(tas, hurs)] # Wet bulb T using Stull's approximation
# ccia_future1[,di := (tas+hurs)/2] # Discomfort index
# ccia_future1$wbgtb <- with(ccia_future1, wbgt.Bernard(tas=tas, dewp=dewpt))[["data"]] # shade WBGT. must reference "data", else list also with Tpwb (psychrometric wet bulb T)
# ccia_future1$wbgtb.max <- with(ccia_future1, wbgt.Bernard(tas=tasmax, dewp=dewpt))[["data"]] # shade WBGT. must reference "data", else list also with Tpwb (psychrometric wet bulb T)
# ccia_future1[,steadman.id := 0.89*tas + 0.382*vp - 2.56] # indoor, no ws or radiation
# ccia_future1[,steadman.id.max := 0.89*tasmax + 0.382*vp - 2.56] # indoor, no ws or radiation
# 
# 
# ## Calculate average T using EHF (lead min)
# ccia_future1[City==shift(City,-1) & Model==shift(Model,-1) & RCP==shift(RCP,-1) & Period==shift(Period,-1), tasmin_lead:=shift(tasmin,-1)] # EHF uses min T from same day
# ccia_future1[City=='Sydney' & Model=='NorESM1-M' & RCP=='rcp85' & Period=='2036-2065' & Date>='2065-12-30', tasmin_lead:=shift(tasmin,-1)] # above misses 2nd to last day, which this addresses
# ccia_future1[, tas_ehf := rowMeans(.SD), .SDcols = c("tasmax","tasmin_lead")]
# 
# ## DMT3 and DMT over last 30days. As only have average relative humidity, EHF_HI calculated using same max/min instead of adjusting min as per T
# time <- proc.time()[3]
# for(a in cities) {
#   for(b in gcm) {
#     for(c in rcp) {
#       for(d in prd) {
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmt3 := rollmean(tas_ehf, 3, fill=NA, align='right')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmtp30 := (shift(tas_ehf,1+2) + shift(tas_ehf,2+2) + shift(tas_ehf,3+2) + shift(tas_ehf,4+2) + shift(tas_ehf,5+2)
#                                                                          + shift(tas_ehf,6+2) + shift(tas_ehf,7+2) + shift(tas_ehf,8+2) + shift(tas_ehf,9+2) + shift(tas_ehf,10+2)
#                                                                          + shift(tas_ehf,11+2) + shift(tas_ehf,12+2) + shift(tas_ehf,13+2) + shift(tas_ehf,14+2) + shift(tas_ehf,15+2)
#                                                                          + shift(tas_ehf,16+2) + shift(tas_ehf,17+2) + shift(tas_ehf,18+2) + shift(tas_ehf,19+2) + shift(tas_ehf,20+2)
#                                                                          + shift(tas_ehf,21+2) + shift(tas_ehf,22+2) + shift(tas_ehf,23+2) + shift(tas_ehf,24+2) + shift(tas_ehf,25+2)
#                                                                          + shift(tas_ehf,26+2) + shift(tas_ehf,27+2) + shift(tas_ehf,28+2) + shift(tas_ehf,29+2) + shift(tas_ehf,30+2))/30] # retrospective (with +2)
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmt3f := rollmean(tas_ehf, 3, fill=NA, align='left')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmtp30f := (shift(tas_ehf,1) + shift(tas_ehf,2) + shift(tas_ehf,3) + shift(tas_ehf,4) + shift(tas_ehf,5)
#                                                                           + shift(tas_ehf,6) + shift(tas_ehf,7) + shift(tas_ehf,8) + shift(tas_ehf,9) + shift(tas_ehf,10)
#                                                                           + shift(tas_ehf,11) + shift(tas_ehf,12) + shift(tas_ehf,13) + shift(tas_ehf,14) + shift(tas_ehf,15)
#                                                                           + shift(tas_ehf,16) + shift(tas_ehf,17) + shift(tas_ehf,18) + shift(tas_ehf,19) + shift(tas_ehf,20)
#                                                                           + shift(tas_ehf,21) + shift(tas_ehf,22) + shift(tas_ehf,23) + shift(tas_ehf,24) + shift(tas_ehf,25)
#                                                                           + shift(tas_ehf,26) + shift(tas_ehf,27) + shift(tas_ehf,28) + shift(tas_ehf,29) + shift(tas_ehf,30))/30] # retrospective (with +2)
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmhi3 := rollmean(heat.index, 3, fill=NA, align='right')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmhip30 := (shift(heat.index,1+2) + shift(heat.index,2+2) + shift(heat.index,3+2) + shift(heat.index,4+2) + shift(heat.index,5+2)
#                                                                           + shift(heat.index,6+2) + shift(heat.index,7+2) + shift(heat.index,8+2) + shift(heat.index,9+2) + shift(heat.index,10+2)
#                                                                           + shift(heat.index,11+2) + shift(heat.index,12+2) + shift(heat.index,13+2) + shift(heat.index,14+2) + shift(heat.index,15+2)
#                                                                           + shift(heat.index,16+2) + shift(heat.index,17+2) + shift(heat.index,18+2) + shift(heat.index,19+2) + shift(heat.index,20+2)
#                                                                           + shift(heat.index,21+2) + shift(heat.index,22+2) + shift(heat.index,23+2) + shift(heat.index,24+2) + shift(heat.index,25+2)
#                                                                           + shift(heat.index,26+2) + shift(heat.index,27+2) + shift(heat.index,28+2) + shift(heat.index,29+2) + shift(heat.index,30+2))/30] # retrospective (with +2)
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmhi3f := rollmean(heat.index, 3, fill=NA, align='left')] # above misses 2 Adelaide dates that it shouldn't, which this addresses
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmhip30f := (shift(heat.index,1) + shift(heat.index,2) + shift(heat.index,3) + shift(heat.index,4) + shift(heat.index,5)
#                                                                            + shift(heat.index,6) + shift(heat.index,7) + shift(heat.index,8) + shift(heat.index,9) + shift(heat.index,10)
#                                                                            + shift(heat.index,11) + shift(heat.index,12) + shift(heat.index,13) + shift(heat.index,14) + shift(heat.index,15)
#                                                                            + shift(heat.index,16) + shift(heat.index,17) + shift(heat.index,18) + shift(heat.index,19) + shift(heat.index,20)
#                                                                            + shift(heat.index,21) + shift(heat.index,22) + shift(heat.index,23) + shift(heat.index,24) + shift(heat.index,25)
#                                                                            + shift(heat.index,26) + shift(heat.index,27) + shift(heat.index,28) + shift(heat.index,29) + shift(heat.index,30))/30] # retrospective (with +2)
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmt95a := quantile(tas_ehf, 0.95, na.rm=T)] # 95th quantile of dmt (average air t). Uses CCiA data to represent adaptation for given RCP and period
#         ccia_future1[City==a & Model==b & RCP==c & Period==d, dmhi95a := quantile(heat.index, 0.95, na.rm=T)] # 95th quantile of dmhi (average hi). Uses CCiA data to represent adaptation
#       }
#     }
#   }
# }
# proc.time()[3]-time
# 
# ## CCiA 2030 EHF thresholds: for testing partial adaptation in 2050. Partial referring to adaptation to a recent but not the current period
# ccia_future1[, lapply(.SD, mean, na.rm=T), by=c('City','Period','RCP','Model'), .SDcols=c('dmt95a','dmhi95a')] # day of month
# ccia_30ehfa <- ccia_future1[Period=='2016-2045', lapply(.SD, mean, na.rm=T), by=c('City','RCP','Model'), .SDcols=c('dmt95a','dmhi95a')] # day of month
# setnames(ccia_30ehfa, old=c('dmt95a','dmhi95a'), new=c('dmt95a30','dmhi95a30'), skip_absent=T)
# 
# ## Add BARRA EHF and CCiA 2030 thresholds
# load(file=paste0(barra.loc,'barra_ehfr.rda')) # DMT95 thresholds from BARRA data, barra_ehfr
# ccia_future2 <- merge(merge(ccia_future1,barra_ehfr, by='City'), ccia_30ehfa, by=c('City','RCP','Model'), all.x=T)
# 
# ## EHIsig, btoh with and without adaptation
# for(a in cities) {
#   ccia_future2[City==a, ehisig := dmt3 - dmt95] # EHIsig
#   ccia_future2[City==a, ehisigf := dmt3f - dmt95] # EHIsig forward
#   ccia_future2[City==a, ehisig.hi := dmhi3 - dmhi95] # EHIsig HI
#   ccia_future2[City==a, ehisigf.hi := dmhi3f - dmhi95] # EHIsig forward HI
#   ccia_future2[City==a, ehisiga := dmt3 - dmt95a] # EHIsig with adaptation
#   ccia_future2[City==a, ehisigfa := dmt3f - dmt95a] # EHIsig forward with adaptation
#   ccia_future2[City==a, ehisiga.hi := dmhi3 - dmhi95a] # EHIsig HI with adaptation
#   ccia_future2[City==a, ehisigfa.hi := dmhi3f - dmhi95a] # EHIsig forward HI with adaptation
#   ccia_future2[City==a, ehisigpa := dmt3 - dmt95a30] # EHIsig with partial adaptation (identical to full adaptation in 2030)
#   ccia_future2[City==a, ehisigfpa := dmt3f - dmt95a30] # EHIsig forward with partial adaptation
#   ccia_future2[City==a, ehisigpa.hi := dmhi3 - dmhi95a30] # EHIsig HI with partial adaptation
#   ccia_future2[City==a, ehisigfpa.hi := dmhi3f - dmhi95a30] # EHIsig forward HI with partial adaptation
# }
# 
# ## EHIaccl
# ccia_future2[, ehiaccl := dmt3 - dmtp30] # EHIaccl
# ccia_future2[, ehiacclf := dmt3f - dmtp30f] # EHIaccl forward
# ccia_future2[, ehiaccl.hi := dmhi3 - dmhip30] # EHIaccl HI
# ccia_future2[, ehiacclf.hi := dmhi3f - dmhip30f] # EHIaccl forward HI
# 
# ## EHF
# ccia_future2[, ehf := ehisig * pmax(1, ehiaccl)] # EHF
# ccia_future2[, ehff := ehisigf * pmax(1, ehiacclf)] # EHF forward
# ccia_future2[, ehf.hi := ehisig.hi * pmax(1, ehiaccl.hi)] # EHF HI
# ccia_future2[, ehff.hi := ehisigf.hi * pmax(1, ehiacclf.hi)] # EHF forward EHF
# 
# ccia_future2[, ehfa := ehisiga * pmax(1, ehiaccl)] # EHF with adaptation
# ccia_future2[, ehffa := ehisigfa * pmax(1, ehiacclf)] # EHF forward with adaptation
# ccia_future2[, ehfa.hi := ehisiga.hi * pmax(1, ehiaccl.hi)] # EHF HI with adaptation
# ccia_future2[, ehffa.hi := ehisigfa.hi * pmax(1, ehiacclf.hi)] # EHF forward EHF with adaptation
# 
# ccia_future2[, ehfpa := ehisigpa * pmax(1, ehiaccl)] # EHF with partial adaptation
# ccia_future2[, ehffpa := ehisigfpa * pmax(1, ehiacclf)] # EHF forward with adaptation
# ccia_future2[, ehfpa.hi := ehisigpa.hi * pmax(1, ehiaccl.hi)] # EHF HI with adaptation
# ccia_future2[, ehffpa.hi := ehisigfpa.hi * pmax(1, ehiacclf.hi)] # EHF forward EHF with adaptation
# 
# ## Define projected outcome series as historical average per day of year, then repeated along same projection period of modelled temperatures series
# ccia_future2[,Month:=month(Date)] # define Month to limit to warm season. much faster to do outside of loop
# ccia_future2[RCP %in% c('rcp26','RCP26'),RCP:='RCP2.6'] # format RCP
# ccia_future2[RCP %in% c('rcp45','RCP45'),RCP:='RCP4.5']
# ccia_future2[RCP %in% c('rcp60','RCP60'),RCP:='RCP6.0']
# ccia_future2[RCP %in% c('rcp85','RCP85'),RCP:='RCP8.5']
# 
# 
# ### Save
# ccia_future2 <- ccia_future2[City %in% c('Melbourne','Sydney')] # for this reproducible example, only Melbourne and Sydney
# ccia_future2 <- ccia_future2[Model %in% c('ACCESS1-0','CanESM2')] # for this reproducible example, only two GCMs
# ccia_future2[,':='(humidex=NULL,humidex.max=NULL,wbt=NULL,di=NULL,wbgtb=NULL,wbgtb.max=NULL,steadman.id=NULL,steadman.id.max=NULL)] # if these variables exist, remove them
# save(ccia_future2, file='ccia_future2.rda')
load(file='ccia_future2.rda')


##################################################
#### Projected population data: https://www.abs.gov.au/statistics/people/population/population-projections-australia/latest-release
##   This step is already done for you and is commented out
##   However, the code is included to demonstrate how the original data was extracted and edited
###################################################

# # abs.loc <- paste0("your_directory_with_the_files") # Location of Australian Bureau of Statistics (ABS) projected population dataset
# # abs.loc <- paste0("/Users/MatthewBorg/Library/CloudStorage/Box-Box/Data original/ABS/") # Location of ABS projected population dataset
# # abs.loc <- paste0("C:/Users/a1210385/Box/Data original/CCIA/future/") # Location of ABS projected population dataset
# 
# projpop.fn <- function(city., cityno.) {
#   pop.proj.by <- c('Category','Series')
#   pop.proj <- as.data.table(read_excel(paste0(abs.loc,'Projected population/32220ds01_2017-2066_summary_statistics.xls'), sheet=cityno., skip=5, col_names=T)) # import
#   colnames(pop.proj) <- c(pop.proj.by,'2017',colnames(pop.proj)[-(1:3)]) # change 1st 3 column names
#   pop.proj1 <- pop.proj[str_detect(Category, 'End population'), lapply(.SD, as.numeric), by=pop.proj.by] # keep total population, make non-by columns numeric
# 
#   pop2030 <- rowMeans(pop.proj1[,mget(paste0(c(2017:2044)))]) # population for 2030. (30 June) 2017 is earliest date in dataset, so exclude last year to compensate
#   pop2050 <- rowMeans(pop.proj1[,mget(paste0(c(2036:2065)))]) # population for 2050
#   pop <- pop.proj1[,get(paste0(2017))] # population for 2050
#   .pop <- data.table(City=rep(city.,3), Series=c('A','B','C'), 'n'=pop, 'n2030'=pop2030, 'n2050'=pop2050)
#   return(.pop)
# }
# 
# ## Calculate population per city and combine
# projpop.sy <- projpop.fn('Sydney',4)
# projpop.me <- projpop.fn('Melbourne',7)
# projpop.br <- projpop.fn('Brisbane',10)
# projpop.ad <- projpop.fn('Adelaide',13)
# projpop.pe <- projpop.fn('Perth',16)
# projpop.ho <- projpop.fn('Hobart',19)
# projpop.da <- projpop.fn('Darwin',22)
# projpop.ca <- projpop.fn('Canberra',24)
# projpop.c <- rbindlist(list(projpop.sy, projpop.me, projpop.br, projpop.ad, projpop.pe, projpop.ho, projpop.da, projpop.ca))
# 
# ## Add total across cities
# projpop.t <- projpop.c[, lapply(.SD, sum, na.rm=T), keyby=.(Series), .SDcols = c('n','n2030','n2050')]
# projpop.t[,City:='Total']
# setcolorder(projpop.t, c('City',colnames(projpop.t)[-length(projpop.t)])) # Move City to be 1st column
# projpop <- rbindlist(list(projpop.c, projpop.t))
# 
# ## Calculate population ratios for projected and base
# projpop[,':='(r2030=n2030/n, r2050=n2050/n)]
# 
# ## Make dataset vertical
# projpop1 <- melt(projpop, id.vars=c('City','Series'), measure.vars=c('r2030','r2050'))
# colnames(projpop1) <- c('City','Series','Period','popratio')
# projpop1[str_detect(Period,'2030'), Period:='2016-2045'] # match projected value
# projpop1[str_detect(Period,'2050'), Period:='2036-2065']
# 
# 
# ## Make Series horizontal
# projpop2 <- merge(merge(projpop1[Series=='A',-"Series"], projpop1[Series=='B',-"Series"], by=c('City','Period')), projpop1[Series=='C',-"Series"], by=c('City','Period'))
# colnames(projpop2) <- c('City','Period','popratioA','popratioB','popratioC')
# 
# 
# ## Save dataset
# save(projpop2, file=paste0('projpop.rda'))
load(paste0('projpop.rda'))



########################################################################################################
# Prepare projected data. No leap days, like in Vicedo-Cabrera example
########################################################################################################

# Stage 3 folders
projresults <- paste0(base,'/Projections/')
dir.create(paste0(base,'/Projections')) # Warning if exists (doesn't replace)
adresults <- paste0(projresults,'Adapt',adapt,'/')
dir.create(paste0(projresults,'Adapt',adapt)) # Warning if exists (doesn't replace)

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
write.csv(dmt, file=paste0(projresults,'Projected DMT.csv'), na='', row.names=F) # create csv file


## Rename exposures to match daily.ds
ccia.newnames <- c('ave_sh','Excess heat factor','Excess heat factor (forward)','Excess heat index factor','Excess heat factor index (forward)')
if(adapt==1) { # use EHF metrics with adaptation
  setnames(ccia_future2, old=c('sh','ehfa','ehffa','ehfa.hi','ehffa.hi'), new=ccia.newnames, skip_absent=T) 
} else if(adapt==0.5) { # use EHF metrics with partial adaptation
  setnames(ccia_future2, old=c('sh','ehfpa','ehffpa','ehfpa.hi','ehffpa.hi'), new=ccia.newnames, skip_absent=T) 
} else { # use EHF metrics without adaptation
  setnames(ccia_future2, old=c('sh','ehf','ehff','ehf.hi','ehfa.hi'), new=ccia.newnames, skip_absent=T) 
}



########################################################################################################
# Loop components / dimensions
########################################################################################################

# gcm <- c('ACCESS1-0','CanESM2','CESM1-CAM5','CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','MIROC5','NorESM1-M') # General circulation models used in this study
gcm <- c('ACCESS1-0','CanESM2') # General circulation models used in this reprdoucible example (to reduce dataset size)
rcp <-  c('RCP4.5','RCP8.5') # not all have rcp26 or rcp60
prd <- c('2016-2045','2036-2065') 
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
  oiproj30[[a]] <- rep(oidoy[[a]], length=length(unique(ccia_future2[City=='Sydney' & Model==gcm[1] & RCP==rcp[1] & Period==prd[1] & Month %in% c(10,11,12,1,2,3),Date]))) # repeat along 2016-45. Choice of city, gcm, rcp and prd doesn't matter
  oiproj50[[a]] <- rep(oidoy[[a]], length=length(unique(ccia_future2[City=='Sydney' & Model==gcm[1] & RCP==rcp[1] & Period==prd[2] & Month %in% c(10,11,12,1,2,3),Date]))) # repeat along 2036-65
}

oiperiod[tot.city[length(tot.city)]] <- sum(oiperiod[-length(tot.city)]) # add total to oiperiod



########################################################################################################
# (1) Extrpaolate exposure-response curve, (2) project and quantify ANNUAL impact and (3) ensemble estimates & quantification of uncertainty
########################################################################################################

# Define objects
ansim <- past_ansim <- array(NA, dim=c(length(prd), length(sim.names), length(absrel), length(gcm), length(rcp), length(tot.city), nsim2+1), 
                             dimnames=list(prd, sim.names, absrel, gcm, rcp, tot.city, c("est",paste0("sim",seq(nsim2))))) # Array to store AN + a dummy for past results

proj.exp.size <- nrow(ccia_future2[City=='Sydney' & RCP==rcp[1] & Model==gcm[1] & Period==prd[1] & Month %in% c(10,11,12,1,2,3)]) # Choice of city, gcm, rcp and prd doesn't matte
projexpmark <- matrix(NA, proj.exp.size, length(sim.names))
colnames(projexpmark) <- sim.names


gc()
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
# loop takes about 5-10 sec on Apple M1 Max (Darwin Kernel Version 21.6.0) per GCM


## Add national results to array, filling in the Total/National category
ansim[,,,,,tot.city[length(tot.city)],] <- apply(ansim, c(1:5,7), sum, na.rm=T) # sum by all dimensions but city (sums city results)


## Estimate AN relative "rel" to baseline (absolute/projected - baseline) for all cities and total, assuming no change in pop
for(a in tot.city) { # fill in dummy ansim with past results, including total/national
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
}



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
write.csv(projaf2, file=paste0(adresults, 'Fractions.csv'), na='', row.names=F) # create csv file



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
write.csv(projan2, file=paste0(adresults, 'Numbers.csv'), na='', row.names=F) # create csv file



######################### END ############################
######################### END ############################
  