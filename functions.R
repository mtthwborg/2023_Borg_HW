

## Calculate qAIC. fqaic.fn is adapted on fqaic from Gasparrini et al. 2015 Lancet
## model. Model outputted from glm() or mgcv::gam(). If glm() is not a quasimodel, qAIC will equal AIC (where dispersion factor=1)
## gam. Specify TRUE if gam() used, FALSE otherwise (default)
fqaic.fn <- function(model, gam=F) {
  .loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE)) # log-likelihood
  .phi <- summary(model)$dispersion # dispersion factor
  if(gam==T) { # if used gam() from mgcv
    .np <- summary(model)$dispersion
  } else { # if used glm from stats package
    .np <- summary(model)$df[3]
  }
  qaic <- -2*.loglik + 2*.np*.phi # QAIC
  return(qaic)
}    


## Perform Wald test with mixmeta or mvmeta object
## fwald2 from Gasparrini et al. 2015 Lancet with additional "short" option (if T, also outputs Wald statistic and df)
fwald2.fn <- function(model, var, short=F) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  p <- 1-pchisq(waldstat,df)
  if(isTRUE(short)) { # return p-value only
    return(p)
  } else { # return statistic, df and p-value
    results <- c(waldstat,df,p)
    names(results) <- c('Wald statistic','df','p-value')
    return(results)
  }
}

## Report mixmeta Cochran Q test, I^2 and AIC statistics from a mixmeta or mvmeta object (mv.)
mv.results.fn <- function(mv.) {
  sum.mv.q <- summary(mv.)[["qstat"]] # Q-stat values. Can also get with qtest()
  .vars <- attr(mv.[['terms']],"term.labels") # fixed variables
  if(identical(.vars,character(0))) { # if no fixed predictors
    .vars.wald <- NULL
  } else {
    .vars.wald <- c(sapply(.vars, function(w) fwald2.fn(mv., w))) # results from multivariate Wald test
    names(.vars.wald) <- paste(rep(.vars, each=3), c('W','df','p'), sep = ".")
  }
  .mv.results <- c('Q-statistic'=sum.mv.q$Q[1], 'df'=sum.mv.q$df[1], 'P-value'=sum.mv.q$pvalue[1], 'I^2'=summary(mv.)[["i2stat"]][1], 'AIC'=AIC(mv.), .vars.wald)
  names(.mv.results) <- str_remove_all(names(.mv.results), '..all') # remove ..all  from mixmeta code
  return(.mv.results) # return combined results
}


## Wrap format() around round(). Enables more control over presentation of rounded results. Converts results to character
# round.fn <- function(x, digits=0) {format(round(x, digits), nsmall=digits, trim=T)}
round.fn <- function(x, round=0,trim=T,big.mark=",",significant=NULL, justify=c("left","right","centre","none"),width=NULL,na.encode=T,scientific=NA,big.interval=3L,small.mark="",small.interval=5L,decimal.mark=getOption("OutDec"),zero.print=NULL,drop0trailing=F,...) {
  .x <- round(x, digits=round)
  .x <- format(.x, nsmall=round, trim=trim, big.mark=big.mark, digits=significant,justify=justify,width=width,na.encode=na.encode,scientific=scientific,big.interval=big.interval,small.mark=small.mark,small.interval=small.interval,decimal.mark=decimal.mark,zero.print=zero.print,drop0trailing=drop0trailing,...)
  return(.x)
}


## Calculate day of year, but 29th Feb assigned 60, 1 Mar assigned 61 regardless of whether leap year and then continue until 366. This ensures leap years are done on same order as non-leap years
# Requires lubridate
# Source: https://stackoverflow.com/questions/54946395/counting-days-of-the-year-with-leap-years
leap.ever.year.fn <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}

