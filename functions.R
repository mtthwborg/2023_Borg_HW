
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
