
source("part_fit.R")

# definition of monotone regression model.

mono_reg <- function (formula = .~., diagonal = TRUE) {
  retval <- new("FLXMC", weighted = TRUE,
                # TODO is dist param necessary?
                formula = formula, dist = "mvnorm",
                name = "my model-based clustering") 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  # fit must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  # ... all must be defined by fit() function
  retval@defineComponent <- function(fit, fitted_pava) {
                  # @logLik: A function(x,y) returning the log-likelihood for observations in matrices x and y
                  logLik <- function(x, y) { 
                    dnorm(y, mean=predict(x, ...), sd=fit$sigma, log=TRUE)
                  }
                  # @predict: A function(x) predicting y given x. 
                  # TODO x must be partitioned into linear and monotone covars
                  predict <- function(x, ...) {
                    dotarg = list(...)
                    if("mon_inc_index" %in% names(dotarg)){
                      inc_ind <- mon_inc_index
                    } 
                    else{
                      inc_ind <- 1
                    }
                    if("mon_dec_index" %in% names(dotarg)){
                      dec_ind <- mon_dec_index
                    } 
                    else{
                      dec_ind <- NULL
                    }
                    
                    p <- (x %*% fit$coef) + get_pred(fitted_pava, x_mon)
                    p
                  }
                  new("FLXcomponent", parameters =
                        list(center = fit$center, cov = fit$cov),
                      df = fit$df, logLik = logLik, predict = predict)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  #TODO inputs are x = covariates, y = DV, w = weights. Must have specification of covars with monotone relationship
  retval@fit <- function(x, y, w, ...) {
                  fit <- part_fit(x, y, w, ...)
                  
                  retval@defineComponent(fit, ...) 
                  }
  retval 
  }


