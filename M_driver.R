
source("part_fit.R")


# allow slots defined for numeric to accept NULL
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))

# Define new classes
setClass(
  "FLX_component_plus",
  contains="FLXcomponent",
  # TODO allow mon_index to take either numeric or NULL
  slots=c(mon_inc_index="numericOrNULL", mon_dec_index="numericOrNULL")
) 

## TODO DEFINE FLXM_monoreg CORRECTLY!!!### give it mon_index attributes
setClass("FLXM_monoreg",
         representation(y="matrix",
                        dist="character"),
         contains = "FLXM")


# definition of monotone regression model.

# TODO include OFFSET as optional argument?
mono_reg <- function (formula = .~., ...) {

  retval <- new("FLXM_monoreg", weighted = TRUE,
                # TODO is dist param necessary?
                formula = formula,
                name = "partially linear monotonic regression") 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  # fit must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  # ... all must be defined by fit() function
  # TODO are extra arguments ... to defineComponent() necessary?
  retval@defineComponent <- function(fit, ...) {
                  # @logLik: A function(x,y) returning the log-likelihood for observations in matrices x and y
                  logLik <- function(x, y) { 
                    dnorm(y, mean=predict(x, ...), sd=fit$sigma, log=TRUE)
                  }
                  # @predict: A function(x) predicting y given x. 
                  # TODO x must be partitioned into linear and monotone covars
                  predict <- function(x) {
                    inc_ind <- fit$mon_inc_index
                    dec_ind <- fit$mon_dec_index
                    
                    # TODO give appropriate prediction whether partial linear model is increasing, 
                    # decreasing, with monotone index at 1 or elsewhere, and whether or not there are
                    # linear components
                    
                    p <-  get_pred(fit$fitted_pava, x[,c(inc_ind, dec_ind)])
                    if(!is.null(fit$coef)){
                      p <- p + (x[,-c(inc_ind, dec_ind)] %*% fit$coef)
                    }
                    p
                  }
                  # TODO what FLXcomponent object NEED?
                  ret <- new("FLXcomponent", parameters =
                        list(coef = fit$coef, sigma = fit$sigma),
                      df = fit$df, logLik = logLik, predict = predict)
                  
                  # ret is cast as subclass object
                  ret <- as(ret,"FLX_plus")
                  ret@mon_inc_index <- fit$mon_inc_index
                  ret@mon_dec_index <- fit$mon_dec_index
                  
                  ret
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  #TODO inputs are x = covariates, y = DV, w = weights. Must have specification of covars with monotone relationship
  retval@fit <- function(x, y, w, component,...) {
    print(list(...))
                  fit <- part_fit(x, y, w, component, ...)
                  
                  retval@defineComponent(fit, ...)
                  }
  retval 
  }


