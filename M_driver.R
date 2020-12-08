
source("part_fit.R")


# allow slots defined for numeric to accept NULL
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))

# Define new classes
setClass(
  "FLX_component_plus",
  contains="FLXcomponent",
  # allow mon_index to take either numeric or NULL
  slots=c(mon_inc_index="numericOrNULL", mon_dec_index="numericOrNULL")
) 

# Define FLXM_monoreg 
setClass("FLXM_monoreg",
         # TODO what does FLXM_monoreg need to inherit?
         contains = "FLXM",
         slots = c(mon_inc_index="numericOrNULL", mon_dec_index="numericOrNULL"))


# definition of monotone regression model.

# TODO include OFFSET as optional argument?
mono_reg <- function (formula = .~., mon_inc_index=NULL, mon_dec_index=NULL) {

  retval <- new("FLXM_monoreg", weighted = TRUE,
                formula = formula,
                name = "partially linear monotonic regression",
                mon_inc_index= mon_inc_index,
                mon_dec_index= mon_dec_index) 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  # fit must have attributes: coef, sigma, cov, df, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  # ... all must be defined by fit() function
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
                  # return new FLX_component_plus object
                  new("FLX_component_plus", parameters =
                        list(coef = fit$coef, sigma = fit$sigma),
                      df = fit$df, logLik = logLik, predict = predict,
                      mon_inc_index = fit$mon_inc_index,
                      mon_dec_index = fit$mon_dec_index)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  retval@fit <- function(x, y, w, component, mon_inc_index = retval@mon_inc_index, 
                         mon_dec_index = retval@mon_dec_index, ...) {
    
                  fit <- part_fit(x, y, w, component, mon_inc_index=mon_inc_index, 
                                  mon_dec_index=mon_dec_index, ...)
                  
                  retval@defineComponent(fit, ...)
                  }
  retval 
  }


