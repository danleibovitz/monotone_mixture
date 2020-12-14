

# allow slots defined for numeric to accept NULL
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))

# Define new classes
setClass(
  "FLX_Cox_component",
  contains="FLXcomponent",
  # allow mon_index to take either numeric or NULL
  # TODO what slots does Cox_component need?
  slots=c(mon_inc_index="numericOrNULL", mon_dec_index="numericOrNULL")
) 

# Define FLXM_Cox
setClass("FLXM_Cox",
         contains = "FLXM",
         # TODO what slots does FLXM_Cox need?
         slots = c(elap_time="numericOrNULL", status="numericOrNULL"))


# definition of monotone regression model.

# TODO include OFFSET as optional argument?
mono_reg <- function (formula = .~., elap_time=NULL, status=NULL) {
  
  retval <- new("FLXM_Cox", weighted = TRUE,
                formula = formula,
                name = "Cox PH regression",
                elap_time= elap_time,
                status= status) 
  
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
    # return new FLX_monoreg_component object
    new("FLX_monoreg_component", parameters =
          list(coef = fit$coef, sigma = fit$sigma),
        df = fit$df, logLik = logLik, predict = predict,
        mon_inc_index = fit$mon_inc_index,
        mon_dec_index = fit$mon_dec_index)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  retval@fit <- function(x, y, w, component, mon_inc_index = retval@mon_inc_index, 
                         mon_dec_index = retval@mon_dec_index, ...) {
    
    fit <- coxph(x, y, w, component, mon_inc_index=mon_inc_index, 
                    mon_dec_index=mon_dec_index, ...)
    
    retval@defineComponent(fit, ...)
  }
  retval 
}


