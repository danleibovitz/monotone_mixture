# Wishlist
# - For cox_reg model, include "+ cluster(id)" option in formula? I imagine this would allow soft groupings (correlation) 
# of ids within component models of the mixture, whereas what we really want is the "| id" grouping of flexmix, 
# since a single id should not be able to belong to multiple clusters. That said, perhaps a "cluster(family)" argument
# might make sense, when, e.g., members of a family can belong to distinct clusters, but members of the same family within
# a cluster should have their family status incorporated in the model estimate

# allow slots defined for numeric to accept NULL
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))

# Define new classes
setClass(
  "FLX_Cox_component",
  contains="FLXcomponent",
  # TODO what slots does Cox_component need?
  # TODO take either start/stop, or time, but not both
  slots=c(start="numericOrNULL", 
          stop="numericOrNULL", 
          time="numericOrNULL",
          status="numericOrNULL")
) 

# Define FLXM_Cox
setClass("FLXM_Cox",
         contains = "FLXM",
         # TODO what slots does FLXM_Cox need?
         # TODO take either start/stop, or time, but not both
         slots = c(start="numericOrNULL", 
                   stop="numericOrNULL", 
                   time="numericOrNULL",
                   status="numericOrNULL"))


# definition of monotone regression model.

# TODO include OFFSET as optional argument?
cox_reg <- function (formula = .~., start=NULL, stop=NULL, time=NULL, status=NULL) {
  
  retval <- new("FLXM_Cox", weighted = TRUE,
                formula = formula,
                name = "Cox PH regression",
                start= start,
                stop= stop,
                time= time,
                status= status) 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  # fit must have attributes: coef, sigma, cov, df, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  # ... all must be defined by fit() function
  retval@defineComponent <- function(fit, ...) {
    # @logLik: A function(x,y) returning the log-likelihood for observations in matrices x and y
    logLik <- function(x, y) { 
      dnorm(y, mean=predict(x, ...), sd=fit$sigma, log=TRUE)
      
      newfit <- coxph(retval@formula, init = fit$coefficients, control=coxph.control(iter.max=0), 
                     data = data.frame(x,y))
      logLik(newfit)
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
    new("FLX_Cox_component", parameters =
          list(coef = fit$coef, sigma = fit$sigma),
        df = fit$df, logLik = logLik, predict = predict,
        mon_inc_index = fit$mon_inc_index,
        mon_dec_index = fit$mon_dec_index)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  retval@fit <- function(retval, component, ...) {
    
    # TODO provide coxph an 'init' argument from last component's coefficients
    fit <- coxph(formula = retval@formula, 
                 data = df, weights = w, ...)
    
    retval@defineComponent(fit, ...)
  }
  
  # retval@fit <- function(x, y, w, component, start = retval@start, 
  #                        stop = retval@stop, status = retval@status, ...) {
  #   
  #   df <- data.frame(y = y, x = x, w = w)
  #   names(df)[start] <- "start"
  #   names(df)[stop] <- "stop"
  #   df$y <- "status"
  #   
  #   # TODO provide coxph an 'init' argument from last component's coefficients
  #   fit <- coxph(formula = Surv(start, stop, status) ~ . -w -start -stop -start, 
  #                data = df, weights = w, ...)
  #   
  #   retval@defineComponent(fit, ...)
  # }
  
  
  
  
  retval 
}


