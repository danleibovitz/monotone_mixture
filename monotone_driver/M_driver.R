

# ## Wishlist:
# - somehow pass design matrix names from flexmix() call back to the construction of monoreg() so that
# names can be interpreted as indices. Perhaps, override "FLXgetModelMatrix" method and add names(design_matrix)
# slot? 
# flexmix formula as y~x|g where g is the grouping variable, but what about a mixed-model component formula?
# - allow factors for non-monotone components of part_fit, and disallow factors for monotone components
# - mono_reg cannot accept ANY na values...
# - returned model should have confidence intervals, with different options for CI 
# construction. Plotting should (automatically?) include CIs


# source("monotone_driver/part_fit.R")


# allow slots defined for numeric to accept NULL
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))
setOldClass("monoreg")
setClassUnion("matrixOrMonoreg", members = c("matrix", "monoreg"))

# Define new classes
setClass(
  "FLX_monoreg_component",
  contains="FLXcomponent",
  # allow mon_index to take either numeric or NULL
  slots=c(mon_inc_index="numericOrNULL", 
          mon_dec_index="numericOrNULL",
          mon_obj="matrix",
          mon_inc_names="characterOrNULL",
          mon_dec_names="characterOrNULL"
          )
) 

# Define FLXM_monoreg 
setClass("FLXM_monoreg",
         # TODO what does FLXM_monoreg need to inherit?
         contains = "FLXM",
         slots = c(mon_inc_index="numericOrNULL", 
                   mon_dec_index="numericOrNULL",
                   mon_inc_names="characterOrNULL",
                   mon_dec_names="characterOrNULL"))






# definition of monotone regression model.

# TODO include OFFSET as optional argument?
mono_reg <- function (formula = .~., mon_inc_names = NULL, 
                      mon_dec_names = NULL, mon_inc_index=NULL, mon_dec_index=NULL, ...) {

  # only names or indices can be indicated, not both
  if((!is.null(mon_inc_names)|!is.null(mon_dec_names)) &
     (!is.null(mon_inc_index)|!is.null(mon_dec_index))) stop("mono_reg() can accept either monotone
                                                             names or indices can be chosen, but not both.")
  
  retval <- new("FLXM_monoreg", weighted = TRUE,
                formula = formula,
                name = "partially linear monotonic regression",
                mon_inc_index= sort(mon_inc_index),
                mon_dec_index= sort(mon_dec_index), 
                mon_inc_names= mon_inc_names,
                mon_dec_names= mon_dec_names) 
  
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
                      p <- p + (as.matrix(x[,-c(inc_ind, dec_ind)]) %*% fit$coef)
                    }
                    p
                  }
                  # return new FLX_monoreg_component object
                  new("FLX_monoreg_component", parameters =
                        list(coef = fit$coef, sigma = fit$sigma, mon_obj = fit$fitted_pava),
                      df = fit$df, logLik = logLik, predict = predict,
                      mon_inc_index = fit$mon_inc_index,
                      mon_dec_index = fit$mon_dec_index,
                      # mon_obj = fit$fitted_pava,
                      mon_inc_names = fit$mon_inc_names,
                      mon_dec_names = fit$mon_dec_names)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  retval@fit <- function(x, y, w, component, mon_inc_index = retval@mon_inc_index, 
                         mon_dec_index = retval@mon_dec_index, 
                         mon_inc_names = retval@mon_inc_names,
                         mon_dec_names = retval@mon_dec_names, ...) {
    
                  
                  if(is.null(mon_inc_index) & is.null(mon_dec_index)){
                    
                    # if not all monotone names are in the design matrix, stop & print the name that is missing
                    if(!all(c(mon_inc_names, mon_dec_names) %in% colnames(x))){
                      stop(paste(setdiff(c(mon_inc_names, mon_dec_names), colnames(x)),
                                 "could not be found in the model matrix. Check your spelling."))
                    } 
                    # Discover correct monotone indices
                    if(any(colnames(x) %in% mon_inc_names)){
                      mon_inc_index <- which(colnames(x) %in% mon_inc_names)
                    }
                    if(any(colnames(x) %in% mon_dec_names)){
                      mon_dec_index <- which(colnames(x) %in% mon_dec_names)
                    }
                  }
                  if(is.null(mon_inc_names) & is.null(mon_dec_names)){
                    # Discover correct monotone names
                    mon_inc_names <- colnames(x)[sort(mon_inc_index)]
                    mon_dec_names <- colnames(x)[sort(mon_dec_index)]
                  }
    
                  # TODO exclude monotone fits of factors
                  # if(any(apply(x, 2, function(x) is.factor(x)))) stop("x cannot have factor columns as monotone components")

                  fit <- part_fit(x, y, w, component, mon_inc_index=mon_inc_index, 
                                  mon_dec_index=mon_dec_index, ...)
                  
                  retval@defineComponent(fit, ...)
                  }
  retval 
  }


