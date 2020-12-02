# define fitting function for partial linear model with 1 monotone-related variable


# TODO define CPAV algorithm for addittive, iterative PAVA
# TODO replace monoreg() with cpav()
cpav <- function(){
  
}


# first, define function for obtaining f(x_new) for monotone regression f()
# TODO enable mr_obj to hold multiple monoreg fits
get_pred <- function(mr_obj, xval){
  xval <- as.vector(xval)
  last <- length(mr_obj$x)
  
  mr_obj$yf[sapply(xval, function(z)
    ifelse( z < mr_obj$x[1], 1,
            ifelse(z >= tail(mr_obj$x, n=1), last, 
                   which.min(mr_obj$x <= z)-1 ))
    )]
}


#  ifelse( xval < mr_obj$x[1], mr_obj$yf[1],
#          ifelse(xval > tail(mr_obj$x, n=1), tail(mr_obj$yf, n=1), 
#                 mr_obj$yf[sapply(testy, function(x) which.min(mr$x <= x)-1)]))
#}
#
#sapply(testy, function(x) which.min(mr$x <= x)-1)
#
#  if(xval < mr_obj$x[1]){ # xval is lower than PAVA range
#    yval <- mr_obj$yf[1]
#  }
#  else{
#    if(xval > tail(mr_obj$x, n=1)){ # xval is higher than PAVA range
#      yval <- tail(mr_obj$yf, n=1)
#    }
#    else{ # xval is within PAVA range
#      yval <- mr_obj$yf[Position(function(x) x >= xval, mr_obj$x)- 1]
#    }
#  }
#  
#  yval
#}


# define partial linear regression of y on x with weights w
part_fit <- function(x, y, wates, ...){
  
  # TODO cast y and wates to matrices ?
  # TODO correct behaviour for if x is ONLY a vector
  x <- as.matrix(x)
  
  # TODO make sure y and wates is not multivariate
  if(length(y) != dim(x)[1] | length(y) != length(wates)) stop("Inputs are not of the same dimension!")
  
  # TODO adapt to include monotonic _decreasing_ regression
  # assume that monotone variable is first column in x and increasing, unless specified otherwise
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
  
  # throw error if the number of indices exceeds columns of x
  if(length(c(inc_ind, dec_ind)) > ncol(x)) stop("Number of proposed monotonic relationships exceeds columns of x.")
  
  # option for fit with no linear independent components
  # TODO option for fit with no linear independent components and multiple monotone components
  if(length(c(inc_ind, dec_ind)) == ncol(x)){
    
    yhat <- monoreg(x = x[,inc_ind], y = y, w = wates)
    
    # get residuals of model
    resids <- y - get_pred(yhat, x[,c(inc_ind, dec_ind)])
    
    # mod must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
    # may have mon_inc_index and mon_dec_index attributes
    mod <- list(coef = NULL, fitted_pava = NULL, sigma = NULL, df = NULL,
                mon_inc_index = NULL, mon_dec_index = NULL, iterations = NULL)
  
    mod$coef <- NULL
    mod$fitted_pava <- yhat
    mod$mon_inc_index <- inc_ind
    mod$mon_dec_index <- dec_ind
    # TODO sigma is the sqrt of __ divided by (nrow(x) - model$rank). For single monoreg, rank is 1...
    # ... but update when implementing CPAV
    # TODO ask matthias : rank of input matrix? or model matrix? different, depending on dummy coding, etc.
    mod$sigma <- sqrt(sum(wates * (resids)^2 /
                                     mean(wates))/ (nrow(x)-qr(x)$rank))
    mod$df <- ncol(x)+1
    
    return(mod)
  }
  else{
    # for starting values, fit a regular lm
    fit <- lm.wfit(x=x, y=y, w=wates)
    betas <- coef(fit)[-c(inc_ind, dec_ind)]
  
    # set maximum iterations for convergence
    if("max_iter" %in% names(dotarg)){
      if(max_iter < 1) stop("max_iter must be positive")
      maxiter <- max_iter
    }
    else{ 
      maxiter <- 10000
    }
      
    # set while loop initial values
    iter <- 0
    delta <- 10
    # iterate between pava and linear model
    # TODO set while loop condition(s). Get appropriate measure of coefficient change
    while(delta > 1e-12 & iter < maxiter){
  
      yhat <- monoreg(x = x[,inc_ind], y = (y - x[,-inc_ind] %*% betas), w = wates)
      
      old_betas <- betas    # save old betas for distance calculation
      # to retrieve old ordering of y for fitted values, we use y[match(x, sorted_x)]
      betas <- coef(lm.wfit(x=x[,-inc_ind], y= (y - yhat$yf[match(x[,inc_ind], yhat$x)] ), w=wates))
      
      # TODO quantify change in yhat vals and beta vals
      # get euclidian distance between betas transformed into unit vectors
      delta <- dist(rbind( as.vector(betas)/norm(as.vector(betas), type ="2"), 
                           as.vector(old_betas)/norm(as.vector(old_betas), type="2")
                           ))
  
      iter <- iter + 1    # iterate maxiter
    }
  }
  
  # get residuals of model
  resids <- y - (get_pred(yhat, x[,c(inc_ind, dec_ind)]) + (x[,-c(inc_ind, dec_ind)] %*% betas))
  
  
  # mod must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  mod <- list(coef = NULL, fitted_pava = NULL, sigma = NULL, df = NULL,
              mon_inc_index = NULL, mon_dec_index = NULL, iterations = NULL)
  
  mod$coef <- betas
  mod$fitted_pava <- yhat
  # TODO keep iterations attribute? 
  mod$iterations <- iter
  mod$mon_inc_index <- inc_ind
  mod$mon_dec_index <- dec_ind
  # TODO sigma is the sqrt of __ divided by (nrow(x) - model$rank). ...
  # ... but update when implementing CPAV
  # TODO ask matthias : rank of input matrix? or model matrix? different, depending on dummy coding, etc.
  mod$sigma <- sqrt(sum(wates * (resids)^2 / mean(wates))/ (nrow(x)-qr(x)$rank))
  mod$df <- ncol(x)+1
  
  
  return(mod)
}






