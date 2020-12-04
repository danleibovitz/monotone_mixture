# define fitting function for partial linear model with 1 monotone-related variable

# Wishlist:
# provide part_fit() return object with a print() method that plots monotone component(s)
# select monotone component by name instead of index

# TODO define CPAV algorithm for addittive, iterative PAVA
# TODO replace monoreg() with cpav()
# TODO mon_inc_index argument not working...
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


# define partial linear regression of y on x with weights w
# inputs are: x, y, wates, mon_inc_index, mon_dec_index, max_iter
part_fit <- function(x, y, wates = NULL, mon_inc_index=NULL, mon_dec_index=NULL, max_iter=NULL, 
                     component = NULL, ...){
  
  print(c(class(mon_inc_index), class(mon_dec_index), class(component)))
  
  # TODO cast y and wates to matrices ?
  # TODO correct behaviour for if x is ONLY a vector
  x <- as.matrix(x)
  
  # set default weights
  if(is.null(wates)) wates <- rep(1, length(y))
  
  # TODO make sure y and wates is not multivariate
  if(length(y) != dim(x)[1] | length(y) != length(wates)) stop("Inputs are not of the same dimension!")
  
  # TODO adapt to include monotonic _decreasing_ regression
  
  # take monotone indices of previous component
  if(!is.null(component)){
    inc_ind <- dotarg$component$mon_inc_index
    dec_ind <- dotarg$component$mon_dec_index
  }
  else{
    # assume that monotone variable is first column in x and increasing, unless specified otherwise
    if(!is.null(mon_inc_index)){
      inc_ind <- mon_inc_index
    } 
    else{
      inc_ind <- 1
    }
    if(!is.null(mon_dec_index)){
      dec_ind <- mon_dec_index
    } 
    else{
      dec_ind <- NULL
    }
  }
  
  # throw error if the number of indices exceeds columns of x
  if(length(c(inc_ind, dec_ind)) > ncol(x)) stop("Number of proposed monotonic relationships exceeds columns of x.")
  
  # option for fit with no linear independent components
  # TODO option for fit with no linear independent components and multiple monotone components
  if(length(c(inc_ind, dec_ind)) == ncol(x)){
    
    yhat <- monoreg(x = x[wates != 0,inc_ind], y = y[wates != 0], w = wates[wates != 0])
    
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
    print(c(class(inc_ind), class(dec_ind)))
    betas <- coef(fit)[-c(inc_ind, dec_ind)]
  
    # set maximum iterations for convergence
    if(!is.null(max_iter)){
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
  
      yhat <- monoreg(x = x[wates != 0,inc_ind], y = (y[wates != 0] - x[wates != 0,-inc_ind] %*% betas), w = wates[wates != 0])
      
      old_betas <- betas    # save old betas for distance calculation
      # to retrieve old ordering of y for fitted values, we use y[match(x, sorted_x)]
      betas <- coef(lm.wfit(x=x[,-inc_ind], y= (y - get_pred(yhat, x[,inc_ind]) ), w=wates))

      # TODO quantify change in yhat vals and beta vals
      # get euclidian distance between betas transformed into unit vectors
      delta <- dist(rbind( as.vector(betas)/norm(as.vector(betas), type ="2"), 
                           as.vector(old_betas)/norm(as.vector(old_betas), type="2")
                           ))
  
      iter <- iter + 1    # iterate maxiter
    }
  }
  
  # get residuals of model
  print(c(class(inc_ind), class(dec_ind)))
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






