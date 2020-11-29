# define fitting function for partial linear model with 1 monotone-related variable


# TODO define CPAV algorithm for addittive, iterative PAVA
# TODO replace monoreg() with cpav()
cpav <- function(){
  
}


# first, define function for obtaining f(x_new) for monotone regression f()
# TODO vectorize get_pred to accept vector arguments for xval
# TODO enable mr_obj to hold multiple monoreg fits
get_pred <- function(mr_obj, xval){
  if(xval < mr_obj$x[1]){ # xval is lower than PAVA range
    yval <- mr_obj$yf[1]
  }
  else{
    if(xval > tail(mr_obj$x, n=1)){ # xval is higher than PAVA range
      yval <- tail(mr_obj$yf, n=1)
    }
    else{ # xval is within PAVA range
      yval <- mr_obj$yf[Position(function(x) x > xval, mr_obj$x)- 1]
    }
  }
  
  yval
}


# define partial linear regression of y on x with weights w
part_fit <- function(x, y, wates, ...){
  
  # cast x, y and wates to matrices
  #x <- as.matrix(x)
  #y <- as.matrix(y)
  #wates <- as.matrix(wates)
  
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
    
    yhat <- monoreg(x = x[inc_ind], y = y, w = wates[inc_ind])
    
    mod$para <- NULL
    mod$fitted_pava <- yhat

    return(mod)
  }
  
  # for starting values, fit a regular lm
  fit <- lm.wfit(x=x, y=y, w=wates)
  betas <- coef(fit)[-c(inc_ind, dec_ind)]

  # set maximum iterations for convergence
  if("max_iter" %in% names(dotarg)){
    if(max_iter < 1) stop("max_iter must be positive")
    maxiter <- max_iter
  }
  else{ 
    maxiter <- 100
  }
    
  # set while loop initial values
  iter <- 0
  delta <- 10
  # iterate between pava and linear model
  # TODO set while loop condition(s). Get appropriate measure of coefficient change
  while(delta > 1e-3 & iter < maxiter){

    yhat <- monoreg(x = x[inc_ind], y = (y - x[-inc_ind] %*% betas), w = wates[inc_ind])
    
    old_betas <- betas    # save old betas for distance calculation
    # to retrieve old ordering of y for fitted values, we use y[match(x, sorted_x)]
    betas <- coef(lm.wfit(x=x[-inc_ind], y= (y - yhat$yf[match(x[inc_ind], yhat$x)] ), w=wates[-inc_ind]))
    
    # TODO quantify change in yhat vals and beta vals
    # get euclidian distance between betas transformed into unit vectors
    delta <- dist(rbind( as.vector(betas)/norm(as.vector(betas), type ="2"), 
                         as.vector(old_betas)/norm(as.vector(old_betas), type="2")
                         ))

    iter <- iter + 1    # iterate maxiter
  }
  
  # TODO how does return object need to be formatted?
  # mod must have coefficients attribute and monoreg fit attribute
  mod$para <- betas
  mod$fitted_pava <- yhat
  
  return(mod)
}


X <- cbind(
  sample(seq(from = -50, to = 50), size = 200, replace = TRUE),
  sample(seq(from = -100, to = 100), size = 200, replace = TRUE),
  sample(seq(from = -100, to = 100), size = 200, replace = TRUE),
  sample(seq(from = -100, to = 100), size = 200, replace = TRUE),
  sample(seq(from = -100, to = 100), size = 200, replace = TRUE)
)

W <- rep(1, 200)

Y <- (X[,1])^3 + 3*X[,2] + 2*X[,3] - 4*X[,4] + X[,5] + rnorm(200, 0, 3)

pseudo_df <- data.frame(Y, X, W)





