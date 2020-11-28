# define fitting function for partial linear model with 1 monotone-related variable


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
  x <- as.matrix(x)
  y <- as.matrix(y)
  wates <- as.matrix(wates)
  
  # assume that monotone variable is first column in x, unless specified otherwise
  dotarg = list(...)
  if("mon_index" %in% names(dotarg)){
    ind <- mon_index
  } 
  else{
    ind <- 1
  }
  # TODO option for monotonic decreasing regression
  
  # throw error if the number of indices exceeds columns of x
  if(length(ind) > ncol(x)) stop("Number of proposed monotonic relationships exceeds columns of x.")
  
  # TODO option for fit with no linear independent components
  # TODO option for fit with no linear independent components and multiple monotone components
  if(length(ind) == ncol(x)){
    
    yhat <- monoreg(x = x[ind], y = y, w = wates[ind])
    
    mod$para <- NULL
    mod$fitted_pava <- yhat
    # TODO how does return object need to be formatted?
    return(mod)
  }
  
  # for starting values, fit a regular lm
  fit <- lm.wfit(x=x, y=y, w=wates)
  betas <- coef(fit)[-ind]

  # set maximum iterations for convergence
  if("max_iter" %in% names(dotarg)){
    if(max_iter < 1) stop("max_iter must be positive")
    maxiter <- max_iter
  }
  else{ 
    maxiter <- 100
  }
    
  # iterate between pava and linear model
  # TODO set while loop condition(s). Get appropriate measure of coefficient change
  iter <- 0
  delta <- 10
  while(delta > 1 & iter < maxiter){
    # yhat <- pava( (y - x[-ind] %*% betas), long.out = T, stepfun = T)
    
    yhat <- monoreg(x = x[ind], y = (y - x[-ind] %*% betas), w = wates[ind])
    
    old_betas <- betas    # save old betas for distance calculation
    #TODO how to get monoreg sorted values?
    betas <- coef(lm.wfit(x=x[-ind], y= (y - yhat$sorted_vals), w=wates[-ind]))
    
    # TODO quantify change in yhat vals and beta vals
    # get euclidian distance
    delta <- dist(rbind(as.vector(betas), as.vector(old_betas)))

    iter <- iter + 1    # iterate maxiter
  }
  
  # TODO how does return object need to be formatted?
  # mod must have coefficients attribute and monoreg fit attribute
  mod$para <- betas
  mod$fitted_pava <- yhat
  
  return(mod)
}



# Example 1.1.1. (dental study) from Robertson, Wright and Dykstra (1988)
age = c(14, 13, 7, 8, 8, 9, 10, 10, 11, 11.5, 12)
size = c(23.5, 25, 21, 23.5, 23, 24, 21, 25, 21.5, 22, 19)
weights = c(0.8, 0.14, 0.5, 1, 1, 0.2, 1, 1, 0.9, 0.5, 0.5)
mr = monoreg(age, size, w = weights)

# sorted x values
mr$x # 8 10 12 14
# weights and merged y values
mr$w  # 3 3 3 2
mr$y #  22.50000 23.33333 20.83333 24.25000
# fitted y values
mr$yf # 22.22222 22.22222 22.22222 24.25000
fitted(mr)
var(residuals(mr))
mr$sd <- var(residuals(mr))
mr$sd




plot(mr)  # this shows the averaged data points
points(age, size, pch=2)  # add original data points

y <- (1:10) + rnorm(10, sd = 3)
x <- 1:10 

mr <- monoreg(x, y)
plot(mr)
points(x, y, pch=2)  # add original data points

ystar <- pava(y, long.out = T, stepfun = T)
plot(y)
lines(ystar$y,type='s')
ystar$h(10)



