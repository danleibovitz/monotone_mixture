# define fitting function for partial linear model with arbitrary monotone-constrained component

# Wishlist:
# - each part_fit() should return a part_fit object, such that plot(part_fit) is defined
# - provide part_fit() return object with a plot() method that plots monotone component(s)
# - select monotone component by name instead of index


# TODO define CPAV algorithm for addittive, iterative PAVA
# TODO replace monoreg() with cpav() in part_fit()
cpav <- function(x_mat, y, weights, inc_index, dec_index){
  
  joint_ind <- c(inc_index, dec_index)
  
  if(!is.matrix(x_mat)) stop("x_mat is not of class matrix, and will be rejected by lm.wfit")
  if(any(weights == 0)) stop("monoreg() cannot take weights of 0!")
  if(length(y) != length(weights)) stop("The dimension of the inputs is not 
                                        equal to the dimension of the weights")
  
  # if there is only 1 monotone component, apply ordinary monotone regression
  if(length(joint_ind) == 1){
    if(length(inc_index) == 1){ # the component is monotone increasing
      return(
        monoreg(x = x_mat[,inc_index], y = y, w = weights)
      )
    }
    else{ # the component is monotone decreasing
      return(
        monoreg(x = x_mat[,dec_index], y = y, w = weights, type = "antitonic")
      )
    }
  }
  
  else{ # the monotone components are multiple, so continue with cyclic algorithm
    
    # fit ordinary lm on x_mat and y
    start_betas <- coef(lm.wfit(x=x_mat[,joint_ind], y=y, w=weights))
    
    # set initial monotone reg estimates by calling each monoreg() against y - lm.predict(all other vars)
    mr_fits <- as.vector(sapply(1:length(joint_ind), function(i) 
      if(joint_ind[i] %in% inc_index){
        # i apologize to anyone trying to read this line, but think: the columns of the x_matrix
        # indicated by join_ind, except the value of joint_ind at the ith place in joint_ind
        monoreg(x = x_mat[,joint_ind[i]], 
                y = (y - x_mat[,joint_ind[-joint_ind[i]]] %*% start_betas[-i] ), w = weights, type = "isotonic")
        }
      else if(joint_ind[i] %in% dec_index){
        monoreg(x = x_mat[,joint_ind[i]], 
                y = (y - x_mat[,joint_ind[-joint_ind[i]]] %*% start_betas[-i] ), w = weights, type = "antitonic")
      }))
    
    return(mr_fits)
    # iterate through mr_fits
  }
  
  
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
# TODO max_iter must be passed through M_driver just like mon_inc_index, or 
# else user will not be able to specify the max_iter param
part_fit <- function(x, y, wates = NULL, mon_inc_index=NULL, mon_dec_index=NULL, max_iter=NULL, 
                     component = NULL, ...){
  

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
    inc_ind <- component$mon_inc_index
    dec_ind <- component$mon_dec_index
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
 
  # throw warning if there are duplicates in inc_ind or dec_ind, and then remove
  if(anyDuplicated(inc_ind) | anyDuplicated(dec_ind)){
    warning("There are duplicate index instructions; Duplicates are being removed.")
    inc_ind <- unique(inc_ind)
    dec_ind <- unique(dec_ind)
  }
  
  # throw error if indices overlap
  if(length(intersect(inc_ind, dec_ind)) > 0) stop("At least one variable was marked as BOTH
                                               monotone increasing and monotone decreasing.")
  
  # throw error if indices are not integers
  if(!is.null(inc_ind)){
    if(any(inc_ind != as.integer(inc_ind))) stop("Monotone increasing indices are not integers.")
  }
  if(!is.null(dec_ind)){
    if(any(dec_ind != as.integer(dec_ind))) stop("Monotone decreasing indices are not integers.")
  }
  
  # throw error if indices are not positive 
  if(any(c(inc_ind, dec_ind) < 1)) stop("all monotone component indices must be positive")
  
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
    # TODO print(c(class(inc_ind), class(dec_ind)))
    betas <- coef(fit)[-c(inc_ind, dec_ind)]
  
    # set maximum iterations for convergence
    if(!is.null(max_iter) & !is.list(max_iter)){
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
  
      yhat <- monoreg(x = x[wates != 0,inc_ind], 
                      y = (y[wates != 0] - x[wates != 0,-inc_ind] %*% betas), w = wates[wates != 0])
      
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
  # TODO print(c(class(inc_ind), class(dec_ind)))
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






