# define fitting function for partial linear model with arbitrary monotone-constrained component

# Wishlist:
# - select monotone component by name instead of index
# - can part_fit take factors as monotone components?
# - Confidence intervals for linear components
# - plot method for linear components


# libraries
library(gridExtra)
library(dplyr)

# TODO apparently, cpav with any more than 2 components overfits infintely if the error of the
# original Y construction is small....
# inc_index and dec_index are the indices of x_mat which are supposed to have an isotonic and antitonic
# relationship (respectively) with dependent variable y
cpav <- function(x_mat, y, weights, inc_index=NULL, dec_index=NULL, max_iters_cpav=NULL, max_delta_cpav=NULL){
  
  joint_ind <- c(inc_index, dec_index)
  
  if(!is.matrix(x_mat)) stop("x_mat is not of class matrix, and will be rejected by lm.wfit")
  if(any(weights == 0)) stop("monoreg(), and therefore cpav(), cannot take weights of 0!")
  if(length(y) != length(weights) | length(y) != dim(x_mat)[1]) stop("The dimension of the inputs is not 
                                        equal to the dimension of the weights")
  
  # if there is only 1 monotone component, apply ordinary monotone regression
  if(length(joint_ind) == 1){
    if(length(inc_index) == 1){ # the component is monotone increasing
      return( # cast the monoreg object as a matrix, with all attributes as rows in the first column
        matrix(suppressWarnings(monoreg(x = x_mat[,inc_index], y = y, w = weights)), 
               dimnames = list(c("x", "y", "w", "yf", "type", "call")))
      )
    }
    else{ # the component is monotone decreasing
      return( # cast the monoreg object as a matrix, with all attributes as rows in the first column
        
        matrix(suppressWarnings(monoreg(x = x_mat[,dec_index], y = y, w = weights, type = "antitonic")), 
               dimnames = list(c("x", "y", "w", "yf", "type", "call")))
      )
    }
  }
  
  else{ # the monotone components are multiple, so continue with cyclic algorithm
    
    # fit ordinary lm on x_mat and y
    start_betas <- coef(lm.wfit(x=x_mat[,joint_ind], y=y, w=weights))
    
    # set initial monotone reg estimates by calling each monoreg() against y - lm.predict(all other vars)
    
    mr_fits <- sapply(1:length(joint_ind), function(i) 
      if(joint_ind[i] %in% inc_index){
        # i apologize to anyone trying to read this line, but think: the columns of the x_matrix
        # indicated by join_ind, except the value of joint_ind at the ith place in joint_ind
        suppressWarnings(monoreg(x = x_mat[,joint_ind[i]], 
                y = (y - (as.matrix(x_mat[,joint_ind[-i]]) %*% start_betas[-i]) ), w = weights, type = "isotonic"))
        }
      else if(joint_ind[i] %in% dec_index){
        suppressWarnings(monoreg(x = x_mat[,joint_ind[i]], 
                y = (y - (as.matrix(x_mat[,joint_ind[-i]]) %*% start_betas[-i]) ), w = weights, type = "antitonic"))
      })
    
    # iterate through mr_fits. each column of mr_fits (e.g., mr_fits[,1]) is a monoreg fitted object,
    # and its attributes can be called (e.g., mr_fits[,1]$yf)
    # TODO provide stopping conditions to iteration step
    # TODO remove iters diagnostic
    iters <- 0
    delta <- 0.5
    if(is.null(max_iters_cpav)){
      max_iters_cpav <- 100
    }
    if(is.null(max_delta_cpav)){
      max_delta_cpav <- 0.00001
    }
    
    while(abs(delta) > max_delta_cpav & iters < max_iters_cpav){
      old_SS <- mean( (y - get_pred(mr_fits, x_mat[,joint_ind]))^2 )
      
      for(i in 1:length(joint_ind)){
        if(joint_ind[i] %in% inc_index){
          # i apologize to anyone trying to read this line, but think: the columns of the x_matrix
          # indicated by join_ind, except the value of joint_ind at the ith place in joint_ind
          mr_fits[,i] <- suppressWarnings(monoreg(x = x_mat[,joint_ind[i]],
                                 y = (y - get_pred(mr_fits[,-i], x_mat[,joint_ind[-i]]) ), w = weights, type = "isotonic"))
        }
        else if(joint_ind[i] %in% dec_index){
          mr_fits[,i] <- suppressWarnings(monoreg(x = x_mat[,joint_ind[i]],
                                 y = (y - get_pred(mr_fits[,-i], x_mat[,joint_ind[-i]]) ), w = weights, type = "antitonic"))
        }
        
      }

      new_SS <- mean( (y - get_pred(mr_fits, x_mat[,joint_ind]))^2 )
      
      # TODO get delta in terms of _mean_ SS, so that number of observations doesnt matter
      delta <- (old_SS - new_SS)/old_SS
      
      iters <- iters + 1
    }
    
    return(mr_fits)
  }
}


# first, define function for obtaining f(x_new) for monotone regression f()
# get_pred returns a vector of length = nrows(xvals), ie, a value for each observation of xvals
get_pred <- function(mr_obj, xvals){
  xvals <- as.matrix(xvals)
  mr_obj <- as.matrix(mr_obj)
  
  # TODO this stop is not catching badly formed arguments...
  if(dim(mr_obj)[2] != dim(xvals)[2]) stop("get_pred() must take an X-matrix with as many columns
                                            as monoreg() objects")
  if(dim(xvals)[1] == 1){ # if xvals is a single observation
    sapply(1:ncol(xvals), function(j)
      mr_obj[,j]$yf[sapply(xvals[,j], function(z)
        ifelse( z < mr_obj[,j]$x[1], 1,
                ifelse(z >= tail(mr_obj[,j]$x, n=1), length(mr_obj[,j]$x), 
                       which.min(mr_obj[,j]$x <= z)-1 )))])
  }
  else{
    apply(sapply(1:ncol(xvals), function(j)
         mr_obj[,j]$yf[sapply(xvals[,j], function(z)
           ifelse( z < mr_obj[,j]$x[1], 1,
            ifelse(z >= tail(mr_obj[,j]$x, n=1), length(mr_obj[,j]$x), 
                   which.min(mr_obj[,j]$x <= z)-1 )))]
         ), 1, function(h) sum(h))
    }

}


# define partial linear regression of y on x with weights w
# inputs are: x, y, wates, mon_inc_index, mon_dec_index, max_iter
# TODO max_iter must be passed through M_driver just like mon_inc_index, or 
# else user will not be able to specify the max_iter param
part_fit <- function(x, y, wates = NULL, mon_inc_index=NULL, mon_dec_index=NULL, max_iter=NULL, 
                     component = NULL, na.rm=T, mon_inc_names = NULL, mon_dec_names = NULL, start_fit = NULL, ...){

  
  # TODO cast y and wates to matrices ?
  # TODO correct behaviour for if x is ONLY a vector
  x <- as.matrix(x)
  
  # set default weights
  if(is.null(wates)) wates <- rep(1, length(y))
  
  # remove incomplete cases
  if(T){ # for now, there is no alternative to na.rm=T.  All incomplete cases are removed.
    cc <- complete.cases(y) & complete.cases(x) & complete.cases(wates)
    y <- y[cc]
    x <- x[cc,, drop=FALSE]
    wates <- wates[cc]
    cc <- NULL
  }
  
  x <- as.matrix(x) # cast again. hacky but necessary?
  
  
  # make sure y and wates is not multivariate
  if(length(y) != dim(x)[1] | length(y) != length(wates)) stop("Inputs are not of the same dimension!")
  
  # take monotone indices of previous component
  if(!is.null(component)){
    inc_ind <- component$mon_inc_index
    dec_ind <- component$mon_dec_index
  }
  else{
    
    if(is.null(mon_inc_index) & is.null(mon_dec_index) & is.null(mon_inc_names) & is.null(mon_dec_names)){
      stop("Some monotone index or name must be specified")
    }
    if(!is.null(mon_inc_names)){ # add names to index list
      mon_inc_index <- unique(c(mon_inc_index, which(colnames(x) %in% mon_inc_names)))
    }
    if(!is.null(mon_dec_names)){
      mon_dec_index <- unique(c(mon_dec_index, which(colnames(x) %in% mon_dec_names)))
    }
    # transfer inc_names to inc_index
      inc_ind <- mon_inc_index
      dec_ind <- mon_dec_index
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
  
  # If there is an intercept but no other linear effects, stop
  if((length(c(inc_ind, dec_ind))+1) == ncol(x) & "(Intercept)" %in% colnames(x)){
    stop("For identifiability purposes, you cannot build a part_fit with only an intercept as a linear component.")
  }

  # If start_fit is specified, make sure it has only part_fit elements with the appropriate dimensions
  if(!is.null(start_fit)){
    
    if(!is(start_fit, "list")) stop("start_fit must be a list of part_fit attributes")
    if("coef" %in% names(start_fit)){
      if(length(start_fit$coef) != (dim(x)[2] - length(c(inc_ind, dec_ind))) ){
        stop("Not the right number of coefficients in starting values")
      }
      if(!all(names(start_fit$coef) %in% colnames(x)[-c(inc_ind, dec_ind)]) ){
        stop("Some coefficient(s) in starting values have incorrect names")
      }
      if(!all(colnames(x)[-c(inc_ind, dec_ind)] %in% names(start_fit$coef)) ){
        stop("Some coefficient(s) in starting values are missing")
      }
      if(!all(names(start_fit$coef) == colnames(x)[-c(inc_ind, dec_ind)])){
        start_fit$coef <- start_fit$coef[match(colnames(x)[-c(inc_ind, dec_ind)],start_fit$coef)]
      }
    }
    # TODO use starting mon_obj
    # if(mon_obj %in% names(start_fit)){
    #   
    # }
  }
  
  
  # TODO get names of monotone components for plotting
  # if(is.null(c(mon_inc_names, mon_dec_names))){
  #   # set names according to monotone indices. First increasing indices, then decreasing indices
  # }
  
  # option for fit with no linear independent components and one or multiple monotone components:
  if(length(c(inc_ind, dec_ind)) == ncol(x)){
    
    yhat <- cpav(x_mat = as.matrix(x[wates != 0,]), y = y[wates != 0], weights = wates[wates != 0], 
                 inc_index = inc_ind, dec_index = dec_ind)

    # get residuals of model
    resids <- y - get_pred(yhat, x[,c(inc_ind, dec_ind)])
    
    # mod must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
    # may have mon_inc_index and mon_dec_index attributes
    mod <- list(coef = NULL, fitted_pava = NULL, sigma = NULL, df = NULL,
                mon_inc_index = NULL, mon_dec_index = NULL, iterations = NULL, 
                mon_inc_names = NULL, mon_dec_names = NULL)
  
    mod$coef <- NULL
    mod$fitted_pava <- yhat
    mod$mon_inc_index <- inc_ind
    mod$mon_dec_index <- dec_ind
    # TODO sigma is the sqrt of __ divided by (nrow(x) - model$rank). For single monoreg, rank is 1...
    # ... but update when implementing CPAV
    mod$sigma <- sqrt(sum(wates * (resids)^2 /
                                     mean(wates))/ (nrow(x)-qr(x)$rank))
    mod$df <- ncol(x)+1
    
    class(mod) <- "part_fit"
    
    return(mod)
  }
  else{
    # for starting values, fit a regular lm or use starting values
    if(!is.null(start_fit) && "coef" %in% names(start_fit)){
      betas <- start_fit$coef
    }
    else{
      fit <- lm.wfit(x=x, y=y, w=wates)
      betas <- coef(fit)[-c(inc_ind, dec_ind)]
    }
  
    # set maximum iterations for convergence
    if(!is.null(max_iter) & !is.list(max_iter)){
      if(max_iter < 1) stop("max_iter must be positive")
      maxiter <- max_iter
    }
    else{ 
      maxiter <- 200
    }
      
    # set while loop initial values
    iter <- 0
    delta <- 10
    # iterate between pava and linear model
    # TODO set while loop condition(s). Get appropriate measure of coefficient change
    while(delta > 1e-6 & iter < maxiter){ # works well enough with delta > 1e-12. Trying 1e-6
  
      yhat <- cpav(x_mat = as.matrix(x[wates != 0,]), y = (y[wates != 0] - (as.matrix(x[wates != 0,-c(inc_ind, dec_ind)]) %*% betas)), 
                   weights = wates[wates != 0], inc_index = inc_ind, dec_index = dec_ind)
      # TODO remove below call to monoreg()
      # monoreg(x = x[wates != 0,inc_ind], 
      #                 y = (y[wates != 0] - x[wates != 0,-inc_ind] %*% betas), w = wates[wates != 0])
      
      old_betas <- betas    # save old betas for distance calculation
      # to retrieve old ordering of y for fitted values, we use y[match(x, sorted_x)]
      betas <- coef(lm.wfit(x=as.matrix(x[,-c(inc_ind, dec_ind)]), y= (y - get_pred(yhat, x[,c(inc_ind, dec_ind)]) ), w=wates))

      # TODO quantify change in yhat vals and beta vals
      # get euclidian distance between betas transformed into unit vectors
      delta <- dist(rbind( as.vector(betas)/norm(as.vector(betas), type ="2"), 
                           as.vector(old_betas)/norm(as.vector(old_betas), type="2")
                           ))
  
      iter <- iter + 1    # iterate maxiter
    }
  }
  
  # get residuals of model
  resids <- y - (get_pred(yhat, x[,c(inc_ind, dec_ind)]) + (as.matrix(x[,-c(inc_ind, dec_ind)]) %*% betas))

  # mod must have: coef attribute, sigma attribute, cov attribute, df attribute, ..., and 
  # may have mon_inc_index and mon_dec_index attributes
  mod <- list(coef = NULL, fitted_pava = NULL, sigma = NULL, df = NULL,
              mon_inc_index = NULL, mon_dec_index = NULL, iterations = NULL, 
              mon_inc_names = NULL, mon_dec_names = NULL)
  
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

  class(mod) <- "part_fit"
  
  return(mod)
}


# write plot method for objects returned from part_fit()

append_suffix <- function(num){
  suff <- case_when(num %in% c(11,12,13) ~ "th",
                    num %% 10 == 1 ~ 'st',
                    num %% 10 == 2 ~ 'nd',
                    num %% 10 == 3 ~'rd',
                    TRUE ~ "th")
  paste0(num, suff)
}

plot.part_fit <- function(z){
  if(dim(as.matrix(z$fitted_pava))[2] > 1){
    temp <- list()
    for(i in 1:dim(as.matrix(z$fitted_pava))[2]){
      temp[[i]] <- ggplotGrob(ggplot() +
        geom_line(aes(x = z$fitted_pava[,i]$x, y = z$fitted_pava[,i]$yf)) + 
        theme_bw() +
        labs(title = paste(append_suffix(i), " Monotone Regression"),
             x = "X",
             y = "Y"))
    }
    return(grid.arrange(grobs=temp, ncol=1))
  }
  else{
    temp <- ggplot() +
      geom_line(aes(x = z$fitted_pava[,1]$x, y = z$fitted_pava[,1]$yf)) + 
      theme_bw() +
      labs(title = "Monotone Regression",
           x = "X",
           y = "Y")
    return(temp)
  }
}
