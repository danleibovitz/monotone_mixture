## Wishlist

# - the override plot method should be able to include the rootogram from the old plot method, e.g. 
# by calling flexmix::plot() within the plot method. This, however, seems to call the current method --
# even though the current method is not within the flexmix package -- initiating infinite recursion
# - predict method

# wrap various flexmix methods

library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(zoo)
library(lemon)
library(knitr)
library(LaplacesDemon)
library(LICORS)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


####


## TODO mon_obj is in *different* location now

# overwrite method for plot.flexmix
setMethod('plot',  signature(x="flexmix", y="missing"),
          function(x, mark=NULL, markcol=NULL, col=NULL, 
                   eps=1e-4, root=TRUE, ylim=NULL, xlim=NULL, main=NULL, xlab=NULL, ylab=NULL,
                   as.table = TRUE, endpoints = c(-0.04, 1.04), rootogram=F, palet = NULL, 
                   root_scale = "unscaled", subplot=NULL, boot_CI = NULL, ...) {
            
            if(is.null(palet)){
              palet <- "Accent"
            }
            
          
  if(is(x@model[[1]], "FLXM_monoreg")){ # check that this is a mixture of part_fits
    # assign appropriate names for graph labelling
    if(is.null( c(x@model[[1]]@mon_inc_names, x@model[[1]]@mon_dec_names) )){
      xnames <- sapply(1:dim(x@components$Comp.1[[1]]@parameters$mon_obj)[2], function(x) paste0("X", x))
      mono_names <- c("Y", xnames)
    }
    else{
      mono_names <- c(x@formula[[2]], c(x@model[[1]]@mon_inc_names, x@model[[1]]@mon_dec_names))
    }
    
    if(!is.null(boot_CI)){ # check that bootstrapped CIs are of the correct dimesion
      if(length(boot_CI$coef_CI) != length(x@components) | length(boot_CI$mono_CI) != length(x@components)){
        stop("The boot_CI object has the wrong number of components")}
      if(dim(boot_CI$coef_CI[[1]])[1] != length(x@components[[1]][[1]]@parameters$coef) + 1 | 
         length(boot_CI$mono_CI[[1]]) != length(c(x@model[[1]]@mon_inc_names, x@model[[1]]@mon_dec_names))){
        stop("The boot_CI object has the wrong number of monotone effects or of linear effects")
      }
    }
            if(dim( x@components$Comp.1[[1]]@parameters$mon_obj )[2] > 1){ # get dimension of monotone components by reading columns of fitted_pava object
              np <- list()
              for(i in 1:dim( x@components$Comp.1[[1]]@parameters$mon_obj )[2]){
                holder <- ggplot()
                
                if(length(x@components) == 1){
                    if(is.null(boot_CI)){
                      holder <- holder + 
                      geom_step(aes(x = x@components$Comp.1[[1]]@parameters$mon_obj[,i]$x,
                                y = x@components$Comp.1[[1]]@parameters$mon_obj[,i]$yf)) +
                      theme_bw() +
                      labs(title = paste(append_suffix(i), " Monotone Regression"),
                       x = mono_names[i+1],
                       y = mono_names[1])}
                    else{ # if there is a boot_CI object, add the CI intervals
                      holder <- holder + geom_step(aes(x = boot_CI$mono_CI[[1]][[i]][,1],
                                                       y = boot_CI$mono_CI[[1]][[i]][,4])) +
                        theme_bw() +
                        labs(title = paste(append_suffix(i), " Monotone Regression"),
                             x = mono_names[i+1],
                             y = mono_names[1]) +
                        geom_stepconfint(aes(x = boot_CI$mono_CI[[1]][[i]][,1],
                                                              ymin = boot_CI$mono_CI[[1]][[i]][,2], 
                                                              ymax = boot_CI$mono_CI[[1]][[i]][,3]), alpha = 0.3)
                    }
                    }
                
                if(length(x@components) > 1){
                  
                  monlist <- list()
                  for(b in 1:length(x@components)){
                    if(!is.null(boot_CI)){
                      monx <- sort(union(x@components[[b]][[1]]@parameters$mon_obj[,i]$x, 
                                         boot_CI$mono_CI[[b]][[i]][,1]))
                      monyf <- rep(NA, length(monx))
                      monyf[match(x@components[[b]][[1]]@parameters$mon_obj[,i]$x, monx)] <- 
                        x@components[[b]][[1]]@parameters$mon_obj[,i]$yf
                      monyf <- na.locf0(na.locf0(monyf), fromLast = TRUE)
                      
                      monlist[[b]] <- data.frame(x = monx, 
                                                 yf = monyf,
                                                 upper = boot_CI$mono_CI[[b]][[i]][,2],
                                                 lower = boot_CI$mono_CI[[b]][[i]][,3],
                                                 mean = boot_CI$mono_CI[[b]][[i]][,4])
                    }
                    else{
                      monlist[[b]] <- data.frame(x = x@components[[b]][[1]]@parameters$mon_obj[,i]$x, 
                                               yf = x@components[[b]][[1]]@parameters$mon_obj[,i]$yf)
                    }
                  }
                  
                  mondf   <- cbind(Cluster=rep(1:length(x@components),sapply(monlist,nrow)),do.call(rbind,monlist))
                  mondf$Cluster <- as.factor(mondf$Cluster)
                  
                  if(is.null(boot_CI)){
                    holder <- holder + geom_step(mondf, mapping = aes(x,yf, color=Cluster)) + 
                    scale_color_brewer(palette=palet) +
                    theme_bw() +
                    labs(title = paste(append_suffix(i), " Monotone Regression"),
                         x = mono_names[i+1],
                         y = mono_names[1])}
                  else{
                    holder <- holder + geom_step(mondf, mapping = aes(x, mean, color=Cluster)) + 
                      theme_bw() +
                      labs(title = paste(append_suffix(i), " Monotone Regression"),
                           x = mono_names[i+1],
                           y = mono_names[1]) +
                      geom_stepconfint(mondf, mapping = aes(x = x, ymin = lower, 
                                                            ymax = upper, fill=Cluster), alpha = 0.3) +
                      scale_color_brewer(aesthetics = c("colour", "fill"), palette=palet)
                  }
                }
                
                if(!is.null(ylim)){
                  if(length(ylim) != dim( x@components[[1]][[1]]@parameters$mon_obj )[2] |
                     length(ylim[[1]]) != 2 ){
                    stop("If you pass a ylim argument, it must have as many element pairs 
                         as the model has monotone components. Try formulating the argument
                         as: ylim = list(c(i,j), c(i,j), ...)")}
                  holder <- holder + ylim(ylim[[i]])
                }
                if(!is.null(xlim)){
                  if(length(xlim) != dim( x@components[[1]][[1]]@parameters$mon_obj )[2] |
                     length(xlim[[1]]) != 2 ){
                    stop("If you pass a xlim argument, it must have as many element pairs 
                         as the model has monotone components. Try formulating the argument
                         as: xlim = list(c(i,j), c(i,j), ...)")}
                  holder <- holder + xlim(xlim[[i]])
                }
                if(!is.null(ylab)){
                  if(length(ylab) != dim( x@components[[1]][[1]]@parameters$mon_obj )[2]){
                    stop("If you pass a ylab argument, it must have as many elements 
                         as the model has monotone components. Try formulating the argument
                         as: ylab = c(\"first\",\"second\",...)")}
                  holder <- holder + ylab(ylab[[i]])
                }
                if(!is.null(xlab)){
                  if(length(xlab) != dim( x@components[[1]][[1]]@parameters$mon_obj )[2]){
                    stop("If you pass a xlab argument, it must have as many elements 
                         as the model has monotone components. Try formulating the argument
                         as: xlab = c(\"first\",\"second\",...)")}
                  holder <- holder + xlab(xlab[[i]])
                }
                if(!is.null(main)){
                  if(length(main) != dim( x@components[[1]][[1]]@parameters$mon_obj )[2]){
                    stop("If you pass a main argument, it must have as many elements 
                         as the model has monotone components. Try formulating the argument
                         as: main = c(\"first\",\"second\",...)")}
                  holder <- holder + ggtitle(main[[i]])
                }
                
                
                
                np[[i]] <- holder
              }
              # return(grid.arrange(grobs=np, ncol=1))
              # return(grid.arrange(grobs=np, ncol=1))
            }
            else{
              np <- ggplot()
              
              if(length(x@components) == 1){
                if(is.null(boot_CI)){
                  np <- np + geom_step(aes(x = x@components[[1]][[1]]@parameters$mon_obj[,1]$x, 
                                         y = x@components[[1]][[1]]@parameters$mon_obj[,1]$yf)) + 
                  theme_bw() +
                  labs(title = "Monotone Component",
                     x = mono_names[2],
                     y = mono_names[1])}
                else{
                  np <- np + geom_step(aes(x = boot_CI$mono_CI[[1]][[1]][,1], 
                                           y = boot_CI$mono_CI[[1]][[1]][,4])) + 
                    theme_bw() +
                    labs(title = "Monotone Component",
                         x = mono_names[2],
                         y = mono_names[1]) +
                    geom_stepconfint(aes(x = boot_CI$mono_CI[[1]][[1]][,1],
                                                  ymin = boot_CI$mono_CI[[1]][[1]][,2], 
                                                  ymax = boot_CI$mono_CI[[1]][[1]][,3]), alpha = 0.3)
                }
                }
              
              if(length(x@components) > 1){
                
                monlist <- list()
                for(b in 1:length(x@components)){
                  if(!is.null(boot_CI)){
                    monx <- sort(union(x@components[[b]][[1]]@parameters$mon_obj[,1]$x, 
                                       boot_CI$mono_CI[[b]][[1]][,1]))
                    monyf <- rep(NA, length(monx))
                    monyf[match(x@components[[b]][[1]]@parameters$mon_obj[,1]$x, monx)] <- 
                      x@components[[b]][[1]]@parameters$mon_obj[,1]$yf
                    monyf <- na.locf0(na.locf0(monyf), fromLast = TRUE)
                    
                    monlist[[b]] <- data.frame(x = monx, 
                                               yf = monyf,
                                               upper = boot_CI$mono_CI[[b]][[1]][,2],
                                               lower = boot_CI$mono_CI[[b]][[1]][,3],
                                               mean = boot_CI$mono_CI[[b]][[1]][,4])
                  }
                  else{
                    monlist[[b]] <- data.frame(x = x@components[[b]][[1]]@parameters$mon_obj[,1]$x, 
                                                   yf = x@components[[b]][[1]]@parameters$mon_obj[,1]$yf)
                  }
                }
                
                mondf   <- cbind(Cluster=rep(1:length(x@components),sapply(monlist,nrow)),do.call(rbind,monlist))
                mondf$Cluster <- as.factor(mondf$Cluster)
                
                if(is.null(boot_CI)){
                  np <- np + geom_step(mondf, mapping = aes(x,yf, color=Cluster)) + 
                  scale_color_brewer(palette=palet) +
                  theme_bw() +
                  labs(title = "Monotone Component",
                       x = mono_names[2],
                       y = mono_names[1])}
                else{
                  np <- np + geom_step(mondf, mapping = aes(x, mean, color=Cluster)) + 
                    theme_bw() +
                    labs(title = "Monotone Component",
                         x = mono_names[2],
                         y = mono_names[1]) +
                    geom_stepconfint(mondf, mapping = aes(x = x, ymin = lower, 
                                                          ymax = upper, fill=Cluster), alpha = 0.3) + 
                    scale_color_brewer(aesthetics = c("colour", "fill"), palette=palet)
                }
              }
              
              
              if(!is.null(ylim)){
                np <- np + ylim(ylim)
              }
              if(!is.null(xlim)){
                np <- np + xlim(xlim)
              }
              if(!is.null(ylab)){
                np <- np + ylab(ylab)
              }
              if(!is.null(xlab)){
                np <- np + xlab(xlab)
              }
              if(!is.null(main)){
                np <- np + ggtitle(main)
              }
              
              
              
              
              # return(np)
            }
  }
            else{ stop("This method is only for FLXM_monoreg objects")}
            # plot and append rootogram
            post <- data.frame(x@posterior$scaled) # collect posteriors
            names(post) <- 1:dim(post)[2] # change columns of posteriors to cluster numbers
            post <- melt(setDT(post), measure.vars = c(1:dim(post)[2]), variable.name = "Cluster")
            rg <- ggplot(post, aes(x=value, fill=Cluster)) + # plot rootogram, with color indicating cluster
              geom_histogram(binwidth = 0.05) + 
              scale_fill_brewer(palette=palet) +
              theme_bw() +
              labs(title = "Rootogram",
                   x = "Posteriors",
                   y = "Count")
            
            if(root_scale == "sqrt"){rg <- rg + 
              scale_y_sqrt() +
              labs(title = "Rootogram (square root scale)",
                   x = "Posteriors",
                   y = "Count (square root)")}
            if(root_scale == "log"){rg <- rg + 
              scale_y_log10() + 
              labs(title = "Rootogram (log scale)",
                   x = "Posteriors",
                   y = "Count (log)")}
            
            if(!is.null(subplot)){
              return(list(rg, np)[[subplot[1]]])
            }
            else{
              multiplot(rg, np)
            }
          }          
)
            
# TODO include in plot method for monoreg flexmix objects a forest plot of coefficients and a forest plot of priors.
# These plots should generate with or without bootstrap; without bootstrap, they should be point estimates (unless we can find a tractable parametric estimate of the coefficients without having to calculate some ridiculous covariance matrix with a kernel method)





# TODO write predict method for monoreg flexmix objects
predict_marginal <- function(fm, dat){
  # this function is not written defensively. It expects properly formatted inputs. You should probably not use it. I apologize for laziness...
  # return a density for each observation. should density be a matrix or function? if a matrix, must specify the granularity of the range
  
  # if(dim(dat)[1] > 1 && observations dont belong to same grouping){
  #   # class(mod6@group)
  #   stop("run predict_marginal for one grouping at a time.")
  # }
  
  if(all.vars(fm@formula)[1] %in% colnames(dat)){
    # return ONE posterior vector for entire group
    retval <- apply()
  }
  
  if(!all.vars(fm@formula)[1] %in% colnames(dat)){    
    # return a marginal distribution FOR EACH observation in the form of a vector of priors, 
    # a vector of mus, and a vector of sigmas
    sigs <- sapply(fm@components, function(i) i[[1]]@parameters$sigma)
    retval <- lapply(1:dim(dat)[1], function(x) {
      list(prior = normalize(fm@prior), mus = sapply(fm@components, function(j) j[[1]]@predict(dat[x,])), # get predicted value at row x for each component 
           sigmas = sigs)
    })
  }
  
  return(retval)
}

            
            


### building a confidence interval for an ordinary bootstrapped, component monotone function
build_CI_trad <- function(fboot, CI = 0.89, HDI = FALSE){
  
  if(!is(fboot, "FLXboot") | !is(fboot@object@model[[1]], "FLXM_monoreg")) stop("Can only run build_CI_trad for FLXboot objects containing FLXM_monoreg objects") # only run for FLXboot objects
  if(CI <= 0 | CI >= 1) stop("CI is incorrect magnitude")
  mon_names <- c(fboot@object@model[[1]]@mon_inc_names, fboot@object@model[[1]]@mon_dec_names)
  R <- length(fboot@parameters)
  
  # Build empirical distribution of coefficients
  cdim <- length(fboot@object@components[[1]][[1]]@parameters$coef)
  coef_list <- rep(list(matrix(nrow = cdim+1, ncol = R+1, 
                               dimnames = list(c(names(fboot@object@components[[1]][[1]]@parameters$coef), "sigma"), 1:(R+1)))), 
                   length(fboot@object@components) ) # make 3-D list: components, monotone indices, bootstrap repitions
  prior_CI <- matrix(nrow = length(fboot@object@components), ncol = R+1, 
                     dimnames = list(1:length(fboot@object@components), 1:(R+1)))
    
    
  for(k in 1:length(fboot@object@components)){ # populate coef_list
    if(cdim > 0){
      coef_list[[k]][1:cdim,1] <- fboot@object@components[[k]][[1]]@parameters$coef}
    coef_list[[k]][cdim+1,1] <- fboot@object@components[[k]][[1]]@parameters$sigma
    prior_CI[k, 1] <- fboot@object@prior[k]
    for(r in 1:R){
      if(cdim > 0){
        coef_list[[k]][1:cdim,r+1] <- fboot@parameters[[r]][[1]][[1]][[k]]$coef}
      coef_list[[k]][cdim+1,r+1] <- fboot@parameters[[r]][[1]][[1]][[k]]$sigma
      prior_CI[k, r+1] <- fboot@priors[[r]][[1]][k]
    }
  }
  
  if(HDI == TRUE){ # calculate HDI intervals
    for(k in 1:length(fboot@object@components)){ # populate coef_list
      coef_list[[k]] <- cbind(t(hdi(t(coef_list[[k]]), CI)), rowMeans(coef_list[[k]]))
      colnames(coef_list[[k]])[3] <- "mean"
    }
    prior_CI <- cbind(t(hdi(t(prior_CI), CI)), rowMeans(prior_CI))
    colnames(prior_CI)[3] <- "mean"
    
  }
  else{ # calculate confidence intervals
    for(k in 1:length(fboot@object@components)){ # populate coef_list
      coef_list[[k]] <- cbind(t(apply(coef_list[[k]], 1, function(j) quantile(j, c((1-CI)/2, 1 - (1-CI)/2 ) ))),
                              rowMeans(coef_list[[k]]))
      colnames(coef_list[[k]]) <- c("lower", "upper", "mean")
    }
    prior_CI <- cbind(t(apply(prior_CI, 1, function(j) quantile(j, c((1-CI)/2, 1 - (1-CI)/2 ) ))),
                            rowMeans(prior_CI))
    colnames(prior_CI) <- c("lower", "upper", "mean")
  }
  
  # Build empirical distribution of monotone fits
  CI_list <- rep(list(  rep( list(NA),length(mon_names))   ), 
                 length(fboot@object@components) ) # make 3-D list: components, monotone indices, bootstrap repitions
  
  for(k in 1:length(fboot@object@components)){ # populate coef_list
    for(j in 1:dim(fboot@object@components[[k]][[1]]@parameters$mon_obj)[2]){
      xrange <- fboot@object@components[[k]][[1]]@parameters$mon_obj[,j]$x
      for(r in 1:R){ # get union of all x vals for a component
        xrange <- union(xrange, fboot@parameters[[r]][[1]][[1]][[k]]$mon_obj[,j]$x)
      }
      xrange <- sort(xrange)
      CI_list[[k]][[j]] <- matrix(nrow = length(xrange), ncol = R+2)
      CI_list[[k]][[j]][,1] <- xrange
      CI_list[[k]][[j]][match(fboot@object@components[[k]][[1]]@parameters$mon_obj[,j]$x, xrange),2] <- 
        fboot@object@components[[k]][[1]]@parameters$mon_obj[,j]$yf
      for(r in 1:R){
        CI_list[[k]][[j]][match(fboot@parameters[[r]][[1]][[1]][[k]]$mon_obj[,j]$x, xrange),r+2] <- 
          fboot@parameters[[r]][[1]][[1]][[k]]$mon_obj[,j]$yf
      }
      # fill out matrix
      CI_list[[k]][[j]] <- apply(CI_list[[k]][[j]], 2, function(col)
        na.locf0(na.locf0(col), fromLast = TRUE))
      
      if(HDI == TRUE){ # calculate HDI intervals
        CI_list[[k]][[j]][,2:4] <- c(t(hdi(t(CI_list[[k]][[j]][,2:(R+1)]), CI)), rowMeans(CI_list[[k]][[j]][,2:(R+1)]))
        CI_list[[k]][[j]] <- CI_list[[k]][[j]][,1:4]
        colnames(CI_list[[k]][[j]]) <- c("x", "lower", "upper", "mean")
      }
      else{ # calculate confidence intervals
        CI_list[[k]][[j]][,2:4] <- c(t(apply(CI_list[[k]][[j]][,2:(R+1)], 1, 
                                             function(row) quantile(row, c((1-CI)/2, 1 - (1-CI)/2 ) ))),
                                     rowMeans(CI_list[[k]][[j]][,2:(R+1)]))
        CI_list[[k]][[j]] <- CI_list[[k]][[j]][,1:4]
        colnames(CI_list[[k]][[j]]) <- c("x", "lower", "upper", "mean")
      }
    }
  }
  

  return(list(coef_CI = coef_list, mono_CI = CI_list, prior_CI = prior_CI, R=R, CI=CI, HDI=HDI))
  
  
}

### building a confidence interval for a conditional bootstrapped, component monotone function
build_CI_cond <- function(fmix, CI = 0.89, HDI = FALSE, R = 1000){
  
  if(!is(fmix@model[[1]], "FLXM_monoreg")) stop("Can only run build_CI_cond for FLXM_monoreg objects") 
  if(CI <= 0 | CI >= 1) stop("CI is incorrect magnitude")
  mon_names <- c(fmix@model[[1]]@mon_inc_names, fmix@model[[1]]@mon_dec_names)
  PF_list <- rep(list(  rep( list(NA),R)   ), 
                 length(fmix@components) ) # make 3-D list of part_fits: components, monotone indices, bootstrap repitions
  
  n <- length(fmix@model[[1]]@y)
  
  for(k in 1:length(fmix@components)){ # populate CI_list with a part_fit for each bootstrap iteration
    PF_list[[k]][[1]] <- list(coef = fmix@components[[k]][[1]]@parameters$coef,
                              sigma = fmix@components[[k]][[1]]@parameters$sigma,
                              fitted_pava = fmix@components[[k]][[1]]@parameters$mon_obj)
    for(r in 2:R){
      # populate via conditional bootstrap
      # TODO provide starting values
      sam <- sample(1:n, n, replace = TRUE, 
                    prob = fmix@posterior$scaled[,k])
      PF_list[[k]][[r]] <- part_fit(x = fmix@model[[1]]@x[sam,,drop=FALSE], y = fmix@model[[1]]@y[sam], 
                                    wates = rep.int(1,n), 
                                    mon_inc_index = fmix@model[[1]]@mon_inc_index, 
                                    mon_dec_index = fmix@model[[1]]@mon_dec_index, 
                                    mon_inc_names = fmix@model[[1]]@mon_inc_names, 
                                    mon_dec_names = fmix@model[[1]]@mon_dec_names,
                                    start_fit = list(coef = fmix@components[[k]][[1]]@parameters$coef,
                                                     mon_obj = fmix@components[[k]][[1]]@parameters$mon_obj)
      )
    }
  }
  
  # Build empirical distribution of coefficients
  cdim <- length(fmix@components[[1]][[1]]@parameters$coef)
  coef_list <- rep(list(matrix(nrow = cdim+1, ncol = R, 
                               dimnames = list(c(names(fmix@components[[1]][[1]]@parameters$coef), "sigma"), 1:R))), 
                   length(fmix@components) ) # make 3-D list: components, monotone indices, bootstrap repitions
  
  for(k in 1:length(fmix@components)){ # populate coef_list
    for(r in 1:R){
      if(cdim > 0){
        coef_list[[k]][1:cdim,r] <- PF_list[[k]][[r]]$coef
      }
      coef_list[[k]][cdim+1,r] <- PF_list[[k]][[r]]$sigma
    }
  }
  
  if(HDI == TRUE){ # calculate HDI intervals
    for(k in 1:length(fmix@components)){ # populate coef_list
      coef_list[[k]] <- cbind(t(hdi(t(coef_list[[k]]), CI)), rowMeans(coef_list[[k]]))
      colnames(coef_list[[k]])[3] <- "mean"
    }
    
  }
  else{ # calculate confidence intervals
    for(k in 1:length(fmix@components)){ # populate coef_list
      coef_list[[k]] <- cbind(t(apply(coef_list[[k]], 1, function(j) quantile(j, c((1-CI)/2, 1 - (1-CI)/2 ) ))),
                              rowMeans(coef_list[[k]]))
      colnames(coef_list[[k]]) <- c("lower", "upper", "mean")
    }
  }
  
  
  # Build empirical distribution of monotone fits
  CI_list <- rep(list(  rep( list(NA),length(mon_names))   ), 
                 length(fmix@components) ) # make 3-D list: components, monotone indices, bootstrap repitions
  
  for(k in 1:length(fmix@components)){ # populate coef_list
    for(j in 1:dim(PF_list[[k]][[1]]$fitted_pava)[2]){
      xrange <- PF_list[[k]][[1]]$fitted_pava[,j]$x
      for(r in 2:R){ # get union of all x vals for a component
        xrange <- union(xrange, PF_list[[k]][[r]]$fitted_pava[,j]$x)
      }
      xrange <- sort(xrange)
      CI_list[[k]][[j]] <- matrix(nrow = length(xrange), ncol = R+1)
      CI_list[[k]][[j]][,1] <- xrange
      for(r in 1:R){
        CI_list[[k]][[j]][match(PF_list[[k]][[r]]$fitted_pava[,j]$x, xrange),r+1] <- 
          PF_list[[k]][[r]]$fitted_pava[,j]$yf
      }
      # fill out matrix
      CI_list[[k]][[j]] <- apply(CI_list[[k]][[j]], 2, function(col)
        na.locf0(na.locf0(col), fromLast = TRUE))
      
      if(HDI == TRUE){ # calculate HDI intervals
        CI_list[[k]][[j]][,2:4] <- c(t(hdi(t(CI_list[[k]][[j]][,2:R]), CI)), rowMeans(CI_list[[k]][[j]][,2:R]))
        CI_list[[k]][[j]] <- CI_list[[k]][[j]][,1:4]
        colnames(CI_list[[k]][[j]]) <- c("x", "lower", "upper", "mean")
      }
      else{ # calculate confidence intervals
        CI_list[[k]][[j]][,2:4] <- c(t(apply(CI_list[[k]][[j]][,2:R], 1, 
                                             function(row) quantile(row, c((1-CI)/2, 1 - (1-CI)/2 ) ))),
                                     rowMeans(CI_list[[k]][[j]][,2:R]))
        CI_list[[k]][[j]] <- CI_list[[k]][[j]][,1:4]
        colnames(CI_list[[k]][[j]]) <- c("x", "lower", "upper", "mean")
      }
    }
  }
  
  
  return(list(coef_CI = coef_list, mono_CI = CI_list, R=R, CI=CI, HDI=HDI))
}





### building a confidence interval table for an ordinary bootstrapped, component monotone function
CI_table <- function(CImat){
  
  knit_print.data.frame <- lemon_print
  
  # CImat attributes: CImat$coef_list, CImat$CI_list, CImat$R, CImat$CI, CImat$HDI
  # df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("lower", "upper", "mean", "cluster"))))
  df <- list()
  
  for(i in 1:length(CImat$coef_CI)){
    # make dataframe
    # d <- cbind(CImat$coef_CI[[i]], rep.int(i, dim(CImat$coef_CI[[i]])[1]))
    # df[(dim(df)[1]+1):(dim(df)[1]+dim(CImat$coef_CI[[i]])[1]), ] <- as.data.frame(d)
    df[[i]] <- as.data.frame(CImat$coef_CI[[i]])
  }
  
  
  # make priors plot
  df[[i+1]] <- CImat$prior_CI
  return(df)
}
            
            