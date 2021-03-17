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

# overwrite method for plot.flexmix
setMethod('plot',  signature(x="flexmix", y="missing"),
          function(x, mark=NULL, markcol=NULL, col=NULL, 
                   eps=1e-4, root=TRUE, ylim=NULL, xlim=NULL, main=NULL, xlab=NULL, ylab=NULL,
                   as.table = TRUE, endpoints = c(-0.04, 1.04), rootogram=F, palet = NULL, 
                   root_scale = "unscaled", subplot=NULL, ...) {
            
            if(is.null(palet)){
              palet <- "Accent"
            }
            
           
  if(is(x@components[[1]][[1]], "FLX_monoreg_component")){ # check that this is a mixture of part_fits
    # assign appropriate names for graph labelling
    if(is.null( c(x@model[[1]]@mon_inc_names, x@model[[1]]@mon_dec_names) )){
      xnames <- sapply(1:dim(x@components[[1]][[1]]@mon_obj)[2], function(x) paste0("X", x))
      mono_names <- c("Y", xnames)
    }
    else{
      mono_names <- c(x@formula[[2]], c(x@model[[1]]@mon_inc_names, x@model[[1]]@mon_dec_names))
    }
    
            if(dim( x@components[[1]][[1]]@mon_obj )[2] > 1){ # get dimension of monotone components by reading columns of fitted_pava object
              np <- list()
              for(i in 1:dim( x@components[[1]][[1]]@mon_obj )[2]){
                holder <- ggplot()
                
                if(length(x@components) == 1){
                    holder <- holder + 
                      geom_line(aes(x = x@components[[1]][[1]]@mon_obj[,i]$x,
                                y = x@components[[1]][[1]]@mon_obj[,i]$yf)) +
                      theme_bw() +
                      labs(title = paste(append_suffix(i), " Monotone Regression"),
                       x = mono_names[i+1],
                       y = mono_names[1])
                    }
                
                if(length(x@components) > 1){
                  
                  monlist <- list()
                  for(b in 1:length(x@components)){
                    monlist[[b]] <- data.frame(x = x@components[[b]][[1]]@mon_obj[,i]$x, 
                                               yf = x@components[[b]][[1]]@mon_obj[,i]$yf)
                  }
                  
                  mondf   <- cbind(Cluster=rep(1:length(x@components),sapply(monlist,nrow)),do.call(rbind,monlist))
                  mondf$Cluster <- as.factor(mondf$Cluster)
                  
                  holder <- holder + geom_line(mondf, mapping = aes(x,yf, color=Cluster)) + 
                    scale_color_brewer(palette=palet) +
                    theme_bw() +
                    labs(title = paste(append_suffix(i), " Monotone Regression"),
                         x = mono_names[i+1],
                         y = mono_names[1])
                }
                
                if(!is.null(ylim)){
                  if(length(ylim) != dim( x@components[[1]][[1]]@mon_obj )[2] |
                     length(ylim[[1]]) != 2 ){
                    stop("If you pass a ylim argument, it must have as many element pairs 
                         as the model has monotone components. Try formulating the argument
                         as: ylim = list(c(i,j), c(i,j), ...)")}
                  holder <- holder + ylim(ylim[[i]])
                }
                if(!is.null(xlim)){
                  if(length(xlim) != dim( x@components[[1]][[1]]@mon_obj )[2] |
                     length(xlim[[1]]) != 2 ){
                    stop("If you pass a xlim argument, it must have as many element pairs 
                         as the model has monotone components. Try formulating the argument
                         as: xlim = list(c(i,j), c(i,j), ...)")}
                  holder <- holder + xlim(xlim[[i]])
                }
                if(!is.null(ylab)){
                  if(length(ylab) != dim( x@components[[1]][[1]]@mon_obj )[2]){
                    stop("If you pass a ylab argument, it must have as many elements 
                         as the model has monotone components. Try formulating the argument
                         as: ylab = c(\"first\",\"second\",...)")}
                  holder <- holder + ylab(ylab[[i]])
                }
                if(!is.null(xlab)){
                  if(length(xlab) != dim( x@components[[1]][[1]]@mon_obj )[2]){
                    stop("If you pass a xlab argument, it must have as many elements 
                         as the model has monotone components. Try formulating the argument
                         as: xlab = c(\"first\",\"second\",...)")}
                  holder <- holder + xlab(xlab[[i]])
                }
                if(!is.null(main)){
                  if(length(main) != dim( x@components[[1]][[1]]@mon_obj )[2]){
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
                np <- np + geom_line(aes(x = x@components[[1]][[1]]@mon_obj[,1]$x, y = x@components[[1]][[1]]@mon_obj[,1]$yf)) + 
                  theme_bw() +
                  labs(title = "Monotone Component",
                     x = mono_names[2],
                     y = mono_names[1]) 
                }
              
              if(length(x@components) > 1){
                
                monlist <- list()
                for(b in 1:length(x@components)){
                  monlist[[b]] <- data.frame(x = x@components[[b]][[1]]@mon_obj[,1]$x, 
                                                   yf = x@components[[b]][[1]]@mon_obj[,1]$yf)
                }
                
                mondf   <- cbind(Cluster=rep(1:length(x@components),sapply(monlist,nrow)),do.call(rbind,monlist))
                mondf$Cluster <- as.factor(mondf$Cluster)
                
                np <- np + geom_line(mondf, mapping = aes(x,yf, color=Cluster)) + 
                  scale_color_brewer(palette=palet) +
                  theme_bw() +
                  labs(title = "Monotone Component",
                       x = mono_names[2],
                       y = mono_names[1])
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
# setMethod('predict',  signature(x="flexmix", y="missing"),
#           function(x, mark=NULL, markcol=NULL, col=NULL, 
#                    eps=1e-4, root=TRUE, ylim=NULL, xlim=NULL, main=NULL, xlab=NULL, ylab=NULL,
#                    as.table = TRUE, endpoints = c(-0.04, 1.04), rootogram=F, ...) {
#             
#             
#           })

            
            

internal_boot <- function(fm, R=100){
  
  # Check that fm is a flexmix object
  if(!is(fm, "flexmix")) stop("This method is for FlexMix objects of type: partially linear monotonic regression.")
  
  # Check that fm has at least one monotone obj
  if(!"mon_obj" %in% names(attributes(fm@components[[1]][[1]]))) stop("This method is for FlexMix objects of type: partially linear monotonic regression.")
  
  Y <- fm@model[[1]]@y # store independent data
  dat <- fm@model[[1]]@x # store dependent data
  if("(Intercept)" %in% colnames(fm@model[[1]]@x)){
    x <- x[, colnames(m3@model[[1]]@x) %in% "(Intercept)"] # remove intercept
  }
  
  k <- length(m2@components) # Store the number of components
  coefs <- rep(list(),k)
  sigmas <- matrix(ncol = k, nrow = R)
  mon_fits <- rep(list(rep(NA, R)),k)
  
  for(i in 1:k){ # populate each parameter
    iter <- 1
    while(iter < R ){
      # TODO give part_fit an optional formula argument
      # TODO give part_fit a starting coef and g() argument
      current <- part_fit(formula = fm@model[[1]]@fullformula, start = fm@components[[i]][[1]]) 
      
      coefs[[i]][[iter]] <- current$coef
      sigmas[iter, i] <- current$sigma
      mon_fits[[i]][[iter]] <-current$fitted_pava
      
      # assign 
      iter <- iter + 1
    }
    
  }
  
  
}
            
            