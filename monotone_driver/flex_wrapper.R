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


# overwrite method for plot.flexmix
setMethod('plot',  signature(x="flexmix", y="missing"),
          function(x, mark=NULL, markcol=NULL, col=NULL, 
                   eps=1e-4, root=TRUE, ylim=NULL, xlim=NULL, main=NULL, xlab=NULL, ylab=NULL,
                   as.table = TRUE, endpoints = c(-0.04, 1.04), rootogram=F, ...) {
           
  if(is(x@components[[1]][[1]], "FLX_monoreg_component")){
            if(dim( x@components[[1]][[1]]@mon_obj )[2] > 1){ # get dimension of monotone components by reading columns of fitted_pava object
              temp <- list()
              for(i in 1:dim( x@components[[1]][[1]]@mon_obj )[2]){
                holder <- ggplot() +
                  geom_line(aes(x = x@components[[1]][[1]]@mon_obj[,i]$x,
                                y = x@components[[1]][[1]]@mon_obj[,i]$yf)) + 
                  theme_bw() +
                  labs(title = paste(append_suffix(i), " Monotone Regression"),
                       x = "X",
                       y = "Y") 
                
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
                
                if(length(x@components) > 1){
                  
                  pal <- brewer.pal(length(x@components), "Paired")
                  
                  for(j in 2:length(x@components)){
                    dat <- data.frame(x = x@components[[j]][[1]]@mon_obj[,i]$x, 
                                      yf = x@components[[j]][[1]]@mon_obj[,i]$yf)
                    colour <- pal[j]
                    holder <- holder + geom_line(data=dat, aes(x = x, 
                                                           y = yf),
                                             color=colour
                    )
                  }
                }
                
                temp[[i]] <- ggplotGrob(holder)
              }
              return(grid.arrange(grobs=temp, ncol=1))
              # return(grid.arrange(grobs=temp, ncol=1))
            }
            else{
              temp <- ggplot() +
                geom_line(aes(x = x@components[[1]][[1]]@mon_obj[,1]$x, y = x@components[[1]][[1]]@mon_obj[,1]$yf)) + 
                theme_bw() +
                labs(title = "Monotone Component",
                     x = "X",
                     y = "Y") 
              
              if(!is.null(ylim)){
                temp <- temp + ylim(ylim)
              }
              if(!is.null(xlim)){
                temp <- temp + xlim(xlim)
              }
              if(!is.null(ylab)){
                temp <- temp + ylab(ylab)
              }
              if(!is.null(xlab)){
                temp <- temp + xlab(xlab)
              }
              if(!is.null(main)){
                temp <- temp + ggtitle(main)
              }
              
              if(length(x@components) > 1){
                
                pal <- brewer.pal(length(x@components), "Paired")
                
                for(i in 2:length(x@components)){
                  dat <- data.frame(x = x@components[[i]][[1]]@mon_obj[,1]$x, 
                                    yf = x@components[[i]][[1]]@mon_obj[,1]$yf)
                  colour <- pal[i]
                  temp <- temp + geom_line(data=dat, aes(x = x, 
                                               y = yf),
                                           color=colour
                                           )
                }
              }
              
              
              return(grid.arrange(temp, ncol=1))
              # temp <- ggplotGrob(temp)
            }
  }
            
            # hh <- ggplotGrob(flexmix::plot(x))

            # groblist <- append(groblist, plot(as(x, "flexmix"))) # append rootogram
            # 
            # grid.arrange(grobs=groblist, ncol=1)
            # return(grid.arrange(grobs=c(temp #,hh
                                        # ), ncol=1))
          }          
)
            
            
            

# TODO write predict method for monoreg flexmix objects
# setMethod('predict',  signature(x="flexmix", y="missing"),
#           function(x, mark=NULL, markcol=NULL, col=NULL, 
#                    eps=1e-4, root=TRUE, ylim=NULL, xlim=NULL, main=NULL, xlab=NULL, ylab=NULL,
#                    as.table = TRUE, endpoints = c(-0.04, 1.04), rootogram=F, ...) {
#             
#             
#           })

            
            
            
            