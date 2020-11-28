
# definition of monotone regression model.

mono_reg <- function (formula = .~., diagonal = TRUE) {
  retval <- new("FLXMC", weighted = TRUE,
                # TODO is dist param necessary?
                formula = formula, dist = "mvnorm",
                name = "my model-based clustering") 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  # para must have: coef attribute, sigma attribute, cov attribute, df attribute, ... 
  # ... all must be defined by fit() function
  retval@defineComponent <- function(para, fitted_pava) {
                  # @logLik: A function(x,y) returning the log-likelihood for observations in matrices x and y
                  logLik <- function(x, y) { 
                    dnorm(y, mean=predict(x, ...), sd=para$sigma, log=TRUE)
                  }
                  # @predict: A function(x) predicting y given x. 
                  # TODO x must be partitioned into linear and monotone covars
                  predict <- function(x, x_mon) {
                    p <- (x %*% para$coef) + get_pred(fitted_pava, x_mon)
                    p
                  }
                  new("FLXcomponent", parameters =
                        list(center = para$center, cov = para$cov),
                      df = para$df, logLik = logLik, predict = predict)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  #TODO inputs are x = covariates, y = DV, w = weights. Must have specification of covars with monotone relationship
  retval@fit <- function(x, y, w, ...) {
                  fit <- part_fit(x, y, w, ...)
                  
                  df <- some_number
                  
                  retval@defineComponent(fit$para, fit$fitted_pava, df = df, ...) 
                  }
  retval 
  }



y <- (1:20) + rnorm(20, sd = 3)
ystar <- pava(y, long.out = T, stepfun = T)
plot(y)
lines(ystar$y,type='s')
# Decreasing order:
z <- NULL
for(i in 4:8) {
  z <- c(z,rep(8-i+1,i)+0.05*(0:(i-1)))
}
zstar <- pava(z,decreasing=TRUE)
plot(z)
lines(zstar,type='s')
# Using the stepfunction:
zstar <- pava(z,decreasing=TRUE,stepfun=TRUE)
plot(z)
plot(zstar,add=TRUE,verticals=FALSE,pch=20,col.points="red")
