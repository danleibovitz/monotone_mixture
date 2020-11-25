
# definition of monotone regression model.

mono_reg <- function (formula = .~., diagonal = TRUE) {
  retval <- new("FLXMC", weighted = TRUE,
                # TODO is dist param necessary?
                formula = formula, dist = "mvnorm",
                name = "my model-based clustering") 
  
  # @defineComponent: Expression or function constructing the object of class FLXcomponent
  retval@defineComponent <- function(para) {
                  # @logLik: A function(x,y) returning the log-likelihood for observations in matrices x and y
                  logLik <- function(x, y) { 
                    dnorm(y, mean=predict(x, ...), sd=para$sigma, log=TRUE)
                  }
                  # @predict: A function(x) predicting y given x. 
                  # TODO x must be partitioned into linear and monotone covars
                  predict <- function(x) {
                    p <- (x %*% para$coef) + pava_density(x)
                    p
                  }
                  new("FLXcomponent", parameters =
                        list(center = para$center, cov = para$cov),
                      df = para$df, logLik = logLik, predict = predict)
  }
  
  # @fit: A function(x,y,w) returning an object of class "FLXcomponent"
  #TODO inputs are x = covariates, y = DV, w = weights. Must have specification of covars with monotone relationship
  retval@fit <- function(x, y, w, ...) {
                  para <- cov.wt(y, wt = w)[c("center", "cov")] 
                  df <- (3 * ncol(y) + ncol(y)^2)/2
                  if (diagonal) {
                    para$cov <- diag(diag(para$cov)) 
                    df <- 2 * ncol(y)
                  }
                  retval@defineComponent(c(para, df = df)) }
  retval }