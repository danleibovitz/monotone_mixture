X <- cbind(
  runif(100, -5, 5),
  runif(100, -10, 10),
  runif(100, -100, 100),
  runif(100, -100, 100),
  runif(100, -100, 100)
)

W <- rep(1, 100)

Y <- (X[,1])^3 + 1.5*X[,2] - 1.5*X[,3] - 2*X[,4] + X[,5] + rnorm(100, 0, 100)

pseudo_df_1 <- data.frame(Y, X, W)


# df2 <- gendata()
# df3 <- gendata()