# data with no latent categories
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
 

# data with 4 latent categories
X <- cbind(
  runif(400, -5, 5),
  runif(400, -10, 10),
  runif(400, -100, 100),
  runif(400, -100, 100),
  runif(400, -100, 100)
)

W <- rep(1, 400)

Y1 <- (X[1:100,1])^3 + 1.5*X[1:100,2] - 1.5*X[1:100,3] - 1*X[1:100,4] + X[1:100,5] + rnorm(100, 0, 50) # category 1
Y2 <- (X[101:200,1])+3 + 3*X[101:200,2] + 2*X[101:200,3] - 2*X[101:200,4] + 2*X[101:200,5] + rnorm(100, 0, 60) # category 2
Y3 <- 2*((X[201:300,1])+5) - 2*X[201:300,2] - 1*X[201:300,3] + 2*X[201:300,4] + 4*X[201:300,5] + rnorm(100, 0, 60) # category 3
Y4 <- 2*((X[301:400,1])-5) - 3*X[301:400,2] - 3*X[301:400,3] - 3*X[301:400,4] + 3*X[301:400,5] + rnorm(100, 0, 70) # category 4

pseudo_df_2 <- data.frame(c(Y1, Y2, Y3, Y4), X, W)
names(pseudo_df_2) <- c("Y", "X1", "X2", "X3", "X4", "X5", "W")

# df3 <- gendata()



