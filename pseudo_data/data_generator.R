# data with no latent categories
################
X <- cbind(
  runif(100, -5, 5),
  runif(100, -10, 10),
  runif(100, -100, 100),
  runif(100, -100, 100),
  runif(100, -100, 100)
)

W <- runif(100, 0, 1)
W[c(10, 20, 30, 40, 50, 60, 70)] <- 0.000000e+00

Y <- (X[,1])^3 - X[,2]^3 - 1.5*X[,3] - 2*X[,4] + X[,5] + rnorm(100, 0, 1000)

pseudo_df_1 <- data.frame(Y, X, W)
names(pseudo_df_1) <- c("Y", "X1", "X2", "X3", "X4", "X5", "W")
# hold <- pseudo_df_1$X1
# pseudo_df_1$X1 <- pseudo_df_1$X2
# pseudo_df_1$X2 <- hold
################


# data with 4 latent categories
################
X <- cbind(
  runif(1000, -5, 5),
  runif(1000, -10, 10),
  runif(1000, -100, 100),
  runif(1000, -100, 100),
  runif(1000, -100, 100)
)

W <- rep(1, 1000)
W[c(10, 20, 30, 40, 50, 60, 70)] <- 0.000000e+00


Y1 <- (X[1:250,1])^3 + 1.5*X[1:250,2] - 1.5*X[1:250,3] - 1*X[1:250,4] + X[1:250,5] + rnorm(250, 0, 3) # category 1
Y2 <- (X[251:500,1])+3 + 3*X[251:500,2] + 2*X[251:500,3] - 2*X[251:500,4] + 2*X[251:500,5] + rnorm(250, 0, 4) # category 2
Y3 <- 2*((X[501:750,1])+5) - 2*X[501:750,2] - 1*X[501:750,3] + 2*X[501:750,4] + 4*X[501:750,5] + rnorm(250, 0, 3) # category 3
Y4 <- 2*((X[751:1000,1])-5) - 3*X[751:1000,2] - 3*X[751:1000,3] - 3*X[751:1000,4] + 3*X[751:1000,5] + rnorm(250, 0, 4) # category 4

pseudo_df_2 <- data.frame(c(Y1, Y2, Y3, Y4), X, W)
names(pseudo_df_2) <- c("Y", "X1", "X2", "X3", "X4", "hello", "W")
hold <- pseudo_df_2$X1
pseudo_df_2$X1 <- pseudo_df_2$X2
pseudo_df_2$X2 <- hold
################


# data with all monotone components
################
X <- cbind(
  runif(1000, -5, 5),
  runif(1000, -10, 10),
  runif(1000, -100, 100),
  runif(1000, -100, 100),
  runif(1000, -100, 100)
)

W <- rep(1, 1000)


Y <- (X[,1])^3 - X[,2]^3 + 2*(X[,3]^3) - X[,4]^3 + 4*X[,5] + rnorm(1000, 0, 100) 
Z <- (X[,1])^5 - X[,2]^3 + rnorm(1000, 0, 1)

mono_df <- data.frame(Y, X, W, Z)
names(mono_df) <- c("Y", "X1", "X2", "X3", "X4", "X5", "W", "Z")
################



# data with 4 latent categories and 2 monotone components
################
Xb <- cbind(
  runif(1000, -5, 5),
  runif(1000, -10, 10),
  runif(1000, -100, 100),
  runif(1000, -100, 100),
  runif(1000, -100, 100)
)

Wb <- rep(1, 1000)
Wb[c(10, 20, 30, 40, 50, 60, 70)] <- 0.000000e+00


Y1b <- (Xb[1:250,1])^3 - 1.5*Xb[1:250,2] - 1.5*Xb[1:250,3] - 1*Xb[1:250,4] + Xb[1:250,5] + rnorm(250, 0, 3) # category 1
Y2b <- (Xb[251:500,1])+3 - Xb[251:500,2]^3 + 2*Xb[251:500,3] - 2*Xb[251:500,4] + 2*Xb[251:500,5] + rnorm(250, 0, 4) # category 2
Y3b <- 2*((Xb[501:750,1])+5) - 2*Xb[501:750,2] - 1*Xb[501:750,3] + 2*Xb[501:750,4] + 4*Xb[501:750,5] + rnorm(250, 0, 3) # category 3
Y4b <- 2*((Xb[751:1000,1])-5) - 3*Xb[751:1000,2] - 3*Xb[751:1000,3] - 3*Xb[751:1000,4] + 3*Xb[751:1000,5] + rnorm(250, 0, 4) # category 4

pseudo_df_3 <- data.frame(c(Y1b, Y2b, Y3b, Y4b), Xb, Wb)
names(pseudo_df_3) <- c("Y", "X1", "X2", "X3", "X4", "hello", "W")
################

# data with 4 latent categories and 1 monotone components; 
# 1 non-monotone component is a factor
# Random NA values in Y, X and W
################

Xb <- cbind(
  runif(1000, -5, 5),
  runif(1000, -10, 10),
  runif(1000, -100, 100),
  runif(1000, -100, 100),
  rbinom(1000, 1, 0.25)
)

Wb <- rep(1, 1000)
Wb[c(10, 20, 30, 40, 50, 60, 70)] <- 0.000000e+00


Y1b <- (Xb[1:250,1])^3 - 1.5*Xb[1:250,2] - 1.5*Xb[1:250,3] - 1*Xb[1:250,4] + Xb[1:250,5] + rnorm(250, 0, 3) # category 1
Y2b <- (Xb[251:500,1])+3 - Xb[251:500,2] + 2*Xb[251:500,3] - 2*Xb[251:500,4] + 2*Xb[251:500,5] + rnorm(250, 0, 4) # category 2
Y3b <- 2*((Xb[501:750,1])+5) - 2*Xb[501:750,2] - 1*Xb[501:750,3] + 2*Xb[501:750,4] + 4*Xb[501:750,5] + rnorm(250, 0, 3) # category 3
Y4b <- 2*((Xb[751:1000,1])-5) - 3*Xb[751:1000,2] - 3*Xb[751:1000,3] - 3*Xb[751:1000,4] + 3*Xb[751:1000,5] + rnorm(250, 0, 4) # category 4

pseudo_df_factor <- data.frame(c(Y1b, Y2b, Y3b, Y4b), Xb, Wb)
names(pseudo_df_factor) <- c("Y", "X1", "X2", "X3", "X4", "fac", "W")

pseudo_df_factor$Y[sample.int(dim(pseudo_df_factor)[1], 3)] <- NA
pseudo_df_factor$X2[sample.int(dim(pseudo_df_factor)[1], 3)] <- NA
pseudo_df_factor$W[sample.int(dim(pseudo_df_factor)[1], 3)] <- NA
################