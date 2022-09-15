
rm(list = ls())

library(wsbackfit)


###############################################
# Gaussian Simulated Sample
###############################################
set.seed(123)
# Define the data generating process
n <- 1000
x1 <- runif(n)*4-2
x2 <- runif(n)*4-2
x3 <- runif(n)*4-2
x4 <- runif(n)*4-2
x5 <- as.numeric(runif(n)>0.6)
f1 <- 2*sin(2*x1)
f2 <- x2^2
f3 <- 0
f4 <- x4
f5 <- 1.5*x5

mu <- f1 + f2 + f3 + f4 + f5
err <- (0.5 + 0.5*x5)*rnorm(n)
y <- mu + err
df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = as.factor(x5), y = y)
# Fit the model with a fixed bandwidth for each covariate
m0 <- sback(formula = y ~ sb(x1, h = 0.1) + sb(x2, h = 0.1)
+ sb(x3, h = 0.1) + sb(x4, h = 0.1), kbin = 30, data = df)
summary(m0)
op <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
plot(m0,composed = FALSE)

# Fit the model with bandwidths selectec using K-fold cross-validation
## Not run:
m0cv <- sback(formula = y ~ x5 + sb(x1) + sb(x2)
+ sb(x3) + sb(x4), kbin = 30, bw.grid = seq(0.01, 0.99, length = 30), KfoldCV = 5,
data = df)
summary(m0cv)
par(mfrow = c(2,2))
plot(m0cv)
## End(Not run)
# Estimate Variance as a function of x5 (which is binary)
resid <- y - m0$fitted.values
sig0 <- var(resid[x5 == 0])
sig1 <- var(resid[x5 == 1])
w <- x5/sig1 + (1-x5)/sig0
m1 <- sback(formula = y ~ x5 + sb(x1, h = 0.1) + sb(x2, h = 0.1)
+ sb(x3, h = 0.1) + sb(x4, h = 0.1), weights = w, kbin = 30, data = df)
summary(m1)
par(mfrow