mu = 13
sigma_sq = 3
big_T = 300
phi = 0.1

x_one = mu

AR_1 <- function(big_T,x_one,sigma_sq,mu,phi)
{
  x_matrix = matrix(0,nrow=big_T,ncol=1)
  x_matrix[1] = x_one
  for (i in 2:big_T)
  {
    e = rnorm(1,0,sigma_sq)
    x = mu + phi*(x_matrix[i-1] - mu) + e
    x_matrix[i] = x
  }
  plot(x_matrix)
  title(c("phi:",phi))
  abline(h=16, col = "blue") # Upper
  abline(h=10, col = "blue") # Lower
  abline(h=13, col = "red") # Mean
  return (x_matrix)
}

# Phi = 0.95
AR_1(big_T, x_one,sigma_sq,mu,phi=0.95)
# Phi = 0.5
AR_1(big_T, x_one,sigma_sq,mu,phi=0.50)
# Phi = 0.1
AR_1(big_T, x_one,sigma_sq,mu,phi=0.10)
# Phi = -0.95
AR_1(big_T, x_one,sigma_sq,mu,phi=-0.95)
# Phi = -0.5
AR_1(big_T, x_one,sigma_sq,mu,phi=-0.50)
# Phi = -0.95
AR_1(big_T, x_one,sigma_sq,mu,phi=-0.10)


# B
x = as.vector(AR_1(big_T, x_one,sigma_sq,mu,phi=0.20))
y = as.vector(AR_1(big_T, x_one,sigma_sq,mu,phi=0.95))

# Treat mu, sigma_sq and phi as unknown. 

library("rstan") # observe startup messages
#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

input_data_y <- list(N = length(y), 
                    y = y)
input_data_x <- list(N = length(y), 
                     y = x)

warmup <- 1000
niter <- 2000

fit_y <- stan(file = 'Lab_3.stan', data = input_data_y, warmup=warmup, iter=niter, chains=4, cores=1)
fit_x <- stan(file = 'Lab_3.stan', data = input_data_x, warmup=warmup, iter=niter, chains=4, cores=1)

#Posterior mean & Credible intervals.
print(fit_y)
print(fit_x)
plot(fit_x)
plot(fit_y)


posterior_samples_y <- as.matrix(fit_y)
posterior_samples_x <- as.matrix(fit_x)

#Effective sample size
library(coda)
effectiveSize(posterior_samples_x)
effectiveSize(posterior_samples_y)


true_muy <-mean(y)
true_mux <-mean(x)
true_phi_y <- 0.95
true_phi_x <- 0.20
true_sigma2_y <- var(y)
true_sigma2_x <- var(x)



plot(posterior_samples_x[,1], type = 'l', main = "posterior for mu, where prev phi = 0.20", col = "blue", ylim = c(7, 14))
abline(h = true_mux, col = "red")
plot(posterior_samples_x[,2], type = 'l', main = "posterior for phi, where prev phi = 0.20", col = "blue")
abline(h = true_phi_x, col = "red")

plot(posterior_samples_y[,1], type = 'l', main = "posterior for mu, where prev phi = 0.95", col = "blue", ylim = c(0, 14))
abline(h = true_muy, col = "red")
plot(posterior_samples_y[,2], type = 'l', main = "posterior for phi, where prev phi = 0.95", col = "blue")
abline(h = true_phi_y, col = "red")
