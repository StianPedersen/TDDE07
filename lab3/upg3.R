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
x = AR_1(big_T, x_one,sigma_sq,mu,phi=0.20)
y = AR_1(big_T, x_one,sigma_sq,mu,phi=0.95)

# Treat mu, sigma_sq and phi as unknown
library("rstan")
Model = '
  data {
    int<lower=0> N;
    vector[N] x;

    
  }
  parameters {
    real mu;
    real phi;
    real<lower=0> sigma;
  }
  model {
    for (n in 2:N)
      y[n] ~ normal(mu + phi * y[n-1], sigma);
    }

  '

