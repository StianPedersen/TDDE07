//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
    int<lower=0> N;
    vector[N] y;
}
parameters {
  real mu;
  real<lower=-1,upper=1> phi;
  real<lower=0> sigma2;
}
model {
  mu ~ normal(13,10);
  phi ~ uniform(-1, 1);
  sigma2 ~ scaled_inv_chi_square(1,1);
  
  for (n in 2:N)
  {
    y[n] ~ normal(mu + phi * y[n-1], sqrt(sigma2));
  }
}
