data <- readRDS("Precipitation.rds")
logdata <- log(data)

#Init Values
u_0 <- 2
t_0 <- 1
sigma_0 <- 1
v_0 <-  1

#Make a function to draw u 
draw_u <- function(u_0, t_0, sigma2, x_mean, N_of_x)
{
  w <- (N_of_x / sigma2) / ( (N_of_x / sigma2) + (1 / t_0) )
  t_n <- 1 / ( (N_of_x / sigma2) + (1 / t_0) )
  u_n <- (w * x_mean) + ((1 - w) * u_0)
  sample <- rnorm(1, u_n, sqrt(t_n))
  return (sample)
}

draw_inv_chi <- function(v, s)
{
  X <- 1 / rchisq(1,v)
  return (X * (v*s))
}

draw_sigma <- function(v_0, sigma_0, x, u)
{
  v_n <- (v_0 + length(x))
  scaler <- ((v_0 * sigma_0) + sum((x - u)^2)) / v_n
  return(draw_inv_chi(v_n, scaler))
}

#NMBR OF SAMPLES
N <- 1000
#Start values
u_samples <- matrix(ncol = 1, nrow = 1)
sigma_samples <- matrix(ncol = 1, nrow = 1)
u_start <- 3
sigma_start <- sigma_0
u_samples[1] <- u_0
sigma_samples[1] <- sigma_start

data_mean <- mean(logdata)
data_length <- length(logdata)

for (i in 2:(N*2))
{
  #If statement adds the old result for pretty printing 
  if(i %% 2 == 0)
  {
     u_samples[i] <- draw_u(u_0, t_0, sigma_samples[i-1], data_mean, data_length)
     sigma_samples[i] <- sigma_samples[i-1]
  }
  else
  {
    u_samples[i] <- u_samples[i-1]
    sigma_samples[i] <- draw_sigma(v_0, sigma_0, logdata, u_samples[i])  
  }
}

plot(u_samples, sigma_samples, type = 'l', main = "Gibbs sample 1000 iterations")
acf(u_samples)
acf(sigma_samples)
acf_u <- acf(u_samples, plot = FALSE)
acf_sigma2 <- acf(sigma_samples, plot = FALSE)
IF_u <- 1 + 2 * sum(acf_u$acf[-1])
IF_sigma2 <- 1 + 2 * sum(acf_sigma2$acf[-1])

plot(sigma_samples, type = 'l', col = 'red', 
     xlim = c(0, 1000), ylab = "Sample value",
     main = "Trajectories")
lines(u_samples, col = 'blue', type = 'l')
legend("topright", legend = c("sigma", "U"),
       col = c("red", "blue"), lty = 1)
#B
#Discard the first 100 samples when sampling using the draws from a)
y_draws <-  matrix(ncol = (N-100), nrow = 1)
for (i in 1:(N - 100))
{
  u <- u_samples[(2 * i) + 100]
  sd <- sqrt(sigma_samples[(2 * i) + 100])
  y_draws[i] <- exp(rnorm(1, u, sd))
}

#Make a historgram, and add density curves for the data and the posterior predictions.
hist(data, breaks = 20,
     freq = FALSE,
     main = "Histogram and Posterior Predictive Density",
     xlab = "Daily Precipitation", ylim = c(0,0.15))
lines(density(data), col = "red")
lines(density(y_draws), col = "blue")
legend("topright", legend = c("Observed Data", "Posterior Predictive Density"),
       col = c("red", "blue"), lty = 1)
