data <- readRDS("Precipitation.rds")
logdata <- log(data)



#Init Values
u_0 <- 1.5
t_0 <- 1
sigma_0 <- 1.2
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

N <- 1000 #NMBR OF SAMPLES

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

plot(u_samples, sigma_samples, type = 'l')
acf(u_samples)
acf(sigma_samples)
plot(sigma_samples, type = 'l', col = 'red')
lines(u_samples, col = 'blue', type = 'l')

#B
u_start <- 3
sigma_start <- sigma_0

u_samples[1] <- u_0
sigma_samples[1] <- sigma_start

#Start by sample u and sigma.
for (i in 2:(N + 1000))
{
    u_samples[i] <- draw_u(u_0, t_0, sigma_samples[i-1], data_mean, data_length)
    sigma_samples[i] <- draw_sigma(v_0, sigma_0, logdata, u_samples[i])  
}

y_draws <-  matrix(ncol = N, nrow = 1)
for (i in 1:N)
{
  u <- u_samples[i + 1000]
  sd <- sqrt(sigma_samples[i + 1000])
  y_draws[i] <- exp(rnorm(1, u, sd))
}

plot(y = as.dataframe(y_draws), x = length(y_draws))
