#install.packages("readxl")
library("readxl")
data <- as.data.frame(readxl::read_xlsx("Linkoping2022.xlsx", col_names = TRUE))
df <- data.frame(data$temp, c(1:365)/ 365)
colnames(df) <- c("temp", "time")


f_of_time <- function(x, beta, epsilon)
{
  return (x%*%beta + epsilon)
}

draw_inv_chi <- function(nDraws, v, s)
{
  draws <- rchisq(nDraws, v)
  return( as.numeric(v * s)/draws )
}
#install.packages("mvtnorm")
library("mvtnorm")
draw_prior_beta <- function(mu_0, sigma2, omega)
{
    return (MASS::mvrnorm(mu = mu_0,Sigma =  sigma2*solve(omega)))
}
#Mean temp for Jan 1st~ 1 degree. 21 Degrees for 31/06.
mu_0 <- c(1, 80, -80) 
#Prior confidence about the intercept is relatively high compared to the other.
omega_0 <- diag(c(2,0.1,0.1),3)
#Scaling parameter when obtaining sigma^2 for beta.
sigma_0 <- 1
#Degrees of freedom when sampling from the invCHISQ.
v_0 <- 1

sigma2 <- draw_inv_chi(nDraws = 1, v_0, sigma_0)


x <- matrix(ncol = 3, nrow = length(df$time))
x[,1] <- 1
x[,2] <- df$time
x[,3] <- df$time^2

temp <- matrix(nrow = 10, ncol = 365)
plot(df$temp)
for (i in 1:10)
{
beta <- draw_prior_beta(mu_0 = mu_0, sigma2 = sigma2, omega = omega_0)
temp[i,] <- t(f_of_time(x, beta, 0))
lines(temp[i,])
}

################################################################
#b)
#Update the parameters to get a posterior as in lecture 5
beta_hat <- solve(t(x)%*%x)%*%t(x)%*%df$temp
mu_n <- solve(t(x)%*%x + omega_0)%*%(t(x)%*%x%*%beta_hat + omega_0%*%mu_0)
omega_n <- t(x)%*%x + omega_0
v_n <- v_0 + length(df$time)
#sigma2_n is used to sample sigma2 from the chi-square
sigma2_n <- (sigma_0*v_0 + ((t(df$temp)%*%df$temp + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n))/v_n)

#Decide #Of draws
nDraws = 10^4

#Sample as many sigmas as there is in ndraws.
#These is then used when sampling from normal distribution to get beta
sigma_draws <- draw_inv_chi(v = v_n, s = sigma2_n, nDraws = nDraws)
#Sample as many betas as there are elements in the sigma vector
sample_betas <- function(u, sigmas, omega)
{
  temp <- matrix(nrow = length(sigmas), ncol = 3)
  for (i in 1:length(sigmas))
  {
    temp[i,] <- MASS::mvrnorm(n= 1, u, (sigmas[i]*solve(omega_n)))
  }
  return (temp)
}
beta_draws <- sample_betas(mu_n, sigma_draws, omega_n)


par(mfrow=c(2,2)) # Create a 2x2 grid of plots to plot the histograms
hist(beta_draws[,1], main="Histogram of beta0 posterior", xlab="beta0")
hist(beta_draws[,2], main="Histogram of beta1 posterior", xlab="beta1")
hist(beta_draws[,3], main="Histogram of beta2 posterior", xlab="beta2")
hist(sigma_draws, main="Histogram of sigma^2 posterior", xlab="sigma^2")
par(mfrow=c(1,1))


#Create a big matrix for the posterior estimations.
Y_big <- matrix(nrow = nDraws, ncol = nrow(df))
for (i in 1:nrow(Y_big))
{
  Y_big[i,] <- f_of_time(x, beta_draws[i,], epsilon = 0)
}

#Create vectors.
median <- vector(mode = "integer", length = nrow(df))
low_perc <- median
up_perc <- median

#Get median, Upper and Lower percentile
for (i in 1:nrow(df))
{
  median[i] <- median(Y_big[,i])
  low_perc[i] <- quantile(Y_big[,i], 0.05)
  up_perc[i] <- quantile(Y_big[,i], 0.95)
}
plot(x = df$time, y=median, type = 'l', col = "blue",
     ylab = "Temperature", xlab = "Time", main = "Plot with median, upper and lower quantile")
points(y = df$temp, x = df$time, col = "grey") 
lines(y = up_perc, x = df$time, col = "red") 
lines(y = low_perc, x = df$time, col = "red") 

#The confidence interval does not need to cover all 
# data points because it describes the uncertainty about f(time).

################################################################
#c)
x_tilde <- -(beta_draws[,2]/beta_draws[,3])*0.5
hist(x_tilde, main="Histogram of x_tilde", xlab="Time")

#Examples of u_0, omega_0 for a greater poly.
#u_10  <- c(1, 80, -80, 0, 0, 0, 0, 0, 0, 0, 0)
#omega_0_10 <- (diag(c(4,1 ,1 , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),11)) 