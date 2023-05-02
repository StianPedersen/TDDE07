######################################################
##                  Task 1                          ##
######################################################
#A
#Mean temp for Jan 1st~ 1 degree. 21 Degrees for 31/06.
u0 <- c(1, 80, -80) 
#Prior confidence about the intercept is relatively high compared to the other.
omega0 <- diag(c(2,0.1,0.1),3)
#Scaling parameter when obtaining sigma^2 for beta.
sigma2o <- 1
#Degrees of freedom when sampling from the invCHISQ.
v0 <- 2
#install.packages("readxl")
df <- readxl::read_xlsx("Linkoping2022.xlsx", col_names = TRUE)
days <- c(1:365)
time <- days/365 #x
df_with_time <- data.frame(df, time) #temp = y, time = x

#Function to calculate time.
calc_temp <- function(x,beta, epsilon)
{
 return (x%*%beta + epsilon) 
}

#Create the matrix for x where x^1 = time
x = matrix( data = c(df_with_time$time^0, df_with_time$time^1, df_with_time$time^2),
            nrow = 365, ncol =3, byrow = FALSE)
colnames(x) <- c('x0', 'x1', 'x2')

#Generate epsilon
epsilon <- rnorm(n = 365, mean = 0, sd = 1)

#Function to generate sigma2 ~ scaled-inv-X^2 (v0,sigma2o)
#if X ~ inv-X^2(v) then 1/X ~ X^2(v).
#if X ~ inv-X^2(v) then X*(V*S) X ~ scaled-inv-X^2(v,sigma2o)
draw_sigma2 <- function(v, s)
{
  X <- 1 / rchisq(1,v)
  return (X * (v*s))
}

plot(x = df_with_time$time, y = df_with_time$temp, 
     col = "grey", type = 'l', main = "Regression curves for the prior",
     xlab = "Time", ylab = "Temperature")



#install.packages("mvtnorm")
for (i in 1:10)
{  
#Sample a sigma2 ~scaled-inv-X^2(v0, sigma2o)
sigma <- draw_sigma2(v0, sigma2o)
#Sample beta ~ N(u0, sigma2*omega0)
beta <- MASS::mvrnorm(n = 1, mu = u0, sigma*solve(omega0))
yyy <- calc_temp(x, beta, epsilon)
lines(x = df_with_time$time, y = yyy, col = "blue")
}

#########################################################################
#B
#Update the parameters to get a posterior as in lecture 5
beta_hat <- solve(t(x)%*%x)%*%t(x)%*%df_with_time$temp
u_n <- solve(t(x)%*%x + omega0)%*%(t(x)%*%x%*%beta_hat + omega0%*%u0)
omega_n <- t(x)%*%x + omega0
v_n <- v0 + length(df_with_time$time)
sigma2_n <- (sigma2o*v0 + ((t(df_with_time$temp)%*%df_with_time$temp + t(u0)%*%omega0%*%u0 - t(u_n)%*%omega_n%*%u_n))/v_n)

#Sample as many sigmas as there is in ndraws.
sample_sigma2s <- function(v, s, ndraws)
{
  temp <- vector(mode = "integer", length = ndraws)
  for (i in 1:ndraws)
  {
    temp[i]<- draw_sigma2(v,s)
  }
  return (temp)
}

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



nDraws <- 10^4 #10 000 samples for the poserior parameters.
sigma_posts <- sample_sigma2s(v_n, sigma2_n, nDraws)
beta_posts <- sample_betas(u_n, sigma_posts, omega_n)


par(mfrow=c(2,2)) # Create a 2x2 grid of plots to plot the histograms
hist(beta_posts[,1], main="Histogram of beta0 posterior", xlab="beta0")
hist(beta_posts[,2], main="Histogram of beta1 posterior", xlab="beta1")
hist(beta_posts[,3], main="Histogram of beta2 posterior", xlab="beta2")
hist(sigma_posts, main="Histogram of sigma^2 posterior", xlab="sigma^2")

#Create a big matrix for the posterior estimations.
Y_big <- matrix(nrow = nDraws, ncol = nrow(df_with_time))
for (i in 1:nrow(Y_big))
{
  Y_big[i,] <- calc_temp(x, beta_posts[i,], epsilon)
}

#Create vectors.
median <- vector(mode = "integer", length = nrow(df_with_time))
low_perc <- median
up_perc <- median

#Get median, Upper and Lower percentile
for (i in 1:nrow(df_with_time))
{
  median[i] <- median(Y_big[,i])
  low_perc[i] <- quantile(Y_big[,i], 0.05)
  up_perc[i] <- quantile(Y_big[,i], 0.95)
}
dev.off()
plot(x = df_with_time$time, y=median, type = 'l', col = "blue",
     ylab = "Temperature", xlab = "Time", main = "Plot with median, upper and lower quantile")
points(y = df_with_time$temp, x = df_with_time$time, col = "grey") 
lines(y = up_perc, x = df_with_time$time, col = "red") 
lines(y = low_perc, x = df_with_time$time, col = "red") 
####################################################################
#C
x_tilde <- -(beta_posts[,2]/beta_posts[,3])*0.5
hist(x_tilde, main="Histogram of x_tilde", xlab="Time")

#Examples of u_0, omega_0 for a greater poly.
#u_10  <- c(1, 80, -80, 0, 0, 0, 0, 0, 0, 0, 0)
#omega_0_10 <- (diag(c(4,1 ,1 , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),11)) 