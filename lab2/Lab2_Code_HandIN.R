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


######################################################
##                  Task 2                          ##
######################################################
# A ######################################
library("mvtnorm")
library("bayestestR")
WomenData = read.table("WomenAtWork.dat",sep = " ",header = T)
Y = WomenData[,1] # Work or No work (0, 1)
X = WomenData[,2:8] # Attributes

tao = 2 # From task
I_matrix = diag(7) # 7 variables, needed for b ~ N(*, *)

# Code from lecture 6, copypaste   
logistic_posterior <- function(initVal, y, X, mu, Sigma)
{
  linPred <- as.matrix(X)%*%initVal;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) ); # Likelihood
  logPrior <- dmvnorm(initVal, mu, Sigma, log=TRUE); # Denisty multivariate normal, prior
  return(logLik + logPrior) 
}

initVal = rnorm(7) # Init the 7 random variable thetas using the normal distribution
mu <- as.matrix(rep(0,7)) # Prior mean vector is 0 from the task description
Sigma = (tao^2)*I_matrix # Sigma = the right value in b ~ N(*, *) from task description

# Send into the function from lecture 6
OptimRes <- optim(initVal,logistic_posterior,gr=NULL,Y,X,mu,Sigma,method=c("BFGS"),
                  control=list(fnscale=-1),hessian=TRUE)

# Posterior and approx posteriod std deviation from optim, same way as in 
# lecture 6 vode
print("Posterior mode is")
OptimRes$par 
approxPostStd <- sqrt(diag(solve(-OptimRes$hessian)))
print('The approximate posterior standard deviation is:')
print(approxPostStd)

x_hat = mean(WomenData$NSmallChild) # Get the mean from the NSmallChild data
# 95% Confidence interval calculations from wikipedia 
ci = (1.96 * approxPostStd[6])/sqrt(nrow(WomenData)) 
ci_upper = x_hat + ci
ci_lower = x_hat - ci
ci_upper
ci_lower

# Comparison
glmModel<- glm(Work ~ 0 + ., data = WomenData, family = binomial)
print("Comparison")
print(glmModel$coefficients)
print(OptimRes$par)

# The confidence interval does not contain 0 which indicates that the variable
# is important. The values are almoost the same which means that the  
# approximation is reasonable

# B ###################################################
# draw from the random multivariate norm using the means from the optim.
# The inverse hessian is the sigma for the approximation for the parameter vector
# b
get_probability<- function(m,s){
  x = c(1,18,11,7,40,1,1)
  beta = rmvnorm(1000, mean =m , sigma = s)
  prob_1 = exp(t(as.matrix(x)) %*% t(beta)) / (1 + exp(t(as.matrix(x)) %*% t(beta)))
  # Logistic regression gives prob = 1, we want prob = 0 (1 - prob=1)
  prob_0 = 1 - prob_1 
}

# Probability logistic regression model from task description
prob_0 = get_probability(OptimRes$par,-solve(OptimRes$hessian) )

plot(density(prob_0))

# C ############################################
get_probability<- function(m,s){
  x = c(1,18,11,7,40,1,1)
  beta = rmvnorm(1000, mean =m , sigma = s)
  prob_first_woman = exp(t(as.matrix(x)) %*% t(beta)) / (1 + exp(t(as.matrix(x)) %*% t(beta)))
  prob_first_woman = 1 - prob_first_woman 
  probability_sum =   prob_first_woman
  for(i in 1:12)
  {
    beta = rmvnorm(1000, mean =m , sigma = s)
    prob_loop = exp(t(as.matrix(x)) %*% t(beta)) / (1 + exp(t(as.matrix(x)) %*% t(beta)))
    prob_loop = 1 - prob_loop 
    probability_sum = probability_sum + prob_loop
    
  }
  ret_value = probability_sum/13
  ret_value
}
beta_c = get_probability(OptimRes$par,-solve(OptimRes$hessian) )
plot(density(beta_c))
