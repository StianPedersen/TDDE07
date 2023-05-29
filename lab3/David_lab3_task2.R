################################################################
#a)
data = read.table('eBayNumberOfBidderData.dat', header = TRUE)
data_shaved = data
#Remove intercept (const = 1)
data_shaved[2] = NULL

#create a poisson regression model, then print it
a_model = glm(nBids~., data = data_shaved, family = "poisson")
a_model
summary(a_model)
# Significant covariates has a low Pr(>|z|) value
# Most impactful variable is MinBidShare

################################################################
#b)
library("mvtnorm")
library("bayestestR")

poisson_posterior <- function(initVal, y, X, mu, Sigma)
{
  lambda <- exp(as.matrix(X)%*%initVal);
  loglik <- sum(dpois(y,lambda, log = TRUE)); 
  logPrior <- dmvnorm(initVal, mu, Sigma, log=TRUE); # Denisty multivariate normal
  return(loglik + logPrior) 
  
}


initval = rnorm(9)
Y  = data$nBids
X  = data[,2:10]
mean = as.matrix(rep(0,9))
Sigma = 100* solve(t(as.matrix(X))%*%as.matrix(X))

OptimRes <- optim(initval,poisson_posterior,gr=NULL,Y,X,mean,Sigma,method=c("BFGS"),
                  control=list(fnscale=-1),hessian=TRUE)
OptimRes$par
a_model$coefficients

################################################################
#c)
hessian = solve(-OptimRes$hessian)
post_mode = OptimRes$par
posterior_beta_theta0 = rmvnorm(1, post_mode, hessian)

LogPostPoissonBeta <- function(beta, Sigma, mu, y, x)
{
  lambda_k  = (y * (beta %*% t(as.matrix(x))))
  e_lambda = (exp(beta %*% t(as.matrix(x)) ))
  k_factorial = (log(factorial(y)) )
  loglik <-sum( lambda_k - e_lambda - k_factorial)
  logprior <- dmvnorm(beta,mu,Sigma, log = TRUE)
  return (loglik + logprior)
}


RWMSampler <- function(funct_object, iterations, c, Hessian, Post_mode, Sigma, mu, Y,X)
{
  sample_matrix = matrix(1, iterations, 9)
  not_acc = 0
  acc = 1
  for (i in 2:iterations)
  {
    logsample_prev_p = funct_object(sample_matrix[i-1,],Sigma, mu,Y,X)
    
    sample = rmvnorm(1, mean = sample_matrix[i-1,], sigma = Hessian * c)
    logsample_p = funct_object(sample, Sigma, mu, Y,X)
    alpha = min(1,exp(logsample_p-logsample_prev_p))
    random = runif(1, min=0, max=1)
    if(alpha < random)
    {
      sample_matrix[i,] = sample_matrix[i-1,]
      not_acc = not_acc + 1
    }
    else
    {
      sample_matrix[i,] = sample
      acc = acc + 1
    }
    
  }
  
  print(acc/(acc+not_acc))
  return (sample_matrix)
}
#Tuning Parameter
c=0.1
iterations = 10000
sample_matrix<- RWMSampler(LogPostPoissonBeta, iterations,c,hessian,post_mode, Sigma, mean, Y,X)

for(i in 1:9){
  plot(sample_matrix[,i], type = "l")
  abline(h=post_mode[i],col = "blue")
}
# Stabelizes after 2000 ish
# Seems to converge to values we got from optim.

################################################################
#d)
characteristics = c(0, 1,0,1,0,1,0,1.2,0.8) 
after_burnout <- sample_matrix[2001:10000,] # Decided from C 
lambda = exp(after_burnout %*% characteristics) 
dist = rpois(1000, lambda) 
no_bidders = 0 
bidders = 0 
for (i in 1:length(dist)) 
{ 
  if(dist[i] == 0) 
  { 
    no_bidders = no_bidders + 1 
  } 
  else 
  { 
    bidders = bidders + 1 
  } 
} 
#Delen genom det hela  -> procent.
print(no_bidders/(no_bidders+bidders))
