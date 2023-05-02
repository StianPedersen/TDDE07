# Assignment 2
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
