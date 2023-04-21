# Assignment 2
# A ######################################
library("mvtnorm")
library("bayestestR")
WomenData = read.table("WomenAtWork.dat",sep = " ",header = T)
Y = WomenData[,1]
X = WomenData[,2:8]

n = 132
tao = 2
I_matrix = diag(7)

# From lecture 6 code     
logistic_posterior <- function(initVal, y, X, mu, Sigma)
{
  linPred <- as.matrix(X)%*%initVal;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  
  # density "d", "mv" multivariate, initval???
  logPrior <- dmvnorm(initVal, mu, Sigma, log=TRUE);
  return(logLik + logPrior)
}

initVal <- rnorm(7)
mu <- as.matrix(rep(0,7)) # Prior mean vector
Sigma = (tao^2)*I_matrix
prior = dnorm(0,(tao^2)*I_matrix)

OptimRes <- optim(initVal,logistic_posterior,gr=NULL,Y,X,mu,Sigma,method=c("BFGS"),
                  control=list(fnscale=-1),hessian=TRUE)

print("Posterior mode is")
OptimRes$par

approxPostStd <- sqrt(diag(solve(-OptimRes$hessian)))
print('The approximate posterior standard deviation is:')
print(approxPostStd)

x_hat = mean(WomenData$NSmallChild)
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

# The confidence interval does not contain 0 which indicates that the var
# is important. 
# The comparison is almoost the same which means that the approximation is 
# reasonable

# B ###################################################
approxPostStd
OptimRes$par

beta = rmvnorm(10000, mean = OptimRes$par, sigma = -solve(OptimRes$hessian))
beta_t = t(beta)
x = c(1,18,11,7,40,1,1)
beta
prob_1 = exp(t(as.matrix(x)) %*% t(beta)) / (1 + exp(t(as.matrix(x)) %*% t(beta)))
prob_0 = 1 - prob_1 
plot(density(prob_0))

# C ############################################
b_i = rmvnorm(10000, mean = OptimRes$par, sigma = -solve(OptimRes$hessian))
p = exp(t(as.matrix(x)) %*% t(b_i)) / (1 + exp(t(as.matrix(x)) %*% t(b_i)))
bino_simulation = rbinom(13, 1,p)
bino_simulation
sum_0 = sum(bino_simulation == 0)
sum_0
sum_1 = sum(bino_simulation == 1)
barplot(c(sum_0, sum_1),col=c("Red","Blue"),legend=c("Not working", "working"),ylim=c(0,15))
print("Posterior probability")
print(sum_0/(sum_0 + sum_1))
