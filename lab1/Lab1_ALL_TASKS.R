###############################################################################
##                  TASK 1                                                   ##
###############################################################################
###############################################################################
##                  TASK 1a                                                  ##
###############################################################################
n = 70
success = 22
failures = n - success
alpha = 8
b = 8

lhs = alpha + success
rhs = b + failures
nDraws = 10000

SD = c()
E = c()
result = c()
for (n in seq(1, nDraws, 1)){
  result[n] = rbeta(n, lhs, rhs) # Accumulating number of n draws
  SD[n] = sd(result) 
  E[n] = mean(result) 
}

plot(E, main = "Expected value as a function of #Draws", 
     xlab = "# of Draws", ylab = "Expected Value E[theta|y]")
# The true value E = alpha / (alpha + beta)
# where our posterior gives (alpha = lhs) and (beta = rhs)
lhs / (lhs+ rhs) # 0.3488

plot(SD, main = "SD for different sample sizes", xlab = "# of Draws",
     ylab = "Standard deviation", ylim = c(0.03,0.08))
# The true value 
# SD = sqrt((alpha * beta) / (((alpha + beta)^2) *(alpha + beta + 1)))
# where our posterior gives (alpha = lhs) and (beta = rhs)
SD_Calc = sqrt((lhs *rhs) / (((lhs + rhs)^2) *(lhs + rhs + 1)))
SD_Calc # = 0.05109
SD[nDraws] # = 0.05116S

###############################################################################
##                  TASK 1b                                                  ##
###############################################################################
# Draw 10000 random values from the posterior to compute the posterior 
# probability Pr(θ > 0.3|y) and compare with the exact value from the Beta 
# posterior. [Hint: use pbeta()].

drawed_values = rbeta(nDraws,lhs,rhs)
over = 0
under = 0
for (value in drawed_values){
  if (value > 0.3){
    over = over +1
  }
  else
  {
    under = under + 1
  }
}
prob_over = over / length(drawed_values)
prob_over

beta_posterior = pbeta(0.3, lhs, rhs)
beta_posterior

###############################################################################
##                  TASK 1c                                                  ##
###############################################################################
# Draw 10000 random values from the posterior of the odds φ = θ / (1−θ) by using
# the previous random draws from the Beta posterior for θ and plot the posterior
# distribution of φ. [Hint: hist() and density() can be utilized].

odds = drawed_values / (1-drawed_values)
hist(odds)
plot(density(odds))

###############################################################################
##                  TASK 2                                                   ##
###############################################################################





###############################################################################
##                  TASK 2c                                                  ##
###############################################################################
monthly_income = c(33, 24, 48, 32, 55, 74, 23, 17)
n = length(monthly_income)
df = n - 1 # Degrees of freedom
mean = 3.6 # From task description
number_of_draws = 10000

# Inverse-chisq is related distribution to rchisq (Inversed)
draws = rchisq(number_of_draws,n)
tao_sq = sum((log(monthly_income)-mean)^2) / n
upper = n * tao_sq
prior_sigma_squared = upper *  (1/draws)
plot(density(prior_sigma_squared), xlim=c(0,3))

###############################################################################
##                  TASK 2b                                                  ##
###############################################################################
s2 = prior_sigma_squared
x = sqrt(s2)/sqrt(2)
cdf = pnorm(sqrt(s2)/sqrt(2))
G = 2 * cdf - 1
plot(density(G),xlim=c(0,.8))

###############################################################################
##                  TASK 2c                                                  ##
###############################################################################
library("bayestestR")
interval = eti(G,0.95)
plot(density(G), col = "blue",xlim=c(0,.8))
abline(v= interval[2], col = "red")
abline(v= interval[3], col = "red")

###############################################################################
##                  TASK 2d                                                  ##
###############################################################################
sorted_G = sort(G)
test = hdi(sorted_G, ci = 0.95)
plot(density(sorted_G), col = "blue",xlim=c(0,.8))
abline(v= test[2], col = "green")
abline(v= test[3], col = "green")

# Comparison
plot(density(sorted_G), col = "blue",xlim=c(0,.8))
abline(v= interval[2], col = "red")
abline(v= interval[3], col = "red")
abline(v= test[2], col = "green")
abline(v= test[3], col = "green")
###############################################################################
##                  TASK 3                                                   ##
###############################################################################
#Given observations Y (in radians)
# pi < y. < pi
Y <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)
#K ~ Exponential(lambda = 0.5)
#k > 0
lambda <- 0.5
#the given u
u <- 2.4
#Fine grid of K:s
grid <- seq(from = 0.01, to = 6, by = 0.01)
#A function for the posterior.
posterior <- function(k)
{
  return ((1/(2 * pi * besselI(k,0)))^10 *
            lambda * exp(k*(sum(cos(y-u)) - lambda)))
}
post <- c()
for (k in grid)
{
  post <- c(post, posterior(k))
}
normalized_post <- post/sum(post)
plot(grid, normalized_post, main = "Posterior distribution of k",
     xlab = "k-value",
     ylab = "P(k | u, y)",
     type ='l')
###############################################################################
##                  TASK 3b                                                  ##
###############################################################################
abline(v = (which.max(normalized_post)/100))