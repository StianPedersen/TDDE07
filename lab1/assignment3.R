############################################################
##                    TASK 3                              ##
############################################################
#Given observations Y (in radians)
# pi < y. < pi
Y <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)
#K ~ Exponential(lambda = 0.5)
#k > 0
lambda <- 0.5
#the given u
u <- 2.4
#Fine grid of K:s
grid <- seq(from = 0.01, to = 10, by = 0.01)
#A function for the posterior.
posterior <- function(k)
{
  return ((1/(2 * pi * besselI(k,0)))^10 * lambda * exp(k*(sum(cos(y-u)) - lambda)))
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
#B
abline(v = (which.max(normalized_post)/100))

