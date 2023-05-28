y <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)
mu <- 2.4
lambda <- 0.5
k <- 6
denumerator <- function(k, n)
{
  return ( (2 * pi * besselI(k, 0))^n )
}
numerator <- function(k, y, mu, lambda)
{
  return (lambda * exp(k * (sum(cos(y- mu)) - lambda) ) )
}
#f(k)
posterior_draw_k <- function(lambda, k, mu, y)
{
  result <- c()
  for (k in seq(from = 0, to = k, by = 0.01))
  {
    result <- c(result, (numerator(k = k, y = y, mu = mu, lambda = lambda) / denumerator(k = k, n = length(y)) ))
  }
  return (result)
}

draws <- posterior_draw_k(k = k, y = y, mu = mu, lambda = lambda)
draws <- (draws * 100) / sum(draws) #Multiply by 100 to cancel out the grid stepsize
plot(x = seq(0, k, 0.01), y = draws,
     xlab = "k - value", ylab = "Density", 
     main = "Plot of the posterior density and mode",
     type = 'l')
posterior_mode <- max(draws)
abline(v = max(which.max(draws) / 100), col = "red")
legend("topright", legend = c("Posterior", "Mode"),
       col = c("black", "red"), lty = 1)
