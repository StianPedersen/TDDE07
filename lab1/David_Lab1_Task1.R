n <- 70
s <- 22
f <- n - s
alfa <- 8
beta <- alfa

nDraws <- 10000
#The joint posterior of BETA(a,b) with likelihood bern(theta)
#is BETA(alfa ´+ s, beta + f).
posterior_alfa <- alfa + s
posterior_beta <- beta + f
#Draw from posterior using rbeta
posterior_draws <- as.vector(matrix(ncol = nDraws))
SD <- as.vector(matrix(ncol = nDraws))
E <- as.vector(matrix(ncol = nDraws))
for (n in 1:nDraws) 
{
  posterior_draws[n] <- rbeta(1, posterior_alfa, posterior_beta)
  SD[n] <- sd(posterior_draws[1:n])
  E[n] <- mean(posterior_draws[1:n])
}


#E[theta] of a beta distribution = alfa/(beta + alfa)
true_Exp <- posterior_alfa/(posterior_beta + posterior_alfa)
#V[theta] of a beta distribution = (alfa * beta) / ((alfa + beta)^2 * (alfa + beta + 1)
true_Var <- (posterior_alfa * posterior_beta) / ((posterior_alfa + posterior_beta)^2 * (posterior_alfa + posterior_beta + 1))
true_SD <- sqrt(true_Var)
plot(SD, main = "Standard deviation after x iterations",
     xlab = "# of iterations",
     col = "blue",
     type = 'l')
abline(h = true_SD, col = "red")
plot(E, main = "Mean value after x iterations",
     xlab = "# of iterations",
     col = "blue",
     type = 'l')
abline(h = true_Exp, col = "red")

################################################################
#b)

#pbeta returns p(theta < 0.3 | alfa, beta) -> use 1 - pbeta
pr_true <- (1 - pbeta(0.3,posterior_alfa, posterior_beta))

over <- 0 
under <- 0
for (n in 1:nDraws)
{
  if (posterior_draws[n] > 0.3)
  {
    over <- (over + 1)
  }
  else
  {
    under <- (under + 1)
  }
}
pr_sim <- (over / length(posterior_draws))
cat(pr_sim, pr_true)

################################################################
#c)
odds <- (posterior_draws / (1 - posterior_draws))
plot(density(odds))
hist(odds)
