---
output:
  pdf_document: default
  html_document: default
---

```{r}
monthly_income = c(33, 24, 48, 32, 55, 74, 23, 17) #Given
y <- monthly_income
n <- length(y)
mu <- 3.6 #Given
tau_2 <- (sum((log(y) - mu)^2) / n) #Given
nDraws <- 10000 #Given
draw_inv_chi <- function(nDraws, v, s)
{
    draws <- rchisq(nDraws, v)
    return( (v * s)/draws )
}
sigma2 <- draw_inv_chi(nDraws = nDraws, v = n, s = tau_2)

plot(density(sigma2), xlim = c(0,1),
     main = "Density function of sigma2")

################################################################
#b)
#pnorm gives the cummulativ distribution function (phi(x))
CDF <- pnorm( sqrt(sigma2) / sqrt(2) )
G <- ((2 * CDF) - 1)
plot(density(G), main = "Density funciton for Gini-Coef")


################################################################
#c)
lower <- quantile(G, 0.025)
upper <- quantile(G, 0.975)

plot(density(G), main = "Density funciton for Gini-Coef")
abline(v = lower, col = "red")
abline(v = upper, col = "red")
legend("topright", legend = c("Gini-distri", "Quantiles"),
       col = c("black", "red"), lty = 1)
################################################################
#d)
########################################################
## COULD BE DONE LIKE THIS
## library("bayestestR")
## interval = eti(G,0.95)
##
## sorted_G = sort(G)
## test = hdi(sorted_G, ci = 0.95)
##
## AND PLOTTED LIKE THIS
## plot(density(sorted_G), col = "blue",xlim=c(0,.8))
## abline(v= interval[2], col = "red")
## abline(v= interval[3], col = "red")
## abline(v= test[2], col = "green")
## abline(v= test[3], col = "green")
## OR LIKE BELOW
########################################################
kernel <- density(G)
kernel <- data.frame(kernel$y, kernel$x)
#sort the DF based on the Density. 
kernel <- kernel[order(kernel$kernel.y, decreasing = TRUE),]
tot_sum <- sum(kernel$kernel.y)
cum_dens <- 0.00000000000001
i <- 0
while ((cum_dens/tot_sum) < 0.95)
{
  i <- (i + 1)
  cum_dens <- cum_dens + kernel$kernel.y[i]
}

abline(v = max(kernel$kernel.x[1:i]), col = "blue")
abline(v = min(kernel$kernel.x[1:i]), col = "blue")

legend("topright", legend = c("Gini-distri", "Quantiles", "HPDI"),
       col = c("black", "red", "blue"), lty = 1)
```

```{r}
plot(density(G))
```
