############################################################
##                  TASK 2                                ##
############################################################
# A 
monthly_income = c(33, 24, 48, 32, 55, 74, 23, 17)
n = length(monthly_income)
df = n - 1 # From lecture 6
mean = 3.6
number_of_draws = 10000

# https://en.wikipedia.org/wiki/Inverse-chi-squared_distribution
# From wikipedia, inverse-chisq is related distribution
draws = rchisq(number_of_draws,n)
tao_sq = sum((log(monthly_income)-mean)^2) / n
upper = n * tao_sq
prior_sigma_squared = upper / draws
plot(density(prior_sigma_squared), xlim=c(0,3))

# B
# "https://en.wikipedia.org/wiki/Log-normal_distribution" 
s2 = prior_sigma_squared
x = sqrt(s2)/sqrt(2)
cdf = pnorm(sqrt(s2)/sqrt(2))
G = 2 * cdf - 1
plot(density(G),xlim=c(0,.8))

# C
library("bayestestR")
interval = eti(G,0.95)
plot(densityG), col = "blue",xlim=c(0,.8))
abline(v= interval[2], col = "red")
abline(v= interval[3], col = "red")

# D
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


