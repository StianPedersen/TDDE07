#1
n = 70
success = 22
failures = n - success
alpha = 8
b = 8
# Posterior = Beta(alpha + success, b + failures)
lhs = alpha + success
rhs = b + failures
nDraws = 10000
######################################################
# a)
SD_list = c()
E_list = c()
result = c()
for (n in seq(1, nDraws, 1)){
  result[n] = rbeta(n, lhs, rhs)
  SD[n] = sd(result)
  E[n] = mean(result)
}

# The true values E = alpha / (alpha + beta)
# where our posterior gives alpha = lhs and beta = rhs
lhs / (lhs+ rhs)

plot(E)

Variance = (lhs *rhs) / (((lhs + rhs)^2) *(lhs + rhs + 1))
SD_Calc = sqrt(Variance)
plot(SD)
SD_Calc
SD[nDraws]

##########################################################
#B
drawed_values = rbeta(nDraws,lhs,rhs)
drawed_values
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
prob_under = under / length(drawed_values)
prob_over
prob_under

temp = pbeta(0.3, lhs, rhs)
temp

##########################################################
#C
odds = drawed_values / (1-drawed_values)
hist(odds)
plot(density(odds))


############################################################
##                  TASK 2                                ##
############################################################
n = 8 
monthly_income = c(33, 24, 48, 32, 55, 74, 23, 17)
mean = 3.6
log(mean(monthly_income))
