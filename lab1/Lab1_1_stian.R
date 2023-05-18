s = 22
n = 70
f = n -s
a = 8
b = 8

## A
draws = 1000
sd = c()
m = c()
for(i in 1:draws)
{
  res = rbeta(i, a + s, b + f)
  sd[i] = sd(res)
  m[i] = mean(res)
}

plot(sd)
plot(m)
mean_real = (a+s) / ((a+s)+(b+f))
print(mean_real) # Converges
alpha = a + s
beta = b + f
sd_real = sqrt((alpha * beta) / (((alpha + beta)^2) * (alpha + beta + 1) ))
print(sd_real) # Convverges

## B
drawed_posterior = rbeta(draws, a+s, b+f)
over = 0
under = 0
for (val in drawed_posterior){
  if (val > 0.3)
  {
    over = over + 1
  }
  else
  {
    under = under + 1
  }
} 

prob_over = over / (over + under)
prob_beta = pbeta(0.3, s + a, f+b)


# C
odds = drawed_posterior / (1-drawed_posterior)
hist(odds)
plot(density(odds))
