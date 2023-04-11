n = 70
success = 22
failures = n - success
alpha = 8
b = 8
# Posterior = Beta(alpha + success, b + failures)
lhs = alpha + success
rhs = b + failures

# a)
for (n in seq(1, 10000, 1)){
  
}