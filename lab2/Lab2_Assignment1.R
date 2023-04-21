######################################################
##                  Task 1                          ##
######################################################
#A
u0 <- t(c(0, 100, -100))
omega0 <- 0.1*I3
sigma2o <- 1
#install.packages("readxl")
df <- readxl::read_xlsx("Linkoping2022.xlsx")

install.packages("mvtnorm")
MASS::mvrnorm(1, 0,1)
