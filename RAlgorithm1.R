library(Rcsdp)
library(Matrix)
library(matrixcalc)
library(condMVNorm)

#Number of discretization points
N = 10^4
#Number of variables
d = 3
alpha = 0.99

#Pareto example1
# omega = c(3, 4, 5)
# sd = sqrt(omega/((omega-1)^2*(omega-2)))

#Pareto example2
# omega=2
# Not defined variance

#Normal example
sd = c(1, 2, 3)

#Correlation matrix
#Cor (X_1, X_2) = 0.2
R <- sparseMatrix(i = c(1,2,3,1), j = c(1,2,3,2), x = c(rep(1,d),0.2), dims = c(d,d), symmetric = TRUE)
D = diag(sd)

fix_row = 1
fix_col = 2

#Covariance matrix
S = D%*%R%*%D

#Dinv = diag(c(1/sx, 1/sy, 1/sz), 3, 3)


if (is.positive.definite(round(as.matrix(S),8)))
{
  set.seed(567)
  
  
  #Q=function(x){(1-x)^(-1/omega)-1}
  # X=auxiliary matrix
  #X=matrix(rep(Q(seq(0,1,length=N+1)[1:N]),d),nrow=N)
  
  
  
  #Pareto example1 
  #u = seq(0,1,length=N+1)[1:N]
  #x1 = (1-u)^(-1/omega[1])-1 
  #x2 = (1-u)^(-1/omega[2])-1 
  #x3 = (1-u)^(-1/omega[3])-1
  #X = matrix(c(x1, x2, x3), nrow = N)
  
  #Pareto example2
  #Q = function(x){(1-x)^(-1/omega)-1}
  #X = matrix(rep(Q(0, 1, length = N+1)[1:N], d), nrow = N)
  
  #Normal
  u = seq(0,1,length=N+1)[2:N]
  
  x1 = qnorm(u, sd = sd[1])
  x2 = qnorm(u, sd = sd[2])
  x3 = qnorm(u, sd = sd[3])
  X = matrix(c(x1, x2, x3), nrow = N-1)
  
  ###
  
  Xmax_Cov_Info = RAalgorithm1(d ,X, bound = "max", epsilon = 0.001, Cov = round(as.matrix(S),8), ind_row = fix_row, ind_col = fix_col)  
  Xmin_Cov_Info = RAalgorithm1(d ,X, bound = "min", epsilon = 0.001, Cov = round(as.matrix(S),8), ind_row = fix_row, ind_col = fix_col)  
  Xmax = RAalgorithm1(d ,X, bound = "max", epsilon = 0.001)  
  Xmin = RAalgorithm1(d ,X, bound = "min", epsilon = 0.001)  
  
} else
  print("The correlations introduced do not concord with the variance of the variables")

#Value of the min and max expected shortfall 
Xmin$ES
Xmin_Cov_Info$ES
Xmax_Cov_Info$ES
Xmax$ES

cor(Xmin$X)
cor(Xmin_Cov_Info$X)
cor(Xmax_Cov_Info$X)
cor(Xmax$X)


