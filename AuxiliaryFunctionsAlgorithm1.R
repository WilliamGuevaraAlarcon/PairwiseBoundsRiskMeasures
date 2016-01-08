ReplaceFixedConstrains <- function(Cov_X_aug, S, ind_row, ind_col)
{
  for (l in 1:length(ind_row))
  {
    Cov_X_aug[ind_row[l],ind_col[l]] = S[ind_row[l],ind_col[l]]
    Cov_X_aug[ind_col[l],ind_row[l]] = S[ind_col[l],ind_row[l]]
  }
  
  return(Cov_X_aug)
}

BuildPairOfIndexes <- function(ind_row, ind_col)
{
  W = matrix( , nrow = length(ind_row), ncol = 2) 
  
  for (l in 1:length(ind_row))
  {
    if (ind_row[l] > ind_col[l]) 
      W[l,] = c(ind_col[l],ind_row[l])
    else
      W[l,] = c(ind_row[l],ind_col[l])
  }
  
  return(W)
}

MakeVariableCovZero <- function(X, k = 1)
{
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (((i == k) || (j == k)) & (i != j))
        X[i,j] = 0
      
      return(X)
}

MakeSubDiagonalFixed <- function(X, D)
{
  d = nrow(D)
  NewDiag = c(diag(D),X[(d+1),(d+1)])
  diag(X) = NewDiag
  
  return(X)
}

BuildConstraintsByEntry <- function(S, k)
{
  d1 <- nrow(S)
  numpairs <- d1*(d1 + 1)/2 - d1 + 1
  
  A <-list()
  b <- rep(0, numpairs)
  
  for (l in 1:(d1))
  {
    A[l] = list(list(simple_triplet_sym_matrix(i = c(l), j = c(l), v = c(1), n = d1)))
    b[l] = S[l,l]
  }
  
  ind = 1:(d1)
  ind_k = ind[!ind %in% k]  
  
  l = d1+1
  for (i in ind_k)
    for (j in ind_k[ind_k > i])
    {
      A[l] = list(list(simple_triplet_sym_matrix(i = c(i,j), j = c(j,i), v = c(1), n = d1)))
      b[l] = 2*S[i,j]
      l = l+1
    }
  
  return(list(A = A, b = b))
}

IsInRows <- function(v, M)
{
  for (i in 1:nrow(M))
    if (identical (v, M[i,]))
    {  
      return(T) 
      stop()
    }
  return(F)
}

BuildConstraintsWithFixedInfo <- function(S, k, row_fix, col_fix)
{
  d1 <- nrow(S)
  numpairs <- d1*(d1 + 1)/2
  
  A <-list()
  
  A[1] = list(list(simple_triplet_sym_matrix(i = c(1), j = c(1), v = c(1), n = d1)))
  b = S[1,1]
  
  for (l in 2:(d1))
  {
    A[l] = list(list(simple_triplet_sym_matrix(i = c(l), j = c(l), v = c(1), n = d1)))
    b = append(b, S[l,l])
  }
  
  ind = as.numeric(1:(d1))
  
  l = d1+1
  
  IndFix = BuildPairOfIndexes(row_fix, col_fix)
  
  
  for (i in ind)
  {
    for (j in ind[ind > i])
    {
      if ((i != k & j!= k) |  IsInRows(c(i,j), IndFix))
      {  
        A[l] = list(list(simple_triplet_sym_matrix(i = c(i,j), j = c(j,i), v = c(1), n = d1)))
        b = append(b, 2*S[i,j])
        
        l = l+1
      }
    }
  }
  
  return(list(A = A, b = b))    
  
}


FindNegativeOpposedConstraint <- function(C)
{
  vneg = rep(1,length(C$v))
  
  for (l in 1:length(C$v))
    if (C$i[l] != C$j[l])
      vneg[l] = -1
    
    return(simple_triplet_sym_matrix(i = C$i, j = C$j, v = vneg))    
}

FindBoundVariance = function (A, b, C)
{
  K <- list(type="s",size = C$n)
  
  C1 = as.matrix(C)
  C2 = as.matrix(FindNegativeOpposedConstraint(C))
  
  Sol1 = csdp(list(C1),A,b,K)
  Sol2 = csdp(list(C2),A,b,K)
  
  X1=Sol1$X[[1]]
  X2=Sol2$X[[1]]
  objvar1 = sum(diag(C1%*%X1))
  objvar2 = sum(diag(C1%*%X2))
  status1 = Sol1$status
  status2 = Sol2$status
  
  if (objvar2 <= objvar1)
    Sol = list(Xmin = X2, Xmax = X1, min = objvar2, max = objvar1, statusmin = status2, statusmax = status1)    
  else 
    Sol = list(Xmin = X1, Xmax = X2, min = objvar1, max = objvar2, statusmin = status1, statusmax = status2)    
  
  return(Sol)
}

RAalgorithm1 <- function (d, X, alpha = 0.99, epsilon = 0.001, bound = "min", Cov = NA, ind_row = NA, ind_col = NA)
{
  ES = Inf
  # delta = difference between two consecutive estimates 
  delta = Inf
  
  N = nrow(X)
  # iteratively rearrange each column of the matrix until the difference between
  # two consecutive estimates of the best ES is less than epsilon
  
  #Iterator for the columns of the matrix
  while(delta > epsilon )
  {
    
    #EStemp = auxiliary variable
    EStemp = ES  
    
    Indices = seq(1,d)
    
    for (k in Indices)
    {
      ind = Indices
      ind = ind[-k]
      
      S_k = rowSums(X[,ind])  
      
      X_aug = cbind(X,S_k)  
      Cov_X_aug = cov(X_aug)  
      
      if (is.na(ind_row))
      {
        Cov_X_aug = MakeVariableCovZero(Cov_X_aug, k)
        #Cov_X_aug = MakeSubDiagonalFixed(Cov_X_aug, S)
      } else  
        Cov_X_aug = ReplaceFixedConstrains(Cov_X_aug, Cov, ind_row, ind_col)
      
      Var_to_min = simple_triplet_sym_matrix(i = c(k, k, d+1, d+1), j = c(k, d+1, k, d+1), v = c(1,1,1,1))
      
      if (is.na(ind_row))
        Const = BuildConstraintsByEntry(Cov_X_aug, k)
      else  
        Const = BuildConstraintsWithFixedInfo(Cov_X_aug, k, ind_row, ind_col)
      
      Solution = FindBoundVariance (Const$A, Const$b, Var_to_min)
      
      Indices1 = c(Indices,d+1)
      Indices1 = Indices1[!Indices1 %in% k]
      
      mean_cond = rep(0,N)  
      mean_marg = colMeans(X_aug)
      
      if (bound == "min")
      {
        var_cond = condMVN(mean = mean_marg, sigma = Solution$Xmin, dependent = k, 
                           given = Indices1[!Indices1 %in% k], X.given = rep(0,d), 
                           check.sigma = FALSE)$condVar
        if (var_cond < 0)
          var_cond = 0
        
        for (l in 1:N)
        {
          mean_cond[l] = condMVN(mean = mean_marg, sigma = Solution$Xmin, dependent = k, 
                                 given = Indices1[!Indices1 %in% k], X.given = X_aug[l,Indices1], 
                                 check.sigma = FALSE)$condMean
        }
      }
      
      if (bound == "max")
      {
        var_cond = condMVN(mean = mean_marg, sigma = Solution$Xmax, dependent = k, 
                           given = Indices1[!Indices1 %in% k], X.given = rep(0,d), 
                           check.sigma = FALSE)$condVar
        if (var_cond < 0)
          var_cond = 0
        
        for (l in 1:N)
        {
          mean_cond[l] = condMVN(mean = mean_marg, sigma = Solution$Xmax, dependent = k, 
                                 given = Indices1[!Indices1 %in% k], X.given = X_aug[l,Indices1], 
                                 check.sigma = FALSE)$condMean
        }
      }
      
      
      u_k = rnorm(N, mean = mean_cond, sd = sqrt(var_cond))
      
      X[,k][order(u_k)] = sort(X[,k])
    }  
    
    Y=sort(rowSums(X))
    
    ES =sum(Y[(floor(N*alpha)+1):N])/N/(1-alpha)
    delta=abs(ES-EStemp) 
  }  
  
  return(list(X = X, ES = ES, EStemp = EStemp))
}