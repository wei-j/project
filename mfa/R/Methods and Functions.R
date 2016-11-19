# ===============================
# print method
# ===============================
print.mfa <- function(object, subtable = NULL, ...) {
  
  # print eignvalues
  cat('Eigenvalues (listed in the order of components) \n')
  e = matrix(round(object$eigenvalues, 3), nrow = 1, ncol = length(object$eigenvalues))
  # change the row names and column names of the matrix
  rownames(e) = ""
  colnames(e) = 1:length(object$eigenvalues)
  print(e)
  
  cat("\n")
  
  # print partial factor scores if subtable index is indicated
  if (class(subtable) == "numeric") {
    cat(paste('Partial Factor Score Matrix for Assessor No.', subtable, '\n', sep = ""))
    
    # change the row names and column names of the matrix
    rownames(object$`PartialFactorScores`[[subtable]]) = 1:length(object$eigenvalues)
    colnames(object$`PartialFactorScores`[[subtable]]) = 1:ncol(object$Loadings)
    # print the matrix
    print(object$`PartialFactorScores`[[subtable]])
  }
  
  # hide the default description of object
  invisible(object)
}

# ===============================
# plot method
# ===============================

# ===============================
# summaries of eigenvalues
# ===============================
eigen <- function(x) UseMethod("eigen")

eigen.mfa <- function(x){
  eigenvalue <- round(as.numeric(x$eigenvalues),3)
  singular <- round(sqrt(eigenvalue),3)
  cumulative <- cumsum(eigenvalue)
  intertia <- round(eigenvalue/sum(eigenvalue)*100,0)
  cumulative.intertia <- cumsum(intertia)
  eigenvalue.table <- t(data.frame(singular,eigenvalue,cumulative, intertia, cumulative.intertia))
  colnames(eigenvalue.table) <- 1:length(eigenvalue)
  return(eigenvalue.table)
}

# ===============================
# contributions
# ===============================

# contribution of an observation to a dimension

# Define the generic function
contri_obs <- function(x, ...) UseMethod("contri_obs")

contri_obs.mfa <- function(x) {
  # Get the eigenvalues
  lambda <- x$eigenvalues[1:ncol(x$Loadings)]
  # Inverse the eigenvalues
  inv_lambda <- lambda^(-1)
  # Matrix of mass
  M <- diag(1/length(x$eigenvalues),nrow=length(x$eigenvalues))
  # Square the factor scores
  F_squared <- (x$FactorScore)^2
  # Contribution of an observation to a dimension
  ctr_obs <- M %*% F_squared * inv_lambda
  
  return(ctr_obs)
}

#----------------------------------------------
# Contribution of a variable to a dimension

# Define the generic function
contri_var <- function(x, ...) UseMethod("contri_var")

contri_var.mfa <- function(x) {
  # Matrix A
  A <- x$MatrixA
  # Square the loadings
  Q_squared <- (x$Loadings)^2
  # Contribution of a variable to a dimension
  ctr_var <- A %*% Q_squared
  
  return (ctr_var)
}

#----------------------------------------------
# Contribution of a table to a dimension

# Define the generic function
contri_table <- function(x, ...) UseMethod("contri_table")

contri_table.mfa <- function(x) {
  # Matrix A
  A <- x$MatrixA
  # Square the loadings
  Q_squared <- (x$Loadings)^2
  # Contribution of a variable to a dimension
  ctr_var <- A %*% Q_squared
  # Get the sets
  set <- x$sets
  ctr_table <- matrix(nrow = length(set), ncol = ncol(ctr_var))
  for (i in 1:length(set)) {
    ctr_table[i,] <- colSums(ctr_var[set[[i]],])
  }
  return(ctr_table)
}

# ===============================
# Rv coefficient
# ===============================
RV <- function(table1,table2){
  table1 <- as.matrix(table1)
  table2 <- as.matrix(table2)
  numerator <- sum(diag((table1 %*% t(table1)) %*% (table2 %*% t(table2))))
  denominator <- sqrt((sum(diag((table1 %*% t(table1)) %*% (table1 %*% t(table1))))) * (sum(diag((table2 %*% t(table2)) %*% (table2 %*% t(table2)))))) 
  return(numerator/denominator)
}


RV_table <- function(dataset, sets){
  k <- length(sets)
  table <- matrix(nrow= k, ncol= k)
  for (i in 1:k){
    for (j in 1:k){
      table[i,j] <- RV(dataset[,sets[[i]]], dataset[,sets[[j]]])
    }
  }
  return(table)
}

# ===============================
# Lg coefficient
# ===============================

Lg <- function(table1,table2,alpha1,alpha2){
  table1 <- as.matrix(table1)
  table2 <- as.matrix(table2)
  s <- sum(diag((table1 %*% t(table1)) %*% (table2 %*% t(table2))))
  return(s*alpha1*alpha2)
}


Lg_table <- function(dataset, sets){
  k <- length(sets)
  table <- matrix(nrow= k, ncol= k)
  for (i in 1:k){
    for (j in 1:k){
      gam1 <- svd(dataset[,sets[[i]]])$d[1]
      gam2 <- svd(dataset[,sets[[j]]])$d[1]
      a1 <- 1/(gam1^2)
      a2 <- 1/(gam2^2)
      table[i,j] <- Lg(dataset[,sets[[i]]], dataset[,sets[[j]]], alpha1 = a1, alpha2 = a2)
    }
  }
  return(table)
}

# ===============================
# Bootstrap
# ===============================

sets <- list(1:6, 7:12, 13:18, 19:23, 24:29, 30:34, 35:38, 39:44, 45:49, 50:53)
#just bootstrap one a time
bootstrapprep=function(userset,dataset,ncomp){
  bootfactorscore=mfa_gen(dataset,sets=userset,ncomps =ncomp)$PartialFactorScores
  XB=sample(1:length(userset),length(userset),replace=TRUE)
  add=matrix(0,nrow = nrow(dataset),ncol=ncomp)
  for (i in 1:length(XB)){
    Factor=bootPfactorscore[[XB[i]]]
    add=add+Factor 
    Fboot=1/length(XB)*add
  }
  return (Fboot)
}


# bootstrap by L times  
bootstrap=function(L,userset,dataset,ncomp){
  Fl=matrix(0,nrow = nrow(dataset),ncol=ncomp)
  Fstar=list()
  for (k in 1:L){
    Fstar[[k]]=bootstrapprep(userset,dataset,ncomp)
    Fl=Fl+Fstar[[k]]
  }
  #mean
  L=1000
  F.starbar=1/L*Fl
  #std 
  F.starstdsum=matrix(0,nrow = nrow(datanum),ncol=2)
  for (k in 1:L){
    F.starstdsum=F.starstdsum+(F.starbar-Fstar[[k]])^2
  }
  F.starstd=sqrt(1/L*F.starstdsum)
  #Ratio
  T.star=F.starbar * F.starstd^(-1)
  return (T.star)
}