#' @title Construct mfa
#' @description Creates an object of class \code{"mfa"}
#' @param data dataset on which MFA is conducted
#' @param sets list of vectors indicating the sets of variables (i.e. the blocks)
#' @param ncomps integer indicating how many number of components (i.e. factors) are to be extracted.
#' @param center either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @param scale either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @return an object of class mfa
#' @export
#' @examples
#' url="https://raw.githubusercontent.com/ucb-stat243/stat243-fall-2016/master/problem-sets/final-project/data/wines.csv" 
#' data=read.csv(file=url, header = TRUE, sep = ",")
#' datanum <- data[,-1]
#'
#' # a list of numeric vectors with the position of the active variables in the data table.
#' sets <- list(1:6, 7:12, 13:18, 19:23, 24:29, 30:34, 35:38, 39:44, 45:49, 50:53)
#'
#' wine <- mfa_gen(data = datanum, sets = sets, ncomps = 2)

mfa_gen <- function(data, sets, ncomps = NULL, center = TRUE, scale = TRUE){
  if(class(data)!="data.frame" & class(data)!="matrix"){
    stop("Data should be data frame or matrix")
  }

  # store the weights
  alpha <- NULL

  # data.pro stores the normalized grand table
  data.pro <- matrix(nrow=nrow(data), ncol=tail(sets[[length(sets)]],n=1))

  # conduct SVD on each sub table
  for (i in 1:length(sets)){
    #decompose the grand table
    dat <- data[,sets[[i]]]
    #normalize subtable
    dat <- scale(dat,center = T, scale = T)/sqrt(nrow(data)-1)
    data.pro[,sets[[i]]] <- dat

    # SVD
    s <- svd(dat)
    alpha[i] <- 1/(s$d[1]^2)
  }

  # diagonal matrix
  A <- diag(rep(alpha,lapply(sets,length)))

  # GSVD
  M <- diag(1/nrow(data),nrow=nrow(data))
  gsvd <- svd(sqrt(M)%*%data.pro%*%sqrt(A))
  P <- solve(sqrt(M))%*%gsvd$u
  D <- diag(gsvd$d)
  Q <- solve(sqrt(t(A)))%*%gsvd$v

  # Full Factor Scores
  F <- P %*% D

  # requested Factor Scores
  Factorscore = P%*%D[,1:ncomps]

  # Partial Factor Scores
  # create a list of matrix to store the PFS for each experts
  partial.factor <- vector("list", length(sets))
  for (i in 1:length(sets)){
    #decompose the processed data table
    dat <- data.pro[,sets[[i]]]
    partial.factor[[i]] <- length(sets)*alpha[i]*dat%*%Q[sets[[i]],1:ncomps]
  }


  # vector of eignvalues
  eigenvalues <- diag(D)^2

  # create the object mfa
  object <- list("eigenvalues" = eigenvalues, "FactorScore" = Factorscore, "PartialFactorScores" = partial.factor, "Loadings" = Q[,1:ncomps], "alpha" = alpha)
  class(object) <- "mfa"
  return(object)

}

