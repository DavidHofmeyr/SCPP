##    This code forms part of the R package SCPP which implements
##    projection pursuit for spectral clustering
##    Copyright (C) 2017,  David P. Hofmeyr
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.


cluster_performance = function(assigned, labels, beta = 1){
  n <- length(labels)
  T <- table(assigned, labels)
  RS <- rowSums(T)
  CS <- colSums(T)

  ## V-measure

  CK <- - sum(apply(T, 1, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
  KC <- - sum(apply(T, 2, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
  K <- - sum(apply(T, 1, function(x) return(sum(x)*log(sum(x)/n))))/n
  C <- - sum(apply(T, 2, function(x) return(sum(x)*log(sum(x)/n))))/n
  if(C!=0){
    h <- 1 - CK/C
  }
  else{
    h <- 0
  }
  if(K!=0){
    c <- 1 - KC/K
  }
  else{
    c <- 0
  }
  if(h==0 && c==0) v.measure <- 0
  else v.measure <- (1+beta)*h*c/(beta*h+c)

  ## Purity

  purity <- sum(apply(T, 1, function(x) return(max(x))))/n

  ## Adjusted Rand Index

  O <- sum(sapply(T, function(t) choose(t, 2)))
  E <- (sum(sapply(RS, function(t) choose(t, 2)))*sum(sapply(CS, function(t) choose(t, 2))))/choose(n, 2)
  M <- (sum(sapply(RS, function(t) choose(t, 2))) + sum(sapply(CS, function(t) choose(t, 2))))/2
  adj.rand <- (O-E)/(M-E)


  ## Normalised Mutual Information

  prod <- RS%*%t(CS)
  Tp <- T
  Tp[which(T==0)] <- 1e-10
  IXY <- sum(T*log(Tp*n/prod))
  HX <- sum(RS*log(RS/n))
  HY <- sum(CS*log(CS/n))
  NMI <- IXY/sqrt(HX*HY)

  c(adj.rand = adj.rand, purity = purity, v.measure = v.measure, nmi = NMI)
}


norm_vec <- function(v) sqrt(sum(v^2))

SCPP_cluster <- function(X, K, v0 = NULL, ndim = NULL, nMicro = NULL, betamax = NULL, betamin = NULL, smult = NULL, minsize = NULL, minprop = NULL, omega = NULL, type = NULL){

  params <- list()

  if(is.null(type)) type <- 'normalised'
  else if(!type%in%c('standard', 'normalised')) stop('type must be either "standard" or "normalised"')

  if(is.null(ndim)) params$ndim <- 2
  else if(is.numeric(ndim) && ndim%%1==0) params$ndim <- ndim
  else stop('ndim must be an integer')

  if(is.null(nMicro)) params$nMicro <- 200
  else if(is.numeric(nMicro) && nMicro%%1==0) params$nMicro <- nMicro
  else stop('nMicro must be an integer')

  if(is.null(betamax)){
    if(type=='standard') params$betamax <- 1.5
    else params$betamax <- 3
  }
  else if(is.numeric(betamax)) params$betamax <- betamax
  else stop('betamax must be numeric')

  if(is.null(betamin)) params$betamin <- 0.5
  else if(is.numeric(betamin)) params$betamin <- betamin
  else stop('betamin must be numeric')

  params$betaskip <- 0.2

  params$kernel <- 'Gaussian'

  if(is.null(smult)) params$smult <- 1
  else if(is.numeric(smult) && length(smult)==1) params$smult <- smult
  else stop('smult must be numeric')

  params$del <- 1e-2

  params$eps <- 1e-5

  if(params$ndim>1){
    if(is.null(omega)) params$om <- 1
    else if(is.numeric(omega)) params$om <- omega
    else stop('omega must be numeric. Select a positive value for approximately orthogonal projections
              and a negative value for correlated ones')
  }

  if(is.null(minprop)){
    if(is.null(minsize)) params$nmin <- nrow(X)/K/5
    else params$nmin <- minsize
  }
  else params$minprop <- minprop

  n <- nrow(X)

  split_indices <- numeric(2*K-1) + Inf

  ixs <- list(1:n)

  split_indices[1] <- 1/n

  tree <- matrix(0, (2*K-1), 2)
  tree[1,] <- c(1, 1)

  c.split <- spectral_split(X, v0, params, split.index, type)

  pass <- list(c.split$cluster)

  pars <- list(c.split$params)

  vs <- list(c.split$v)


  while(length(ixs)<(2*K-1)){
    id <- which.min(split_indices)

    split_indices[id] <- Inf

    n.clust <- length(ixs)

    ixs[[n.clust+1]] <- ixs[[id]][pass[[id]]]

    ixs[[n.clust+2]] <- ixs[[id]][-pass[[id]]]

    c.split <- spectral_split(X[ixs[[n.clust+1]],], v0, params, split.index, type)

    pass[[n.clust+1]] <- c.split$cluster

    vs[[n.clust+1]] <- c.split$v

    pars[[n.clust+1]] <- c.split$params

    tree[n.clust+1,] <- c(tree[id,1] + 1, 2*tree[id,2]-1)

    c.split <- spectral_split(X[ixs[[n.clust+2]],], v0, params, split.index, type)

    pass[[n.clust+2]] <- c.split$cluster

    vs[[n.clust+2]] <- c.split$v

    pars[[n.clust+2]] <- c.split$params

    tree[n.clust+2,] <- c(tree[id,1] + 1, 2*tree[id,2])


    split_indices[n.clust+1] <- 1/length(ixs[[n.clust+1]])

    split_indices[n.clust+2] <- 1/length(ixs[[n.clust+2]])
  }

  a <- numeric(n)
  for(i in 1:(K-1)) a[ixs[[2*i]]] <- i

  loci <- tree
  for(i in 1:max(tree[,1])){
    rows <- which(tree[,1]==i)
    loci[rows,2] <- rank(tree[rows,2])
  }

  Nodes <- list()
  for(i in 1:length(ixs)) Nodes[[i]] <- list(ixs = ixs[[i]], split = pass[[i]], v = vs[[i]], params = pars[[i]], node = tree[i,], location = loci[i,])

  list(cluster = a, model = tree, Nodes = Nodes, data = X, args = list(v0 = v0, ndim = ndim, nMicro = nMicro, betamax = betamax, betamin = betamin, kernel = kernel, smult = smult, minsize = minsize, minprop = minprop, omega = omega, type = type))


}




f_spectral <- function(Th, X, P, type){

  dm <- ncol(X)
  v <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim){
    CP <- cumprod(c(1, sin(Th[((i-1)*(dm-1)+1):(i*(dm-1))])))
    COS <- c(cos(Th[((i-1)*(dm-1)+1):(i*(dm-1))]), 1)
    v[,i] <- COS*CP
  }

  #### compute projection matrix and project the data
  n <- nrow(X)
  x <- X%*%v

  #### compute constraint set and transformation of the projected data
  nn <- sum(P$x_size)
  sigs <- sapply(1:P$ndim, function(i) sqrt(v[,i]%*%P$COV%*%v[,i]))
  p <- x
  for(i in 1:P$ndim){
    ixlo <- which(x[,i]<(-P$beta*sigs[i]))
    ixhi <- which(x[,i]>(P$beta*sigs[i]))
    p[ixlo,i] <- -P$beta*sigs[i] - P$del*(-P$beta*sigs[i] - x[ixlo,i] + (P$del*(1-P$del))^(1/P$del))^(1-P$del) + P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
    p[ixhi,i] <- P$beta*sigs[i] + P$del*(x[ixhi,i]-P$beta*sigs[i]+(P$del*(1-P$del))^(1/P$del))^(1-P$del)-P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
  }

  #### compute distances between transformed data and associated affinity matrix
  ds <- as.matrix(dist(p))
  w <- exp(-ds^2/2/P$kpar$s^2)

  #### compute Laplacian matrix (L) and first eigenvector (u1)

  if(type=='standard'){
    sqs <- sqrt(P$x_size)
    L <- (diag(rowSums(t(t(w)*P$x_size)))-sqs%*%t(sqs)*w)/nn
    u1 <- sqs
  }
  else if(type=='normalised'){
    w <- P$x_size%*%t(P$x_size)*w
    d <- rowSums(w)
    sqd <- sqrt(d)
    L <- diag(n) - (1/sqd)%*%t(1/sqd)*w
    u1 <- sqd
  }
  else stop('only types "standard" and "normalised" are acceptable for spectral projection pursuit')

  #### compute second eigenvalue of the Laplacian

  f <- rARPACK::eigs_sym(L + 100*u1%*%t(u1), 1, sigma = 1e-10)$values[1]

  #### if more than one column in the projection matrix, add the orthogonality/correlation term
  if(P$ndim>1){
    for(i in 1:(P$ndim-1)){
      for(j in (i+1):P$ndim){
        f <- f + (v[,i]%*%v[,j])^2*P$om
      }
    }
  }
  f
}


df_spectral <- function(Th, X, P, type){

  TH <- matrix(Th, ncol = P$ndim)

  dm <- ncol(X)
  v <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim){
    CP <- cumprod(c(1, sin(TH[,i])))
    COS <- c(cos(TH[,i]), 1)
    v[,i] <- COS*CP
  }

  #### compute projection matrix and project the data
  n <- nrow(X)
  x <- X%*%v

  #### compute constraint set and transformation of the projected data
  nn <- sum(P$x_size)
  sigs <- sapply(1:P$ndim, function(i) sqrt(v[,i]%*%P$COV%*%v[,i]))
  dsigs <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim) dsigs[,i] <- 1/sigs[i]*P$COV%*%v[,i]
  p <- x
  for(i in 1:P$ndim){
    ixlo <- which(x[,i]<(-P$beta*sigs[i]))
    ixhi <- which(x[,i]>(P$beta*sigs[i]))
    p[ixlo,i] <- -P$beta*sigs[i] - P$del*(-P$beta*sigs[i] - x[ixlo,i] + (P$del*(1-P$del))^(1/P$del))^(1-P$del) + P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
    p[ixhi,i] <- P$beta*sigs[i] + P$del*(x[ixhi,i]-P$beta*sigs[i]+(P$del*(1-P$del))^(1/P$del))^(1-P$del)-P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
  }

  #### compute distances between transformed data and associated affinity matrix
  ds <- as.matrix(dist(p))
  w <- exp(-ds^2/2/P$kpar$s^2)

  #### compute Laplacian matrix (L) and first eigenvector (u1)

  if(type=='standard'){
    sqs <- sqrt(P$x_size)
    L <- (diag(rowSums(t(t(w)*P$x_size)))-sqs%*%t(sqs)*w)/nn
    u1 <- sqs
  }
  else if(type=='normalised'){
    w <- P$x_size%*%t(P$x_size)*w
    d <- rowSums(w)
    sqd <- sqrt(d)
    L <- diag(n) - (1/sqd)%*%t(1/sqd)*w
    u1 <- sqd
  }
  else stop('only types "standard" and "normalised" are acceptable for spectral projection pursuit')

  e <- rARPACK::eigs_sym(L + 100*u1%*%t(u1), 1, sigma = 1e-10)

  c.size <- sum(P$x_size[which(e$vectors[,1]<0)])
  if(min(c.size, nn-c.size)<P$nmin) return(numeric(length(Th)))

  #### compute the gradient for each column of the projection matrix
  D <- matrix(0, dm, P$ndim)

  if(type == 'standard'){
    u <- e$vectors[,1]/sqrt(P$x_size)
    distu <- as.matrix(dist(u))
    sizemat <- P$x_size%*%t(P$x_size)
    coeffs <- distu^2*sizemat*w/P$kpar$s^2/nn
  }
  else{
    u <- e$vectors[,1]/sqrt(d)
    lam <- e$values[1]
    distu <- as.matrix(dist(u))
    sumu <- u^2 + matrix(rep(u^2, n), n, n, byrow = TRUE)
    coeffs <- (distu^2-lam*sumu)*w/P$kpar$s^2
  }

  #### compute the constant factors for the sums in derivative formulation

  for(j in 1:P$ndim){

    ixlo <- which(x[,j]<(-P$beta*sigs[j]))
    ixhi <- which(x[,j]>(P$beta*sigs[j]))

    #### DT is the derivative of the transformed projected data w.r.t. the normalised column of the projection matrix
    DT <- X
    DT[ixlo,] <- t(-P$beta*dsigs[,j] + t((P$del*(1-P$del)/(-P$beta*sigs[j]-x[ixlo,j]+(P$del*(1-P$del))^(1/P$del))^P$del)*t(P$beta*dsigs[,j]+t(DT[ixlo,]))))
    DT[ixhi,] <- t(P$beta*dsigs[,j] + t((P$del*(1-P$del)/(-P$beta*sigs[j]+x[ixhi,j]+(P$del*(1-P$del))^(1/P$del))^P$del)*t(t(DT[ixhi,])-P$beta*dsigs[,j])))

    #### dlp is the derivative of the eigenvP$alue w.r.t. corresponding dimension of the transformed projected data
    dlp <- sapply(1:n, function(k){
      sum(coeffs[k,]*(p[,j]-p[k,j]))
    })

    D[,j] <- dlp%*%DT
  }

  #### if more than one column in the projection matrix, add the derivative of the orthogonP$ality/correlation term
  if(P$ndim>1){
    for(i in 1:(P$ndim-1)){
      for(j in (i+1):P$ndim){
        ortho <- (v[,i]%*%v[,j])[1]
        D[,i] <- D[,i] + P$om*2*ortho*v[,j]
        D[,j] <- D[,j] + P$om*2*ortho*v[,i]
      }
    }
  }
  Dfinal <- matrix(0, dm-1, P$ndim)
  for(j in 1:P$ndim){
    Dth <- matrix(0, dm, dm-1)
    Dth[1,1] <- - sin(TH[1,j])
    Dth[2,1] <- cos(TH[1,j])*cos(TH[2,j])
    Dth[2,2] <- -sin(TH[1,j])*sin(TH[2,j])
    if(dm>2){
      COS <- cos(TH[,j])
      SIN <- sin(TH[,j])
      CP <- prod(SIN)
      Dth[dm,] <- COS*CP/SIN
      if(dm>3){
        for(i in 3:(dm-1)){
          COS <- cos(TH[1:i,j])
          SIN <- sin(TH[1:(i-1),j])
          CP <- prod(SIN)
          Dth[i,1:(i-1)] <- COS[i]*COS[1:(i-1)]*CP/SIN
          Dth[i,i] <- -sin(TH[i,j])*CP
        }
      }
    }
    Dfinal[,j] <- D[,j]%*%Dth
  }
  c(Dfinal)
}



spectral_min_check <- function(Th, X, P, type){

  TH <- matrix(Th, ncol = P$ndim)

  dm <- ncol(X)
  v <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim){
    CP <- cumprod(c(1, sin(TH[,i])))
    COS <- c(cos(TH[,i]), 1)
    v[,i] <- COS*CP
  }

  #### compute projection matrix and project the data
  n <- nrow(X)
  x <- X%*%v

  #### compute constraint set and transformation of the projected data
  nn <- sum(P$x_size)
  sigs <- sapply(1:P$ndim, function(i) sqrt(v[,i]%*%P$COV%*%v[,i]))
  dsigs <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim) dsigs[,i] <- 1/sigs[i]*P$COV%*%v[,i]
  p <- x
  for(i in 1:P$ndim){
    ixlo <- which(x[,i]<(-P$beta*sigs[i]))
    ixhi <- which(x[,i]>(P$beta*sigs[i]))
    p[ixlo,i] <- -P$beta*sigs[i] - P$del*(-P$beta*sigs[i] - x[ixlo,i] + (P$del*(1-P$del))^(1/P$del))^(1-P$del) + P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
    p[ixhi,i] <- P$beta*sigs[i] + P$del*(x[ixhi,i]-P$beta*sigs[i]+(P$del*(1-P$del))^(1/P$del))^(1-P$del)-P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
  }

  #### compute distances between transformed data and associated affinity matrix
  ds <- as.matrix(dist(p))
  w <- exp(-ds^2/2/P$kpar$s^2)

  #### compute Laplacian matrix (L) and first eigenvector (u1)

  if(type=='standard'){
    sqs <- sqrt(P$x_size)
    L <- (diag(rowSums(t(t(w)*P$x_size)))-sqs%*%t(sqs)*w)/nn
    u1 <- sqs
  }
  else if(type=='normalised'){
    w <- P$x_size%*%t(P$x_size)*w
    d <- rowSums(w)
    sqd <- sqrt(d)
    L <- diag(n) - (1/sqd)%*%t(1/sqd)*w
    u1 <- sqd
  }
  else stop('only types "standard" and "normalised" are acceptable for spectral projection pursuit')

  e <- rARPACK::eigs_sym(L + 100*u1%*%t(u1), 2, sigma = 1e-10)

  if(abs(e$values[1]-e$values[2])>(e$values[1]*1e-5)) return(TRUE)
  else{
    L <- L+100*(u1%*%t(u1) + e$vectors[,1]%*%t(e$vectors[,1]) + e$vectors[,2]%*%t(e$vectors[,2]))
    repeat{
      enext <- rARPACK::eigs_sym(L, 1, sigma = 1e-10)
      if(abs(enext$values[1]-e$values[1])>(e$values[1]*1e-5)) break
      else{
        e$vectors = cbind(e$vectors, enext$vectors)
        L <- L + 100*enext$vectors[,1]%*%t(enext$vectors[,1])
      }
    }
    n.eigen <- ncol(e$vectors)
    for(id1 in 1:n.eigen){
      for(id2 in id1:n.eigen){
        u <- e$vectors[,id1]
        v <- e$vectors[,id2]
        D <- matrix(0, dm, P$ndim)

        if(type == 'standard'){
          u <- u/sqrt(P$x_size)
          v <- v/sqrt(P$x_size)
          distu <- u - matrix(rep(u, n), n, n, byrow = TRUE)
          distv <- v - matrix(rep(v, n), n, n, byrow = TRUE)
          sizemat <- P$x_size%*%t(P$x_size)
          coeffs <- distu*distv*sizemat*w/P$kpar$s^2/nn
        }
        else{
          u <- u/sqrt(d)
          v <- v/sqrt(d)
          lam <- e$values[1]
          distu <- u - matrix(rep(u, n), n, n, byrow = TRUE)
          distv <- v - matrix(rep(v, n), n, n, byrow = TRUE)
          sumuv <- u*v + matrix(rep(u*v, n), n, n, byrow = TRUE)
          coeffs <- (distu^2-lam*sumu)*w/P$kpar$s^2
        }

        #### compute the constant factors for the sums in derivative formulation

        for(j in 1:P$ndim){

          ixlo <- which(x[,j]<(-P$beta*sigs[j]))
          ixhi <- which(x[,j]>(P$beta*sigs[j]))

          #### DT is the derivative of the transformed projected data w.r.t. the normalised column of the projection matrix
          DT <- X
          DT[ixlo,] <- t(-P$beta*dsigs[,j] + t((P$del*(1-P$del)/(-P$beta*sigs[j]-x[ixlo,j]+(P$del*(1-P$del))^(1/P$del))^P$del)*t(P$beta*dsigs[,j]+t(DT[ixlo,]))))
          DT[ixhi,] <- t(P$beta*dsigs[,j] + t((P$del*(1-P$del)/(-P$beta*sigs[j]+x[ixhi,j]+(P$del*(1-P$del))^(1/P$del))^P$del)*t(t(DT[ixhi,])-P$beta*dsigs[,j])))

          #### dlp is the derivative of the eigenvP$alue w.r.t. corresponding dimension of the transformed projected data
          dlp <- sapply(1:n, function(k){
            sum(coeffs[k,]*(p[,j]-p[k,j]))
          })

          D[,j] <- dlp%*%DT
        }

        #### if more than one column in the projection matrix, add the derivative of the orthogonP$ality/correlation term
        if(P$ndim>1){
          for(i in 1:(P$ndim-1)){
            for(j in (i+1):P$ndim){
              ortho <- (v[,i]%*%v[,j])[1]
              D[,i] <- D[,i] + P$om*2*ortho*v[,j]
              D[,j] <- D[,j] + P$om*2*ortho*v[,i]
            }
          }
        }
        Dfinal <- matrix(0, dm-1, P$ndim)
        for(j in 1:P$ndim){
          Dth <- matrix(0, dm, dm-1)
          Dth[1,1] <- - sin(TH[1,j])
          Dth[2,1] <- cos(TH[1,j])*cos(TH[2,j])
          Dth[2,2] <- -sin(TH[1,j])*sin(TH[2,j])
          if(dm>2){
            COS <- cos(TH[,j])
            SIN <- sin(TH[,j])
            CP <- prod(SIN)
            Dth[dm,] <- COS*CP/SIN
            if(dm>3){
              for(i in 3:(dm-1)){
                COS <- cos(TH[1:i,j])
                SIN <- sin(TH[1:(i-1),j])
                CP <- prod(SIN)
                Dth[i,1:(i-1)] <- COS[i]*COS[1:(i-1)]*CP/SIN
                Dth[i,i] <- -sin(TH[i,j])*CP
              }
            }
          }
          Dfinal[,j] <- D[,j]%*%Dth
        }
        if(max(abs(Dfinal))>1e-7) return(FALSE)
      }
    }
    return(TRUE)
  }
}


spectral_assign <- function(Th, X, P, type){

  dm <- ncol(X)
  v <- matrix(0, dm, P$ndim)
  for(i in 1:P$ndim){
    CP <- cumprod(c(1, sin(Th[((i-1)*(dm-1)+1):(i*(dm-1))])))
    COS <- c(cos(Th[((i-1)*(dm-1)+1):(i*(dm-1))]), 1)
    v[,i] <- COS*CP
  }

  #### compute projection matrix and project the data
  n <- nrow(X)
  x <- X%*%v

  #### compute constraint set and transformation of the projected data
  nn <- sum(P$x_size)
  sigs <- sapply(1:P$ndim, function(i) sqrt(v[,i]%*%P$COV%*%v[,i]))
  p <- x
  for(i in 1:P$ndim){
    ixlo <- which(x[,i]<(-P$beta*sigs[i]))
    ixhi <- which(x[,i]>(P$beta*sigs[i]))
    p[ixlo,i] <- -P$beta*sigs[i] - P$del*(-P$beta*sigs[i] - x[ixlo,i] + (P$del*(1-P$del))^(1/P$del))^(1-P$del) + P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
    p[ixhi,i] <- P$beta*sigs[i] + P$del*(x[ixhi,i]-P$beta*sigs[i]+(P$del*(1-P$del))^(1/P$del))^(1-P$del)-P$del*(P$del*(1-P$del))^((1-P$del)/P$del)
  }

  #### compute distances between transformed data and associated affinity matrix
  ds <- as.matrix(dist(p))

  w <- exp(-ds^2/2/P$kpar$s^2)

  #### compute Laplacian matrix (L) and first eigenvector (u1)

  if(type=='standard'){
    sqs <- sqrt(P$x_size)
    L <- (diag(rowSums(t(t(w)*P$x_size)))-sqs%*%t(sqs)*w)/nn
    u1 <- sqs
  }
  else if(type=='normalised'){
    w <- P$x_size%*%t(P$x_size)*w
    d <- rowSums(w)
    sqd <- sqrt(d)
    L <- diag(n) - (1/sqd)%*%t(1/sqd)*w
    u1 <- sqd
  }
  else stop('only types "standard" and "normalised" are acceptable for spectral projection pursuit')

  #### compute second eigenvector of the Laplacian

  f <- rARPACK::eigs_sym(L, 2, sigma = 1e-10)$vectors

  if(type=='standard'){
    k <- kmeansw(f, 2, weight = P$x_size)$cluster
  }
  else if(type=='normalised'){
    for(i in 1:n) f[i,] <- f[i,]/norm_vec(f[i,])
    k <- kmeansw(f, 2, weight = P$x_size)$cluster
  }
  k
}


pp_spectral <- function(Th, Xu, P, type){

  f <- function(Th) f_spectral(Th, Xu, P, type)
  df <- function(Th) df_spectral(Th, Xu, P, type)

  rad <- 1e-1

  while(rad>=1e-6){
    # find optimal projection angle assuming differentiability

    TH <- optim(Th, f, df, method = 'BFGS')$par

    is.minim <- spectral_min_check(TH, Xu, P, type)

    if(is.minim) return(TH)

    DS <- matrix(0, 2*length(TH)+1, length(TH))
    DS[1,] <- df(TH)
    for(i in 1:(2*length(TH))) DS[i+1,] <- df(TH + rad*(runif(length(TH))-.5))
    H <- chull(DS)
    dir <- QP(t(DS[H,]))
    ndir <- norm_vec(dir)
    if(ndir < 1e-7) rad <- rad/10
    else{
      dir <- dir/ndir
      stp <- 1
      fval <- f(TH)
      fnew <- f(TH-dir*stp)
      while(fnew > (fval - stp/10*ndir)){
        stp <- stp/2
        fnew <- f(TH-dir*stp)
        if(stp<1e-6) break
      }
      if(stp < 1e-6) rad <- rad/10
      else{
        TH <- TH - stp*dir
      }
    }
  }

  TH
}


spectral_split <- function(X, v0, P, split.index, type){

  n <- nrow(rbind(c(),X))
  if(n<=2) return(list(Inf, 1:n))

  P$COV <- cov(X)
  if(ncol(X)>3) evals <- rARPACK::eigs_sym(P$COV, min(20, ncol(X)))$values
  else evals <- eigen(P$COV)$values
  intr <- sum(evals>=1)

  P$kpar <- list()

  if(!is.null(P$s.function)) P$kpar$s <- P$s.function(X)
  else if(!is.null(P$smult)) P$kpar$s <- sqrt(mean(evals[1:intr]))*P$smult*(4/3/n)^(1/(4+intr))
  else P$kpar$s <- P$s

  if(!is.null(P$minprop)) P$nmin <- ceiling(P$minprop*n)

  #### check that data are centralised

  mns <- colMeans(X)
  if(max(abs(mns))>1e-7) X <- t(t(X)-mns)


  km <- MicroClust(X, P$nMicro, P$nMicro, 1000)
  Xu <- km$centers
  P$x_size <- km$size

  if(is.null(v0)){
    if(ncol(X)>3) E <- cbind(c(), c(rARPACK::eigs_sym(P$COV, P$ndim)$vectors))
    else E <- cbind(c(), c(eigen(P$COV)$vectors[,1:P$ndim]))
  }
  else if(is.function(v0)) E <- cbind(c(), v0(X))
  else E <- cbind(c(), v0)


  f.vals <- numeric(ncol(E))

  for(i in 1:ncol(E)){

    P$beta <- P$betamax

    #### compute projection and clustering until the desired balance is met
    repeat{

      v_init <- matrix(E[,i], ncol = P$ndim)
      theta <- matrix(0, ncol(X)-1, P$ndim)
      for(j in 1:P$ndim) theta[,j] <- MakeTheta(v_init[,j])

      #### compute optimal projection
      TH <- c(pp_spectral(theta, Xu, P, type))

      vv <- MakeV(TH, P$ndim)
      for(j in 1:P$ndim) vv[,j] <- vv[,j]/norm_vec(vv[,j])
      if(P$ndim==1) vv <- cbind(vv, eigen(cov(X-X%*%vv%*%t(vv)))$vectors[,1])
      assgn <- spectral_assign(TH, Xu, P, type)
      clusters <- numeric(sum(P$x_size))
      for(j in 1:length(P$x_size)) clusters[which(km$cluster==j)] <- assgn[j]


      #### if balance of clusters is not met decrease alpha and try again
      if((min(sum(clusters==1), sum(clusters==2))<P$nmin) && (P$beta>P$betamin)){
        P$beta <- max(P$beta - P$betaskip, P$betamin)
      }
      else break
    }
    f.vals[i] <- f_spectral(TH, Xu, P, type)
    if(f.vals[i]<=min(f.vals[1:i])){
      pass_opt <- which(clusters==1)
      v_opt <- vv
      params_opt <- P
    }
  }
  list(cluster = pass_opt, v = v_opt, params = params_opt)
}

MicroClust <- function(X, k, kthresh, CLthresh){
  if(nrow(X)<=kthresh){
    return(list(centers = X, size = rep(1, nrow(X)), cluster = 1:nrow(X)))
  }
  else if(FALSE){
    h <- hclust(dist(X), method = 'complete')
    mem <- cutree(h, k=k)
    siz <- numeric(k)
    for(i in 1:k) siz[i] <- sum(mem==i)
    cent <- matrix(0, k, ncol(X))
    for(i in 1:k) cent[i,] <- colMeans(matrix(X[which(mem==i),], nrow = siz[i]))
    return(list(centers = cent, size = siz, cluster = mem))
  }
  else{
    mn = colMeans(X)
    ds = pracma::distmat(X, mn)
    C = which.max(ds)
    ds = pracma::distmat(X, X[C,])
    for(i in 2:k){
      C = c(C, which.max(ds))
      dsnew = pracma::distmat(X, X[C[i],])
      ds = apply(cbind(ds, dsnew), 1, min)
    }
    kmeans(X, X[C,])
  }
}

MakeTheta <- function(v){
  Th <- numeric(length(v)-1)
  Th[1] <- acos(v[1])
  for(j in 2:(length(v)-1)){
    Th[j] <- acos(v[j]/(prod(sin(Th[1:(j-1)])+1e-10)))
  }
  Th
}

MakeV <- function(Th, ndim){

  TH <- matrix(Th, ncol = ndim)

  dm <- length(Th)/ndim+1
  v <- matrix(0, dm, ndim)
  for(i in 1:ndim){
    CP <- cumprod(c(1, sin(TH[,i])))
    COS <- c(cos(TH[,i]), 1)
    v[,i] <- COS*CP
  }
  v
}


QP = function(G){
  n = length(G[1,])
  d = numeric(n)
  A = rbind(-rep(1, n), diag(n))
  b = c(-1, numeric(n))

  # adding a diagonal PSD perturbation 1e-5*diag(n) helps with some
  # potential numerical problems without (substantially) affecting
  # the result
  c(G%*%solve.QP(t(G)%*%G + 1e-5*diag(n), d, t(A), b, meq = 1)$solution)
}

kmeansw <- function(X, k, wts){
  mn = colMeans(X)
  ds = pracma::distmat(X, mn)
  C = which.max(ds)
  ds = pracma::distmat(X, X[C,])
  for(i in 2:k){
    C = c(C, which.max(ds))
    dsnew = pracma::distmat(X, X[C[i],])
    ds = apply(cbind(ds, dsnew), 1, min)
  }
  FactoClass::kmeansW(X, unique(X[C,]), weight = wts)  
}
