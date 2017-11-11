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

SCPP_plot <- function(sol, labels = NULL, node.numbers = TRUE){
  op <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0))
  X <- sol$data
  plot(0, xlim = c(0, 1), ylim = c(0, 1), cex = 0, xaxt = 'n', yaxt = 'n')

  # determine the geometry of the tree for plotting

  d <- max(sol$model[,1])
  ns <- sapply(1:d, function(i) sum(sol$model[,1]==i))
  w <- max(ns)
  width <- 1/w
  height <- .8/d

  # if labels are supplied, ensure they are integers 1:K (K the number of classes)

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  # add the plots of each partition in the hierarchy recursively using add_subtree

  add_subtree(sol, X, 1, 0, 1, 1, height, width*2/3, labels, node.numbers)
  par(op)
}


add_subtree <- function(sol, X, node, L, U, y, h, w, labels, node.numbers){
  is.leaf <- 1-sum((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==sol$model[node,2]*2))
  if(is.leaf){
    Xp <- X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v
    Xp[,1] <- w*.7*(Xp[,1]-min(Xp[,1]))/(max(Xp[,1])-min(Xp[,1])) + (L+U-w*.7)/2
    Xp[,2] <- h*(Xp[,2]-min(Xp[,2]))/(max(Xp[,2])-min(Xp[,2])) + y - h
    if(is.null(labels)){
      if(sol$model[node,2]%%2) points(Xp, cex = .5, col = SCPP_Palette[4])
      else points(Xp, col = SCPP_Palette[2], cex = .5)
    }
    else points(Xp, col = SCPP_Palette[labels[sol$Nodes[[node]]$ixs]], cex = .5)

    if(node.numbers) text((L+U-w*.8)/2, y-h*.9, as.character(node), cex = .7)
  }
  else{
    Xp <- X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v
    Xp[,1] <- w*(Xp[,1]-min(Xp[,1]))/(max(Xp[,1])-min(Xp[,1])) + (L+U-w)/2
    Xp[,2] <- h*(Xp[,2]-min(Xp[,2]))/(max(Xp[,2])-min(Xp[,2])) + y - h
    if(is.null(labels)){
      cols <- numeric(nrow(Xp)) + 4
      cols[sol$Nodes[[node]]$split] <- 2
      points(Xp, col = SCPP_Palette[cols], cex = .5)
    }
    else points(Xp, col = SCPP_Palette[labels[sol$Nodes[[node]]$ixs]], cex = .5)

    if(node.numbers) text((L+U-w*1.1)/2, y-h*.9, as.character(node), cex = .7)

    k1 <- which(((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==(sol$model[node,2]*2-1)))==1)
    k2 <- which(((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==sol$model[node,2]*2))==1)
    w1 <- subtree_width(sol, k1)
    w2 <- subtree_width(sol, k2)
    if(w1>w2) w1 <- w1*1.2
    else if(w1<w2) w2 <- w2*1.2
    M <- L + (U-L)*w1/(w1+w2)
    segments((L+U)/2, y-h, (L+U)/2, y-1.125*h)
    segments((L+M)/2, y-1.125*h, (M+U)/2, y-1.125*h)
    segments((L+M)/2, y-1.125*h, (L+M)/2, y-1.25*h)
    segments((U+M)/2, y-1.125*h, (U+M)/2, y-1.25*h)
    add_subtree(sol, X, k1, L, M, y - 1.25*h, h, (M-L)/w1, labels, node.numbers)
    add_subtree(sol, X, k2, M, U, y - 1.25*h, h, (U-M)/w2, labels, node.numbers)
  }
}


subtree_width <- function(sol, node){
  upper <- lower <- sol$model[node,2]
  depth <- sol$model[node,1]
  rows <- c(node)
  for(i in (depth+1):max(sol$model[,1])){
    upper <- 2*upper
    lower <- 2*lower - 1
    add.rows <- (sol$model[,1]==i)*(sol$model[,2]>=lower)*(sol$model[,2]<=upper)
    rows <- c(rows, which(add.rows==1))
  }
  if(length(rows)==1) return(1)
  sub.model <- sol$model[rows,]
  ns <- sapply(unique(sub.model[,1]), function(i) sum(sub.model[,1]==i))
  max(ns)
}

SCPP_Palette <- c("#0072B2", "#000000", "#808080", "#e79f00", "#009E73", "#F0E442", "#660066", "#9ad0f3", "#D55E00",
             "#CC79A7")

