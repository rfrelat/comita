# createRand -------------------------------------------
#' Create a random dataset with similar characteristics than original data
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @param ar whether to keep the autocorrelation structure (default = TRUE)
#' @param tr whether to add a trend similar to the original one (default =FALSE)
#' @return a randomly generated matrix of similar size than \code{dat}
#' @seealso \code{\link{ita.rand}}, \code{\link{ita.nrand}}
#' @export
#'
createRand <- function(dat, ar=TRUE, tr=FALSE){
  rand <- c()
  ny <- nrow(dat)
  if (ar){
    for (j in 1:ncol(dat)){
      rand <- cbind(rand, as.numeric(surrogate(dat[,j],method="phase")))
    } 
  } else {
    for (j in 1:ncol(dat)){
      rand <- cbind(rand, rnorm(ny, mean(dat[,j]), sd(dat[,j])))
    }
  }
  if (tr){
    rand[,j] <- rand[,j]+1:ny*trend(dat[,j])
  }
  return(rand)
}

# ita.rand -------------------------------------------
#' Run Integrated Trend Analysis on randomly generated dataset
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @param npc number of new dimensions
#' @param met selected ITA method
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @param ar whether to keep the autocorrelation structure in randomly generated data (default = TRUE)
#' @param tr whether to add a trend in randomly generated data similar to the original one (default = FALSE)
#' @return percentage variance explained by the successive new dimensions
#' @seealso \code{\link{createRand}}, \code{\link{ita.nrand}}
#' @export
#'
ita.rand <- function(dat, npc=2, met=1, sca=TRUE, logt=FALSE, ar=TRUE, tr=FALSE){
  rand <- createRand(dat, ar, tr)
  if (logt){
    rand <- log(rand+ifelse(min(rand)<=0, abs(min(rand))+1, 0))
  }
  if (sca){
    rand <- scale(rand) #rescale variables
  }
  rand <- as.matrix(rand)
  
  #multivariate method
  mvar <- ita(rand, npc, met)
  return(mvar$eig/sum(mvar$eig)*100)
}

# ita.nrand -------------------------------------------
#' @describeIn ita.rand extension to generate multiple datasets
#' @param nrep number of randomly generated dataset
#' @export
#'
ita.nrand <- function(dat, npc=2, met=1, sca=TRUE, logt=FALSE, ar=TRUE, tr=FALSE, nrep=10){
  eigmat <- c()
  for (i in 1:nrep){
    eigmat <- rbind(eigmat, ita.rand(dat, npc=npc, met=met, sca=sca, logt=logt, ar=ar, tr=tr))
  }
  #removed arguments: plot=FALSE, col="blue"
  # if (plot){
  #   mvar <- eta(dat, npc=npc, met = met, scale = sca, logt = log, detrend = dif)
  #   plot.nrand(mvar, eigmat, col=col)
  # }
  return(eigmat)
}
