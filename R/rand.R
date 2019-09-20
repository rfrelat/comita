# createRand -------------------------------------------
#' Create a random dataset with similar characteristics than original data
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @param metrand define the method to create random dataset (either "phase", "rand", or "trend")
#' @return a randomly generated matrix of similar size than \code{dat}
#' @details
#' Variable \code{metrand} can be "phase", "rand" or "trend": 
#' \itemize{
#' \item \code{phase} keep similar phase structure in the dataset
#' \item \code{rand} randomly generated data with similar mean and sd
#' \item \code{trend} randomly generated data with similar trend
#' }
#' @seealso \code{\link{ita.rand}}, \code{\link{ita.nrand}}, \code{\link{singleRand}}
#' @examples 
#' data(baltic)
#' createRand(baltic, metrand="phase")
#' @export
#'
createRand <- function(dat, metrand=c("phase", "rand", "trend")){
  rand <- c()
  for (j in 1:ncol(dat)){
      # get time series j
      rand <- cbind(rand, singleRand(dat[,j], metrand=metrand))
  } 
  return(rand)
}


# singleRand ----------------------------------------------
#' Create a random time series with similar characteristics than original data
#'
#' @param ts single time series
#' @param metrand define the method to create random dataset (phase for keeping phase structure, random data, or random + trend (default = TRUE)
#' @return a randomly generated time series
#' @details
#' Variable \code{metrand} can be "phase", "rand" or "trend": 
#' \itemize{
#' \item \code{phase} keep similar phase structure in the dataset
#' \item \code{rand} randomly generated data with similar mean and sd
#' \item \code{trend} randomly generated data with similar trend
#' }
#' @seealso \code{\link{createRand}}
#' @examples 
#' data(baltic)
#' singleRand(baltic[,3], metrand="phase")
#' @export
#'
singleRand <- function(ts, metrand=c("rand", "phase", "trend")){
  tmp <- as.numeric(ts)
  xseq <- 1:length(tmp)
  if (length(metrand)>1){
    metrand <- metrand[1]
  }
  if (metrand=="rand"){
    randts <- rnorm(length(tmp), mean(tmp), stats::sd(tmp))
  }
  if (metrand=="phase"){
    # remove trend and keep residuals of lm
    lmj <- lm(tmp~xseq)
    resj <- residuals(lmj)
    # compute surrogate of detrended data
    surj <- as.numeric(fractal::surrogate(resj, method="phase"))
    # add trend in the dataset
    randts <- surj + xseq*lmj$coefficients[2] + lmj$coefficients[1]
  } 
  if (metrand=="trend"){
    # remove trend and keep residuals of lm
    lmj <- lm(tmp~xseq)
    resj <- residuals(lmj)
    # create random data
    randj <- rnorm(length(tmp), mean(resj), sd(resj))
    # add trend in the dataset
    randts <- randj + xseq*lmj$coefficients[2] + lmj$coefficients[1]
  }
  return(randts)
}

# ita.rand -------------------------------------------
#' Run Integrated Trend Analysis on randomly generated dataset
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @param npc number of new dimensions
#' @param met selected ITA method
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @param metrand define the method to create random dataset (either "phase", "rand", or "trend")
#' @return percentage variance explained by the successive new dimensions
#' @examples 
#' data(baltic)
#' #Eigen values with random data of similar size than 'baltic'
#' mrand <- ita.rand(baltic)
#' @seealso \code{\link{createRand}}, \code{\link{ita.nrand}}
#' @export
#'
ita.rand <- function(dat, npc=2, met=1, sca=TRUE, logt=FALSE, metrand="phase"){
  rand <- createRand(dat, metrand=metrand)
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
#' Extension of ita.rand to generate multiple datasets
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @param npc number of new dimensions
#' @param met selected ITA method
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @param metrand define the method to create random dataset (either "phase", "rand", or "trend")
#' @param nrep number of randomly generated dataset
#' @return percentage variance explained by the successive new dimensions
#' @examples 
#' data(baltic)
#' #matrice of eigen values with multiple randomly generated dataset 
#' mrand <- ita.nrand(baltic)
#' 
#' # To be compared with eigen values with real observations
#' mvar <- ita(baltic, npc=2)
#' plot_nrand(mvar, mrand)
#' @export
#'
ita.nrand <- function(dat, npc=2, met=1, sca=TRUE, logt=FALSE, metrand="phase", nrep=10){
  eigmat <- c()
  for (i in 1:nrep){
    eigmat <- rbind(eigmat, ita.rand(dat, npc=npc, met=met, sca=sca, logt=logt, metrand=metrand))
  }
  return(eigmat)
}
