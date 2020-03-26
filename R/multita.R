#multita --------------------------------------------
#' Run multiple Integrated Trend Analysis
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @param met selected ITA methods (can be a vector of multiple methods)
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @details
#' Variable \code{met} are either string or number indicating: 
#' \itemize{
#' \item \code{ita.pca} eeee
#' \item \code{ita.pcaV} scores of the variables on the new dimensions (if any)
#' \item \code{ita.} eigen values (for PCA associated methods)
#' \item \code{ita.} dataset used for ITA
#' \item \code{ita.} number of new dimensions
#' \item \code{ita.} method used for integrated trend analysis 
#' }
#' @return A list of output from \code{ita} with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} scores of the variables on the new dimensions (if any)
#' \item \code{eig} eigen values (for PCA associated methods)
#' \item \code{dat} dataset used for ITA
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' Additionnally, information about computing time is saved in \code{cptime}
#' @keywords multivariate methods
#' @note The computation can be long, depending on the selected methods.
#' @export
#'
multita <- function(dat, npc, met=c("pca","nmds"), sca =TRUE, logt=FALSE){
  multimvar <- list()
  cptime <- c()
  for (i in met){
    strt <- Sys.time()
    newmvar <- ita(dat, npc, met=i, sca = sca, logt=logt)
    cptime <- c(cptime, as.numeric(Sys.time()-strt))
    if (!is.null(newmvar$ts)){
      newmvar <- negPC(newmvar) #negative trend
    }
    multimvar[[length(multimvar)+1]] <- newmvar
    #test correlation of successive PC with PCA, and inverse if needed
    if (!is.null(newmvar$ts)){
      for (i in 1:npc){
        if(stats::cor(multimvar[[1]]$ts[,i], newmvar$ts[,i])<0){
          newmvar$ts[,i] <- -newmvar$ts[,i]
          if (!is.null(newmvar$co)){
            newmvar$co[,i] <- -newmvar$co[,i]
          }
        }
      }
      multimvar[[length(multimvar)]] <- newmvar
      names(multimvar)[length(multimvar)] <- newmvar$mita
    } else {
      print(paste("Warning:",newmvar$mita, "not computed"))
      multimvar[[length(multimvar)]] <- NULL
      cptime <- cptime[-length(cptime)]
    }
  }
  names(cptime) <- names(multimvar)
  print(cptime)
  attr(multimvar, "cptime") <- cptime 
  return(multimvar)
}


#multitaP --------------------------------------------
#' Run multiple Integrated Trend Analysis with parallele computing
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @param met selected ITA methods (can be a vector of multiple methods)
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @param nCores number of cores for parallele processing (default=2)
#' @details
#' Variable \code{met} are either string or number indicating: 
#' \itemize{
#' \item \code{ita.pca} eeee
#' \item \code{ita.pcaV} scores of the variables on the new dimensions (if any)
#' \item \code{ita.} eigen values (for PCA associated methods)
#' \item \code{ita.} dataset used for ITA
#' \item \code{ita.} number of new dimensions
#' \item \code{ita.} method used for integrated trend analysis 
#' }
#' @return A list of output from \code{ita} with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} scores of the variables on the new dimensions (if any)
#' \item \code{eig} eigen values (for PCA associated methods)
#' \item \code{dat} dataset used for ITA
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' Additionnally, information about computing time is saved in \code{cptime}
#' @keywords multivariate methods
#' @note The computation can be long, depending on the selected methods.
#' @export
#'
multitaP <- function(dat, npc, met=c("pca","mds"), sca =TRUE, logt=FALSE, nCores=2){
  multimvar <- list()
  cores <- parallel::detectCores()
  if (nCores>cores) {nCores=cores}
  cl<-parallel::makeCluster(nCores)
  doParallel::registerDoParallel(nCores)
  i <- 1
  multimvar <- foreach::foreach(i=seq_along(met)) %dopar% {
    newmvar <- ita(dat, npc, met=met[i], sca = sca, logt=logt)
    if (!is.null(newmvar$ts)){
      newmvar <- negPC(newmvar) #negative trend
    }
    return(newmvar)
  }
  metname <- sapply(multimvar, FUN = function(dat){dat$mita})
  names(multimvar) <- metname 
  parallel::stopCluster(cl)
  return(multimvar)
}

