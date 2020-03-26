# cor_rand ----------------------------------------------
#' Compute p-values for correlation between time series
#'
#' @param ts1 single time series
#' @param ts2 single time series
#' @param metrand define the method to create random dataset (phase for keeping phase structure, random data, or random + trend (default = TRUE)
#' @param nrep number of repetition
#' @param ... elements sent to function cor()
#' @return p-value of correlation coefficient
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
#' # test the correlation between cod and herring biomass
#' # using surrogates data conserving phase structure
#' cor_rand(baltic$Cod, baltic$Herring)
#' 
#' #Compare to a traditional test
#' cor.test(baltic$Cod, baltic$Herring)
#' @export
#'
cor_rand <- function(ts1, ts2, nrep=100, metrand="phase", ...){
  scor <- cor(ts1, ts2, use="pairwise.complete.obs", ...)
  ncor <- c()
  for (i in 1:nrep){
    ncor <- c(ncor, cor(singleRand(ts1, metrand=metrand), 
                        singleRand(ts2, metrand=metrand), 
                        use="pairwise.complete.obs", ...))
  }
  pval <- 1-sum(scor>ncor)/length(ncor)
  res <- list("cor"=scor, "pvalue"=pval)
  return(res)
}


# pc_test ----------------------------------------------
#' Compute p-values for principal component eigen values
#'
#' @param dat matrix with 
#' @param nrep number of randomly generated dataset
#' @param metrand define the method to create random dataset (either "phase", "rand", or "trend")
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @return p-values of eigen values
#' @details
#' Variable \code{met} can be "phase", "rand" or "trend": 
#' \itemize{
#' \item \code{phase} keep similar phase structure in the dataset
#' \item \code{rand} randomly generated data with similar mean and sd
#' \item \code{trend} randomly generated data with similar trend
#' }
#' @seealso \code{\link{createRand}}
#' @examples 
#' data(baltic)
#' # p-values of eigen values using randomly generated dataset
#' pc_test(baltic)
#' @export
#'
pc_test <- function(dat, nrep=100, metrand="rand", sca = TRUE, logt = FALSE){
  mvar <- ita(dat, npc=3, met="pca", sca=sca, logt=logt)
  randeig <- ita.nrand(dat, npc=3, met="pca", sca=sca, logt=logt, metrand=metrand, nrep=nrep)
  mateig <- matrix(rep(mvar$eig/sum(mvar$eig) * 100, each=nrep), nrow=nrep)
  pval <- 1-apply(mateig>randeig, 2, sum)/nrep
  res <- list("eig"=mvar$eig, "pvalue"=pval)
  return(res)
}
