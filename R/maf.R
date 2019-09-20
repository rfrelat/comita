#===============================================================================
# Min/Max Autocorrelation Factors (MAF) in time (or 1D)
# EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez, J. Rivoirard, Mines-ParisTech
# Last update : 01 march 2008
#
# Description:
# The MAF are a multivariate statistical method, that allows the set of
# initial variables to be decomposed into factors. The autocorrelation of
# which decreases, or the variogram increases, when going from the first
# factors to the last ones. Hence the very first factors extract the part of
# the variables which is the most continuous in time.
#
# Details:
# The present routine makes twice use of PCA. Here 'prcomp' is used for PCA
# instead of 'princomp', as it allows less samples (e.g. years) than variables.
# After original variables are turned into normalized PC, MAF are obtained as
# the PC of their increments; their order is then reversed to begin with the
# smallest variograms (i.e. the half variances of increments).
#
# Note:
# The number of MAF cannot exceed the number of variables, nor the number of
# year increments (number of sampled years - 1). If the number of variables
# tends to be larger than the number of sampled years, the MAF no.i (i=1, ...)
# tends to have period (nb years - 1)*2/i. In particular there will be
# evidence of a high continuity with period (nb years -1)*2 for MAF1,
# (nb years -1) for MAF2, etc, whatever the data, which may not be significant
# out of this series e.g. for additional years. The number of variables should
# then be better reduced well below the number of sampled years.
#
# Inputs:        
#   formula   A formula with no response variable, referring only 
#             to numeric variables.
#   data      A data frame containing the variables in the formula.
#

#maf ---------------------------------------
#' Min/Max Autocorrelation Factors (MAF) in time from  Min/Max Autocorrelation Factors Analysis (MAFA)
#' 
#' The MAF are a multivariate statistical method, that allows the set of
#' initial variables to be decomposed into factors. The autocorrelation of
#' which decreases, or the variogram increases, when going from the first
#' factors to the last ones. Hence the very first factors extract the part of
#' the variables which is the most continuous in time.
#' 
#' @param formula A formula with no response variable, referring only to numeric variables.
#' @param data A data frame containing the variables in the formula.
#' @return The present routine makes twice use of PCA. Here 'prcomp' is used for PCA
#'  instead of 'princomp', as it allows less samples (e.g. years) than variables.
#'  After original variables are turned into normalized PC, MAF are obtained as
#'  the PC of their increments; their order is then reversed to begin with the
#'  smallest variograms (i.e. the half variances of increments).
#'  The output of \code{maf} is  structured similarly than the output from prcomp. 
#' @author M.Woillez and J. Rivoirard
#' @keywords Min/Max Autocorrelation Factors Analysis
#' @note The number of MAF cannot exceed the number of variables, nor the number of
#' year increments (number of sampled years - 1). If the number of variables
#' tends to be larger than the number of sampled years, the MAF no.i (i=1, ...)
#' tends to have period (nb years - 1)*2/i. In particular there will be
#' evidence of a high continuity with period (nb years -1)*2 for MAF1,
#' (nb years -1) for MAF2, etc, whatever the data, which may not be significant
#' out of this series e.g. for additional years. The number of variables should
#' then be better reduced well below the number of sampled years.
#' @references 
#' Solow AR (1994) Detecting Change in the Composition of a Multispecies Community. Biometrics 50:556
#' Woillez M, Rivoirard J, Petitgas P (2009) Using min/max autocorrelation factors of survey-based indicators to follow the evolution of fish stocks in time. Aquat Living Resour 22:193–200
#' @seealso \code{\link{ita.mafa}}
#' @export
#' 
maf <- function(formula,data){
  if(length(data[,1])==1) stop("Cannot be run with only one sample.")
  
  # PCA of original variables
  bidpr<-stats::prcomp(formula,data,scale=TRUE)
  bidn<-sum(bidpr$sdev>1e-10)
  if(bidn<2) stop("stop: strictly less than two PCs")
  bidpr$sdev<-bidpr$sdev[1:bidn]
  bidpr$rotation<-bidpr$rotation[,1:bidn]
  bidpr$x<-bidpr$x[,1:bidn]
  
  # normalization of PCs
  bidim<-dim(bidpr$x)
  bidprsdev<-matrix(rep(bidpr$sdev,each=bidim[1]),ncol=bidim[2])
  bidprx<-bidpr$x/bidprsdev
  bidim<-dim(bidpr$rotation)
  bidprsdev<-matrix(rep(bidpr$sdev,each=bidim[1]),ncol=bidim[2])
  bidpr$rotation<-bidpr$rotation/bidprsdev
  
  # PCA of increments of normalized PCs
  biddtab<- diff(bidprx)
  biddpr<-prcomp(biddtab,center=F,scale=F)
  biddnmaf<-sum(biddpr$sdev>1e-10)
  if(biddnmaf<2) stop("stop: strictly less than two MAFs")
  biddpr$sdev<-biddpr$sdev[1:biddnmaf]
  biddpr$rotation<-biddpr$rotation[,1:biddnmaf]
  biddpr$x<-biddpr$x[,1:biddnmaf]
  biddx<-biddpr$x
  
  # store complete transformation as a PCA output from original data
  maf<-bidpr
  maf$sdev <- biddpr$sdev
  biddim<-dim(biddpr$x)
  maf$rotation<-bidpr$rotation%*%biddpr$rotation
  maf$x<-bidprx%*%biddpr$rotation
  
  # reverse the order of MAFs
  # in order to begin with lowest eigenvalues i.e. lowest variograms
  
  maf$sdev<-.5*(maf$sdev^2)*(biddim[1]-1)/biddim[1]
  names(maf$sdev)<-rev(names(maf$sdev))
  maf$sdev<-rev(maf$sdev)
  lng<-length(dimnames(maf$rotation)[[2]])
  maf$rotation <- maf$rotation[,seq(lng,1)]
  dimnames(maf$rotation)[[2]]<-paste("MAF",1:lng,sep="")
  lng<-length(dimnames(maf$x)[[2]])
  maf$x<-maf$x[,seq(lng,1)]
  dimnames(maf$x)[[2]]<-paste("MAF",1:lng,sep="")
  maf
}
