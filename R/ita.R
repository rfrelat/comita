# Define the multiple methods to compute Integrative Trend Analysis
# ita.pca ---------------------------------------
#' Internal function to run a Principal Component analysis
#'
#' @param dat matrix with 
#' @param npc number of selected principal components
#' @return Result of multivariate analysis with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @export
#'
ita.pca<- function(dat, npc){
  #use ade4 package to compute pca #require(ade4)
  mvar <- ade4::dudi.pca(dat, nf=npc, scannf = FALSE, center = FALSE, scale = FALSE)
  pvar <- getpvar(mvar$li, mvar$co, dat, npc)
  res <- list("ts"=mvar$li, "co"=mvar$co, "eig"= mvar$eig, 
              "pvar"=pvar, "dat"=dat, "npc"=npc, "mita"="PCA")
  return(res)
}

# ita.pcaV ---------------------------------------
#' Internal function to run a Principal Component analysis followed by Varimax
#'
#' @param dat matrix with 
#' @param npc number of selected principal components
#' @return Result of multivariate analysis with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @export
#'
ita.pcaV<- function(dat, npc){
  mvar <- ade4::dudi.pca(dat, nf=npc, scannf = FALSE, center = TRUE, scale = TRUE)
  rawLoadings     <- as.matrix(mvar$co[,1:npc]) #%*% diag(sqrt(mvar$eig), npc, npc)
  vari <- stats::varimax(rawLoadings)
  ts <- scale(mvar$li[,1:npc]) %*% vari$rotmat
  co <- unclass(vari$loadings)
  pvar <- getpvar(ts, co, dat, npc) 
  res <- list("ts"=ts, "co"=co, "eig"= mvar$eig, "pvar"=pvar, 
              "dat"=dat, "npc"=npc, "mita"="PCA+Varimax")
  return(res)
}

#ita.dfa ----------------------------------------
#' Run a Dynamic Factor Analysis with diagonal and unequal matrix model
#'
#' @param dat matrix with 
#' @param npc number of processes
#' @return Result of multivariate analysis with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note R code inspired from https://nwfsc-timeseries.github.io/atsa-labs/sec-dfa.html
#' Holmes, E. E., M. D. Scheuerell, and E. J. Ward. Applied time series analysis for fisheries and environmental data. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112.
#' @references Zuur AF, Fryer RJ, Jolliffe IT, Dekker R, Beukema JJ (2003) Estimating common trends in multivariate time series using dynamic factor analysis. Environmetrics 14:665–685
#' 
ita.dfa <- function(dat, npc){
  #parameters 
  mm <- npc ## number of processes
  BB <- "identity"  # diag(mm) ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  uu <- "zero"  # matrix(0,mm,1)  ## 'uu' is a column vector of 0's
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1) 
  cc <- "zero"  # matrix(0,1,wk_last)
  QQ <- "identity"  # diag(mm) ## 'QQ' is identity
  aa <- "zero" ## 'aa' is the offset/scaling
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  RR <- "diagonal and unequal"## 'RR' is var-cov matrix for obs errors
  ## 'ZZ' is loadings matrix
  #tab <- tab[,6:10]
  Z_vals <- c()
  for (i in 1:ncol(dat)){
    for (j in 1:npc){
      Z_vals <- c(Z_vals, ifelse(j<=i, paste0("z",i, j), 0))
    }
  }
  ZZ <- matrix(Z_vals, nrow=ncol(dat), ncol=npc, byrow=TRUE)
  ## list with specifications for model vectors/matrices
  mod_list <- list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZ, A=aa, D=DD, d=dd, R=RR)
  ## list with model inits
  init_list <- list(x0=matrix(rep(0,mm),mm,1))
  ## list with model control parameters
  con_list <- list(maxit=3000, allow.degen=TRUE)
  
  ## fit MARSS
  ret <- as.matrix(t(dat))
  mvar <- MARSS::MARSS(y=ret, model=mod_list, inits=init_list, control=con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(mvar, type="matrix")$Z
  ## get the inverse of the rotation matrix
  H_inv <- stats::varimax(Z_est)$rotmat
  ## rotate factor loadings
  Z_rot = Z_est %*% H_inv   
  ## rotate processes
  proc_rot = solve(H_inv) %*% mvar$states
  
  ts=t(proc_rot)
  co=Z_rot
  #compute variance explained
  pvar <- getpvar(ts, co, dat, npc) 
  res <- list("ts"=t(proc_rot), "co"=Z_rot, "eig"= NULL, 
              "pvar"=pvar, "dat"=dat, "npc"=npc, "mita"="DFA")
  return(res)
}


#ita.dpca ----------------------------------------
#' Run a Dynamic Principal Component Analysis
#' @param dat matrix with 
#' @param npc number of processes
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note lag set to 1, which correspond to AR1 processes
#' Function based on dpca from freqdom package
#' @references 
#' Ku W, Storer RH, Georgakis C (1995) Disturbance detection and isolation by dynamic principal component analysis. Chemom Intell Lab Syst 30:179–196
#' Ketelaere B De, Hubert M, Schmitt E (2015) Overview of PCA-based statistical process-monitoring methods for time-dependent, high-dimensional data. J Qual Technol 47:318–335
#' @export
#'
ita.dpca<-function(dat, npc){
  mvar <- freqdom::dpca(dat, q=1, Ndpc = npc)
  ts <- mvar$scores
  co <- apply(mvar$filters$operators, c(2,1), mean)
  #or lag0: mvar$filters$operators[,,2]
  eig <- mvar$var
  pvar <- getpvar(ts, co, dat, npc)
  res <- list("ts"=ts, "co"=co, "eig"= eig, "pvar"=pvar,
              "dat"=dat, "npc"=npc, "mita"="DPCA")
}

#ita.mafa ---------------------------------------
#' Run a Min/Max Autocorrelation Factors Analysis (MAFA)
#' @param dat matrix with 
#' @param npc number of processes
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use maf function created by M.Woillez and J. Rivoirard.
#' @references 
#' Solow AR (1994) Detecting Change in the Composition of a Multispecies Community. Biometrics 50:556
#' Woillez M, Rivoirard J, Petitgas P (2009) Using min/max autocorrelation factors of survey-based indicators to follow the evolution of fish stocks in time. Aquat Living Resour 22:193–200
#' @export
#' @seealso \code{\link{maf}}
#'
ita.mafa<-function(dat, npc){
  oldnames <- colnames(dat) 
  #replace column names with space
  colnames(dat)<- gsub(" ", "", oldnames)
  #create formula from names
  fo<-stats::as.formula(paste("~",sapply(list(colnames(dat)), paste, collapse="+"),sep=""))
  #compute MAF
  mvar <- maf(fo,as.data.frame(dat))
  ts <- as.matrix(mvar$x[,1:npc])
  co <- as.matrix(mvar$rotation[,1:npc])
  row.names(co) <- oldnames
  eig <- mvar$sdev**2
  pvar <- getpvar(ts, co, dat, npc)
  res <- list("ts"=ts, "co"=co, "eig"= eig, "pvar"=pvar, 
              "dat"=dat, "npc"=npc, "mita"="MAFA")
}

#ita.tsfa ---------------------------------------
#' Run a Time-series Factor Analysis
#' @param dat matrix with 
#' @param npc number of processes
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use function estTSF.ML from tsfa package.
#' @references 
#' Gilbert PD, & Meijer E (2005) Time series factor analysis with an application to measuring money. University of Groningen, Research School SOM Research Report 05F10.
#'
ita.tsfa<-function(dat, npc){
  if(ncol(dat)>nrow(dat)){
    warning("higher number of variables than number of years")
    return(list("ts"=NULL, "co"=NULL, "eig"= NULL, 
                "pvar"=NULL, "dat"=dat, "npc"=npc, "mita"="TSFA"))
  } else {
    mvar <- tsfa::estTSF.ML(dat, p = npc, diff. = F)
    pvar <- getpvar(mvar$f, mvar$loadings, dat, npc)
    res <- list("ts"=mvar$f, "co"=mvar$loadings, "eig"= NULL, 
                "pvar"=pvar, "dat"=dat, "npc"=npc, "mita"="TSFA")
  }
}

#ita.mds ----------------------------------------
#' Run Multidimensional scaling
#' @param dat matrix with 
#' @param npc number of processes
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} NULL - no scores of the variables on the new dimensions 
#' \item \code{eig} NULL 
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @export
#'
ita.mds <- function(dat, npc){
  mvar<- stats::cmdscale(stats::dist(dat), k=2, eig=TRUE)
  co <- crossprod(as.matrix(dat), as.matrix(mvar$points))
  #co <- wascores(mvar$points, dat) #from vegan package
  pvar <- getpvar(mvar$points, co, dat, npc) #compute variance explained
  res <- list("ts"=mvar$points, "co"=co, "eig"= mvar$eig, 
              "pvar"=pvar,"dat"=dat, "npc"=npc, "mita"="MDS")
}

#ita.lle ----------------------------------------
#' Run Local Linear Embedding
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} NULL - no scores of the variables on the new dimensions 
#' \item \code{eig} NULL 
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use function lle from the lle package.
#' @references 
#' Roweis ST (2000) Nonlinear Dimensionality Reduction by Locally Linear Embedding. Science (80- ) 290:2323–2326
#' @export
#'
ita.lle<-function(dat, npc){
  nk <- lle::calc_k(dat,npc, kmin=1, kmax=ncol(dat)-1, plotres = FALSE)
  mvar <- lle::lle(dat, m=npc, k=which.min(nk$rho), reg=2, ss=FALSE, id=TRUE, v=0.99)
  eig <- nk$rho # Not really eig !!
  ts <- as.matrix(mvar$Y)[,1:npc]
  #co <- NULL # LLE do not give variables scores
  #but one can compute the average scores on PC
  co <- as.matrix(crossprod(as.matrix(dat), as.matrix(ts)))
  #wascores(ts, dat) #from vegan package
  pvar <- getpvar(ts, co, dat, npc) #compute variance explained
  res <- list("ts"=ts, "co"=co, "eig"= eig, "pvar"=pvar, 
              "dat"=dat, "npc"=npc, "mita"="LLE")
}

#ita.ica ----------------------------------------
#' Run Independent Component Analysis
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} NULL - no scores of the variables on the new dimensions 
#' \item \code{eig} NULL 
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use function fastICA from the fastICA package.
#' @references 
#' Back AD, Weigend AS (1997) A First Application of Independent Component Analysis to Extracting Structure from Stock Returns. Int J Neural Syst 08:473–484
#' Blaschke T, Berkes P, Wiskott L (2006) What Is the Relation Between Slow Feature Analysis and Independent Component Analysis? Neural Comput 18:2495–2508
#' @export
#'
ita.ica<- function(dat, npc){
  mvar <- fastICA::fastICA(dat, n.comp=npc)
  tdat <- mvar$X
  co <- as.data.frame(t(mvar$A))
  #mvar$A change a lot depending on the initialization
  row.names(co) <- colnames(dat) 
  ts <- mvar$S
  pvar <- getpvar(ts, co, tdat, npc) #compute variance explained
  res <- list("ts"=ts, "co"=NULL, "eig"= NULL, "pvar"=pvar,
              "dat"=tdat, "npc"=npc, "mita"="ICA")
  return(res)
}

#ita.sfa ----------------------------------------
#' Run Slow Feature Analysis
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} NULL - no scores of the variables on the new dimensions 
#' \item \code{eig} NULL 
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use function from a package.
#' @references 
#' Wiskott L, Sejnowski TJ (2002) Slow Feature Analysis: Unsupervised Learning of Invariances. Neural Comput 14:715–770
#' Blaschke T, Berkes P, Wiskott L (2006) What Is the Relation Between Slow Feature Analysis and Independent Component Analysis? Neural Comput 18:2495–2508
#' @export
#'
ita.sfa<- function(dat, npc){
  # would be better to use rSFA package, sfa1 or sfa2 function
  # mvar <- sfa1(dat)
  mvar <- ForeCA::sfa(dat, n.comp = npc)#, spectrum.control = list(method = "wosa"))
  ts <- mvar$scores[,1:npc]
  co <- mvar$weightvectors[,1:npc] #or mvar$loadings or mvar$loadings.normalized
  pvar <- getpvar(ts, co, dat, npc)
  res <- list("ts"=ts, "co"=co, "eig"= mvar$sdev**2, 
              "pvar"=pvar,"dat"=dat, "npc"=npc, "mita"="SFA")
  return(res)
}

#ita.fca ----------------------------------------
#' Run Forecastable Component Analysis
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} NULL - no scores of the variables on the new dimensions 
#' \item \code{eig} NULL 
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note Use function foreca from the ForeCA package.
#' Limited to matrices with more time steps than variables
#' @references 
#' Goerg GM (2013) Forecastable Component Analysis. In: Proceedings of the 30th International Conference on Machine Learning.p 64–72
#' @export
#'
ita.fca<- function(dat, npc){
  if(ncol(dat)>nrow(dat)){
    warning("ForeCA can only be calculated with matrices with more time steps than variables.")
    return(list("ts"=NULL, "co"=NULL, "eig"= NULL, 
                "pvar"=NULL, "dat"=dat, "npc"=npc, "mita"="FCA"))
  } else {
    ret <- ts(dat) #no diff here, ok?
    mvar <- ForeCA::foreca(ret, n.comp = ncol(dat)-1)#, spectrum.control = list(method = "wosa"))
    ts <- as.matrix(mvar$scores[,1:npc])
    co <- as.matrix(mvar$weightvectors[,1:npc]) #or mvar$loadings or mvar$loadings.normalized
    pvar <- getpvar(ts, co, dat, npc)
    res <- list("ts"=ts, "co"=co, "eig"= mvar$sdev**2, 
                "pvar"=pvar,"dat"=dat, "npc"=npc, "mita"="FCA")
    return(res)
  }
}


#ita --------------------------------------------
#' Run Integrated Trend Analysis
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @param met selected ITA method
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the culumns are log-transformed (default=FALSE)
#' @details
#' Variable \code{met} are either string or number indicating: 
#' \itemize{
#' \item \code{ita.pca} scores of the time series on the new dimensions 
#' \item \code{ita.pcaV} scores of the variables on the new dimensions (if any)
#' \item \code{ita.} eigen values (for PCA associated methods)
#' \item \code{ita.} dataset used for ITA
#' \item \code{ita.} number of new dimensions
#' \item \code{ita.} method used for integrated trend analysis 
#' }
#' @return Output of ITA with 
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions 
#' \item \code{co} scores of the variables on the new dimensions (if any)
#' \item \code{eig} eigen values (for PCA associated methods)
#' \item \code{dat} dataset used for ITA
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @export
#' @examples
#' # Load data
#' data(baltic)
#' #  Run ita
#' mvar <- ita(baltic, npc=2)
#' # Plot eigen values
#' plot_eig(mvar)
#' # Plot scores of years
#' plot_timeseries(mvar)
#' # Plot scores of variables
#' plot_var(mvar)
#'
ita <- function(dat, npc, met=1, sca=TRUE, logt=FALSE){
  #Pre-process the dataset
  if (logt){
    dat <- log(dat+ifelse(min(dat)<=0, abs(min(dat))+1, 0))
  }
  if (sca){
    dat <- scale(dat) #rescale variables
    #remove scaled parameters
    attributes(dat)$`scaled:center` <- attributes(dat)$`scaled:scale` <- NULL
  }
  dat <- as.matrix(dat)
  
  #Multivariate analysis
  #Principal Component Analysis
  if(tolower(met)%in%c("1", "pca")){
    mvar <- ita.pca(dat, npc)
  }
  
  #Dynamic Factor Analysis
  if(tolower(met)%in%c("2", "dfa")){
    mvar <- ita.dfa(dat, npc)
  }
  
  #Dynamic PCA
  if(tolower(met)%in%c("3", "dpca")){
    mvar <- ita.dpca(dat, npc)
  }
  
  #Min/Max Autocorrelation Factors Analysis
  if(tolower(met)%in%c("4", "mafa", "maf")){
    mvar <- ita.mafa(dat, npc)
  }
  
  # Independent Component Analysis
  if(tolower(met)%in%c("5", "ica")){
    mvar <- ita.ica(dat, npc)
  }  
  
  # Slow Feature Analysis
  if(tolower(met)%in%c("6", "sfa")){
    mvar <- ita.sfa(dat, npc)
  }
  
  # Time-series Factor Analysis
  if(tolower(met)%in%c("7", "tsfa")){
    mvar <- ita.tsfa(dat, npc)
  }  
  
  #Forecastable Component Analysis
  if (tolower(met)%in%c("8", "fca", "foreca")){
    mvar <- ita.fca(dat, npc)
  }
  
  #Local Linear Embedding
  if (tolower(met)%in%c("9", "lle")){
    mvar <- ita.lle(dat, npc)
  }
  
  #PCA + Varimax
  if (tolower(met)%in%c("10", "pcav", "varimax", "pcavar")){
    mvar <- ita.pcaV(dat, npc)
  }
  
  #Multidimensional scaling 
  if (tolower(met)%in%c("11", "mds")){
    mvar <- ita.mds(dat, npc)
  }
  
  return(mvar)
}

