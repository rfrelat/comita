# Define the multiple methods to compute Integrative Trend Analysis
# ita.pca ---------------------------------------
#' Run a Principal Component analysis
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
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
ita.pca <- function(dat, npc) {
  #use ade4 package to compute pca #require(ade4)
  mvar <- ade4::dudi.pca(
    dat,
    nf = npc,
    scannf = FALSE,
    center = FALSE,
    scale = FALSE
  )
  pvar <- getpvar(mvar$li, mvar$co, dat, npc)
  res <- list(
    "ts" = mvar$li,
    "co" = mvar$co,
    "eig" = mvar$eig,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "PCA"
  )
  return(res)
}

# ita.pcaV ---------------------------------------
#' Run a Principal Component analysis followed by Varimax
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
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
ita.pcaV <- function(dat, npc) {
  mvar <- ade4::dudi.pca(
    dat,
    nf = npc,
    scannf = FALSE,
    center = FALSE,
    scale = FALSE
  )
  rawLoadings <- as.matrix(mvar$co[, 1:npc])
  vari <- stats::varimax(rawLoadings)
  ts <- scale(mvar$li[, 1:npc]) %*% vari$rotmat
  co <- unclass(vari$loadings)
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    "ts" = ts,
    "co" = co,
    "eig" = mvar$eig,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "PCA+Varimax"
  )
  return(res)
}

#ita.dpca ----------------------------------------
#' Run a Dynamic Principal Component Analysis
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of processes
#' @param q window size for the kernel estimator (default q=1)
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
#' Function based on dpca() from freqdom package
#' @references
#' Ku W, Storer RH, Georgakis C (1995) Disturbance detection and isolation by dynamic principal component analysis. Chemom Intell Lab Syst 30:179–196
#' Ketelaere B De, Hubert M, Schmitt E (2015) Overview of PCA-based statistical process-monitoring methods for time-dependent, high-dimensional data. J Qual Technol 47:318–335
#' @export
#'
ita.dpca <- function(dat, npc, q = 1) {
  mvar <- freqdom::dpca(dat, q = q, Ndpc = npc)
  ts <- mvar$scores
  co <- apply(mvar$filters$operators, c(2, 1), mean)
  #or lag0: mvar$filters$operators[,,2]
  eig <- mvar$var
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    "ts" = ts,
    "co" = co,
    "eig" = eig,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "DPCA"
  )
  return(res)
}

# ita.rpca ----------------------------------------
#' Run a Robust Principal Component Analysis (ROBPCA)
#'
#' Fits a robust principal component analysis with the Hubert et al.
#' (ROBPCA) algorithm, which down-weights outlying years so that a few
#' unusual observations do not dominate the components.
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of principal components to retain.
#'
#'
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the years on the robust components
#' \item \code{co} loadings of the variables on the robust components
#' \item \code{eig} robust eigenvalues
#' \item \code{pvar} cumulative proportion of variance explained
#' \item \code{dat} original dataset
#' \item \code{npc} number of retained components
#' \item \code{mita} method name, \code{"RPCA"}
#' }
#'
#' @note Requires the \pkg{rrcov} package. Accessor names
#'   (\code{getScores}, \code{getLoadings}, \code{getEigenvalues}) are taken
#'   from rrcov and may change between versions.
#'
#' @references
#' Hubert M, Rousseeuw PJ, Vanden Branden K (2005) ROBPCA: a new approach to
#' robust principal component analysis. Technometrics 47:64-79
#'
#' @seealso [ita.pca()], [ita()], [plot_eig()]
#' @export
ita.rpca <- function(dat, npc) {
  mvar <- rrcov::PcaHubert(as.matrix(dat), k = npc, scale = FALSE)
  ts <- rrcov::getScores(mvar)
  co <- rrcov::getLoadings(mvar)
  eig <- rrcov::getEigenvalues(mvar)
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = ts,
    co = co,
    eig = eig,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "RPCA"
  )
  return(res)
}

#ita.kpca ----------------------------------------
#' Run a Kernel Principal Component Analysis
#'
#' Non-linear PCA: data are mapped to a feature space with a kernel
#' function, and principal components are computed there.
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of components to retain.
#' @param kernel kernel name passed to [kernlab::kpca()] (default `"rbfdot"`).
#' @param sigma kernel width. If `NULL`, estimated with [kernlab::sigest()].
#'
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the years on the kernel components
#' \item \code{co} loadings of the variables, from \code{crossprod(dat, ts)}
#' \item \code{eig} eigenvalues in feature space
#' \item \code{pvar} cumulative proportion of variance explained
#' \item \code{dat} original dataset
#' \item \code{npc} number of retained components
#' \item \code{mita} method name, \code{"KPCA"}
#' }
#'
#' @details Eigenvalues are in feature space, so their percentages are not
#'   directly comparable with [ita.pca()]. Randomisation tests with
#'   [ita.rand()] should be interpreted with care.
#'
#' @note Requires the \pkg{kernlab} package and the functions
#'   \code{kpca}, \code{sigest}, \code{rotated}, \code{eig}).
#' @export
ita.kpca <- function(dat, npc, kernel = "rbfdot", sigma = NULL) {
  dat <- as.matrix(dat)
  if (is.null(sigma)) {
    sigma <- kernlab::sigest(dat, frac = 0.5)[2]
  } # median heuristic
  kp <- kernlab::kpca(
    dat,
    kernel = kernel,
    kpar = list(sigma = sigma),
    features = npc
  )
  ts <- kernlab::rotated(kp)
  eig <- kernlab::eig(kp)
  co <- as.matrix(crossprod(dat, ts))
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = ts,
    co = co,
    eig = eig,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "KPCA"
  )
  return(res)
}

#ita.pcoa ----------------------------------------
#' Run a Principal coordinates analysis (PCoA) with a Bray curtis distance
#' also called Multidimensional scaling (MDS)
#' @param dat matrix or data frame of time series (years in rows, variables
#'   in columns).
#' @param npc number of principal components to retain.
#' @param dist method to calculate distance (default = "bray")
#' @return Output of ITA with
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
ita.pcoa <- function(dat, npc, dist = "bray") {
  d <- vegan::vegdist(as.matrix(dat), method = dist)
  mvar <- stats::cmdscale(d, k = npc, eig = TRUE)
  ts <- mvar$points
  co <- as.matrix(crossprod(as.matrix(dat), ts))
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = ts,
    co = co,
    eig = mvar$eig[seq_len(npc)],
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "PCoA"
  )
  return(res)
}


#ita.nmds ----------------------------------------
#' Run Nonmetric Multidimensional scaling
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of processes
#' @param dist distance method (by default `bray`)
#' @param trymax maximum number of try
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} NULL - no scores of the variables on the new dimensions
#' \item \code{eig} NULL
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis
#' }
#' @note Use vegan::metaMDS and MASS::isoMDS
#' @references
#' Faith, D. P, Minchin, P. R. and Belbin, L. (1987). Compositional dissimilarity as a robust measure of ecological distance. Vegetatio 69, 57–68.
#' Kruskal, J.B. (1964). Nonmetric multidimensional scaling: a numerical method. Psychometrika 29, 115–129.
#' @export
#'
ita.nmds <- function(dat, npc, dist = "bray", trymax = 20) {
  ddat <- vegan::vegdist(dat, method = dist)
  # mvar <- MASS::isoMDS(ddat, k = npc, trace = FALSE)
  # co <- crossprod(as.matrix(dat), as.matrix(mvar$points))
  fit <- vegan::metaMDS(
    ddat,
    k = npc,
    trymax = trymax,
    trace = FALSE,
    autotransform = FALSE,
    engine = "isoMDS",
    wascores = FALSE
  )
  if (!isTRUE(fit$converged)) {
    warning("NMDS did not converge after ", trymax, " tries; check the stress")
  }

  ts <- as.matrix(fit$points)
  co <- as.matrix(crossprod(dat, ts))
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    "ts" = ts,
    "co" = co,
    "eig" = NULL,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "NMDS"
  )
  return(res)
}


#ita.dfa ----------------------------------------
#' Run a Dynamic Factor Analysis with diagonal and unequal matrix model
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
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
#' @note Use MARSS package, the R code is inspired from https://nwfsc-timeseries.github.io/atsa-labs/sec-dfa.html
#' Holmes, E. E., M. D. Scheuerell, and E. J. Ward. Applied time series analysis for fisheries and environmental data. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112.
#'
#' @references Zuur AF, Fryer RJ, Jolliffe IT, Dekker R, Beukema JJ (2003) Estimating common trends in multivariate time series using dynamic factor analysis. Environmetrics 14:665–685
#' @export
#'
ita.dfa <- function(dat, npc) {
  #parameters
  mm <- npc ## number of processes
  BB <- "identity" # diag(mm) ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  uu <- "zero" # matrix(0,mm,1)  ## 'uu' is a column vector of 0's
  ## 'CC' and 'cc' are for covariates
  CC <- "zero" # matrix(0,mm,1)
  cc <- "zero" # matrix(0,1,wk_last)
  QQ <- "identity" # diag(mm) ## 'QQ' is identity
  aa <- "zero" ## 'aa' is the offset/scaling
  ## 'DD' and 'd' are for covariates
  DD <- "zero" # matrix(0,mm,1)
  dd <- "zero" # matrix(0,1,wk_last)
  RR <- "diagonal and unequal" ## 'RR' is var-cov matrix for obs errors
  ## 'ZZ' is loadings matrix
  #tab <- tab[,6:10]
  Z_vals <- c()
  for (i in 1:ncol(dat)) {
    for (j in 1:npc) {
      Z_vals <- c(Z_vals, ifelse(j <= i, paste0("z", i, j), 0))
    }
  }
  ZZ <- matrix(Z_vals, nrow = ncol(dat), ncol = npc, byrow = TRUE)
  ## list with specifications for model vectors/matrices
  mod_list <- list(
    B = BB,
    U = uu,
    C = CC,
    c = cc,
    Q = QQ,
    Z = ZZ,
    A = aa,
    D = DD,
    d = dd,
    R = RR
  )
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 3000, allow.degen = TRUE)

  ## fit MARSS
  ret <- as.matrix(t(dat))
  mvar <- MARSS::MARSS(
    y = ret,
    model = mod_list,
    inits = init_list,
    control = con_list,
    silent = TRUE
  )

  ## get the estimated ZZ
  Z_est <- coef(mvar, type = "matrix")$Z
  ## get the inverse of the rotation matrix
  if (ncol(Z_est) > 1) {
    H_inv <- stats::varimax(Z_est)$rotmat
    ## rotate factor loadings
    Z_rot <- Z_est %*% H_inv
    ## rotate processes
    proc_rot <- solve(H_inv) %*% mvar$states
  } else {
    Z_rot <- Z_est
    proc_rot <- mvar$states
  }

  ts <- t(proc_rot)
  co <- Z_rot
  #compute variance explained
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    "ts" = t(proc_rot),
    "co" = Z_rot,
    "eig" = NULL,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "DFA"
  )
  return(res)
}


#ita.mafa ---------------------------------------
#' Run a Min/Max Autocorrelation Factors Analysis (MAFA)
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
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
ita.mafa <- function(dat, npc) {
  oldnames <- colnames(dat)
  #replace column names with space
  colnames(dat) <- gsub(" ", "", oldnames)
  #create formula from names
  fo <- stats::as.formula(paste(
    "~",
    sapply(list(colnames(dat)), paste, collapse = "+"),
    sep = ""
  ))
  #compute MAF
  mvar <- maf(fo, as.data.frame(dat))
  ts <- as.matrix(mvar$x[, 1:npc])
  co <- as.matrix(mvar$rotation[, 1:npc])
  row.names(co) <- oldnames
  eig <- mvar$sdev**2
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    "ts" = ts,
    "co" = co,
    "eig" = eig,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "MAFA"
  )
  return(res)
}


#ita.lle ----------------------------------------
#' Run Local Linear Embedding
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of new dimensions
#' @param k number of nearest neighbors
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} NULL - no scores of the variables on the new dimensions
#' \item \code{eig} NULL
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis
#' }
#' @note Use function embed from the dimRed package.
#' @references
#' Roweis ST (2000) Nonlinear Dimensionality Reduction by Locally Linear Embedding. Science (80- ) 290:2323–2326
#' @export
#'
ita.lle <- function(dat, npc, k = 10) {
  dat <- as.matrix(dat)
  k <- min(k, nrow(dat) - 1)
  emb <- dimRed::embed(dat, "HLLE", knn = k, ndim = npc)
  ts <- emb@data@data
  co <- as.matrix(crossprod(dat, ts))
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = as.matrix(ts),
    co = co,
    eig = NULL,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "LLE"
  )
  return(res)
}

#ita.ica ----------------------------------------
#' Run Independent Component Analysis
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
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
ita.ica <- function(dat, npc) {
  mvar <- fastICA::fastICA(dat, n.comp = npc)
  tdat <- mvar$X
  # mvar$A change a lot depending on the initialization
  # co <- as.data.frame(t(mvar$A))
  # row.names(co) <- colnames(dat)
  ts <- mvar$S
  # pvar <- getpvar(ts, co, tdat, npc) #compute variance explained
  res <- list(
    "ts" = ts,
    "co" = NULL,
    "eig" = NULL,
    "pvar" = NULL,
    "dat" = tdat,
    "npc" = npc,
    "mita" = "ICA"
  )
  return(res)
}

#ita.sfa ----------------------------------------
#' Run Slow Feature Analysis
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of new dimensions
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} scores of the variables on the new dimensions
#' \item \code{eig} NULL
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis
#' }
#' @note Use function `sfa1()` from package rSFA
#' @references
#' Wiskott L, Sejnowski TJ (2002) Slow Feature Analysis: Unsupervised Learning of Invariances. Neural Comput 14:715–770
#' Blaschke T, Berkes P, Wiskott L (2006) What Is the Relation Between Slow Feature Analysis and Independent Component Analysis? Neural Comput 18:2495–2508
#' @export
#'
ita.sfa <- function(dat, npc) {
  # would be better to use rSFA package, sfa1 or sfa2 function
  mvar <- rSFA::sfa1(dat, )
  ts <- as.matrix(mvar$y[, seq_len(npc), drop = FALSE]) # slow features (years x npc)
  co <- as.matrix(crossprod(dat, ts)) # same convention as ita.mds
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = ts,
    co = co,
    eig = NULL,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "SFA"
  )
  return(res)
}

#ita.fca ----------------------------------------
#' Run Forecastable Component Analysis
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of new dimensions
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} scores of the variables on the new dimensions
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
ita.fca <- function(dat, npc) {
  if (ncol(dat) > nrow(dat)) {
    warning(
      "ForeCA can only be calculated with matrices with more time steps than variables."
    )
    res <- list(
      "ts" = NULL,
      "co" = NULL,
      "eig" = NULL,
      "pvar" = NULL,
      "dat" = dat,
      "npc" = npc,
      "mita" = "FCA"
    )
  } else {
    ret <- ts(dat) #no diff here, ok?
    mvar <- ForeCA::foreca(ret, n.comp = ncol(dat) - 1) #, spectrum.control = list(method = "wosa"))
    ts <- as.matrix(mvar$scores[, 1:npc])
    co <- as.matrix(mvar$weightvectors[, 1:npc]) #or mvar$loadings or mvar$loadings.normalized
    pvar <- getpvar(ts, co, dat, npc)
    res <- list(
      "ts" = ts,
      "co" = co,
      "eig" = mvar$sdev**2,
      "pvar" = pvar,
      "dat" = dat,
      "npc" = npc,
      "mita" = "FCA"
    )
  }
  return(res)
}


#ita.tsne ----------------------------------------
#' Run t-distributed Stochastic Neighbor Embedding
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of new dimensions
#' @param p perplexity parameter (not bigger than 3 * perplexity < nrow(X) - 1)
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} NULL - no scores of the variables on the new dimensions
#' \item \code{eig} NULL
#' \item \code{dat} original dataset
#' \item \code{npc} number of new dimensions
#' \item \code{mita} method used for integrated trend analysis
#' }
#' @note Use function Rtsne from the Rtsne package.
#' @references
#' Maaten, L. Van Der, 2014. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research, 15, p.3221-3245.
#' van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
#' @export
#'
ita.tsne <- function(dat, npc, p = 2) {
  # perplexity must satisfy 3 * perplexity < nrow(dat) - 1
  maxp <- floor((nrow(dat) - 2) / 3)
  if (p > maxp) {
    warning(
      "perplexity ",
      p,
      " too large for ",
      nrow(dat),
      " rows; using ",
      maxp
    )
    p <- max(1, maxp)
  }
  mvar <- Rtsne::Rtsne(
    dat,
    dims = npc,
    perplexity = p,
    verbose = FALSE,
    check_duplicates = FALSE,
    max_iter = 500
  )
  ts <- mvar$Y
  co <- as.matrix(crossprod(as.matrix(dat), as.matrix(ts)))
  pvar <- getpvar(ts, co, dat, npc)
  res <- list(
    ts = ts,
    co = co,
    eig = NULL,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "TSNE"
  )
  return(res)
}

#ita.isomap ----------------------------------------
#' Run Isometric Feature Mapping Ordination
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of processes
#' @param k Number of shortest dissimilarities retained for a point (default k=3)
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} scores of the time series on the new dimensions
#' \item \code{co} scores of the variables on the new dimensions
#' \item \code{eig} NULL
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis
#' }
#' @export
#'
ita.isomap <- function(dat, npc, k = 3) {
  ddat <- vegan::vegdist(dat, method = "gower")
  mvar <- try(
    vegan::isomap(dist = ddat, ndim = npc, fragmentedOK = TRUE, k = k),
    silent = TRUE
  )
  if (inherits(mvar, "try-error")) {
    warning("isomap failed: ", conditionMessage(attr(mvar, "condition")))
    return(list(
      ts = NULL,
      co = NULL,
      eig = NULL,
      pvar = NULL,
      dat = dat,
      npc = npc,
      mita = "imap"
    ))
  }
  co <- crossprod(as.matrix(dat), as.matrix(mvar$points))
  pvar <- getpvar(mvar$points, co, dat, npc)
  res <- list(
    "ts" = mvar$points,
    "co" = co,
    "eig" = mvar$eig,
    "pvar" = pvar,
    "dat" = dat,
    "npc" = npc,
    "mita" = "imap"
  )

  return(res)
}

#ita.umap ----------------------------------------
#' Run Uniform Manifold Approximation and Projection (UMAP)
#'
#' Non-linear embedding of the years into a low-dimensional space, using
#' [uwot::umap()].
#'
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of dimensions to retain.
#' @param n_neighbors number of nearest neighbours. Keep small for short series
#'   (default 15, capped at `nrow(dat) - 1`).
#' @param min_dist minimum distance between embedded points (default 0.1).
#'
#' @return Output of ITA with
#' \itemize{
#' \item \code{ts} coordinates of the years in the embedding
#' \item \code{co} loadings of the variables, from \code{crossprod(dat, ts)}
#' \item \code{eig} NULL (no eigenvalues)
#' \item \code{pvar} cumulative proportion of variance explained
#' \item \code{dat} original dataset
#' \item \code{npc} number of dimensions
#' \item \code{mita} method name, \code{"UMAP"}
#' }
#'
#' @details Axes of a UMAP embedding have no variance ranking and are not
#'   defined up to a fixed rotation, so \code{co} is a projection rather than
#'   a property of the method.
#'
#' @note Requires the \pkg{uwot} package. With only a few dozen years, use
#'   small \code{n_neighbors} values; large values make the embedding
#'   uninformative.
#'
#' @references
#' McInnes L, Healy J, Melville J (2018) UMAP: Uniform Manifold Approximation
#' and Projection for Dimension Reduction. arXiv:1802.03426
#'
#' @export
ita.umap <- function(dat, npc, n_neighbors = 15, min_dist = 0.1) {
  ts <- uwot::umap(
    as.matrix(dat),
    n_components = npc,
    n_neighbors = min(n_neighbors, nrow(dat) - 1),
    min_dist = min_dist
  )
  co <- as.matrix(crossprod(as.matrix(dat), ts))
  pvar <- getpvar(ts, co, dat, npc)
  list(
    ts = ts,
    co = co,
    eig = NULL,
    pvar = pvar,
    dat = dat,
    npc = npc,
    mita = "UMAP"
  )
}


#ita --------------------------------------------
#' Run Integrated Trend Analysis
#' @param dat matrix or data frame of time series (years in rows, variables in columns).
#' @param npc number of new dimensions
#' @param met selected ITA method
#' @param sca whether columns of the original data are scaled (default =TRUE)
#' @param logt whether the columns are log-transformed (default=FALSE)
#' @param ... additional information
#' @details
#' Variable \code{met} is one of the method abbreviation among
#'   pca, dfa, dpca, mafa, ica, fca, varimax, mds, nmds, tsne, isomap
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
ita <- function(dat, npc, met = "pca", sca = TRUE, logt = FALSE, ...) {
  #Check method
  met <- tolower(met[1])

  if (!met %in% names(.itamet)) {
    stop(
      sprintf(
        "'%s' is not a recognised method. Use one of: %s",
        met,
        paste(unique(.itamet), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  umet <- as.character(.itamet[met])

  #Pre-process the dataset
  if (logt) {
    dat <- log(dat + ifelse(min(dat) <= 0, abs(min(dat)) + 1, 0))
  }
  if (sca & !umet %in% c("pcoa", "nmds")) {
    dat <- scale(dat) #rescale variables
    #remove scaled parameters
    attributes(dat)$`scaled:center` <- attributes(dat)$`scaled:scale` <- NULL
  }
  dat <- as.matrix(dat)

  #Multivariate analysis
  #Principal Component Analysis
  if (umet == "pca") {
    mvar <- ita.pca(dat, npc)
  }

  #PCA + Varimax
  if (umet == "pcaV") {
    mvar <- ita.pcaV(dat, npc)
  }

  #Dynamic PCA
  if (umet == "dpca") {
    mvar <- ita.dpca(dat, npc, ...)
  }

  #Robust PCA
  if (umet == "rpca") {
    mvar <- ita.rpca(dat, npc)
  }

  #Kernel PCA
  if (umet == "kpca") {
    mvar <- ita.kpca(dat, npc, ...)
  }

  #PCoA = MDS
  if (umet == "pcoa") {
    mvar <- ita.pcoa(dat, npc, ...)
  }

  #Nonmetric multidimensional scaling
  if (umet == "nmds") {
    mvar <- ita.nmds(dat, npc, ...)
  }

  #Dynamic Factor Analysis
  if (umet == "dfa") {
    mvar <- ita.dfa(dat, npc)
  }

  #Min/Max Autocorrelation Factors Analysis
  if (umet == "mafa") {
    mvar <- ita.mafa(dat, npc)
  }

  # Local Linear Embedding
  if (umet == "lle") {
    mvar <- ita.lle(dat, npc, ...)
  }

  # Independent Component Analysis
  if (umet == "ica") {
    mvar <- ita.ica(dat, npc)
  }

  # Slow Feature Analysis
  if (umet == "sfa") {
    mvar <- ita.sfa(dat, npc)
  }

  #Forecastable Component Analysis
  if (umet == "fca") {
    mvar <- ita.fca(dat, npc)
  }

  #t-distributed Stochastic Neighbor Embedding
  if (umet == "tsne") {
    mvar <- ita.tsne(dat, npc, ...)
  }

  #isomap
  if (umet == "isomap") {
    mvar <- ita.isomap(dat, npc, ...)
  }

  #umap
  if (umet == "umap") {
    mvar <- ita.umap(dat, npc, ...)
  }
  return(mvar)
}
