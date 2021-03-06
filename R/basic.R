#Toolbox of simple functions that makes life easier.

# ar1 -------------------------------------------
#' Compute first autocorrelation (AR1) of a time series
#'
#' @param ts time series
#' @return value of autocorrelation of lag 1
#' @examples 
#' data(baltic)
#' ar1(baltic$Cod)
#' @seealso \code{\link{infoTS}}
#' @export
#'
ar1 <- function(ts){
  return(stats::acf(ts, plot = FALSE)$acf[2])
}

# trend -------------------------------------------
#' Compute the linear trend of a time series
#'
#' @param ts time series
#' @return value of the linear trend
#' @examples
#' data(baltic)
#' trend(baltic$Cod)
#' @export
#'
trend <- function(ts){
  xseq <- 1:length(ts)
  lm1 <- lm(ts~xseq)
  return(lm1$coefficients[2])
}

# infoTS -------------------------------------------
#' Return a set of statistics on the original dataset
#'
#' @param dat matrix of time series (time in rows, variables in columns) 
#' @return characteristics of the dataset, including
#' \describe{
#'   \item{nbio}{number of biological time series}
#'   \item{ntime}{length of the time series}
#'   \item{meanAR}{mean AR1 of the time series}
#'   \item{maxAR}{maximum AR1 among the different time series}
#'   \item{meanCor}{mean pairwize correlation among time series}
#'   \item{maxCor}{maximum pairwize correlation among time series}
#' }
#' @examples
#' data(baltic)
#' infoTS(baltic)
#' @seealso \code{\link{ar1}}
#' @export
#'
infoTS <- function(dat){
  info <- c()
  #Characteritics of the matrix
  info$nbio <- ncol(dat) 
  info$ntime <- nrow(dat)
  #Autocorrelation 1
  autocorr <- apply(dat,2,ar1)
  info$meanAR <- mean(autocorr)
  info$maxAR <- max(autocorr)
  #Pairwize correlation
  paircor <- stats::cor(dat)
  #remove diagonal
  diag(paircor) <- NA
  info$meanCor <- mean(abs(paircor), na.rm=TRUE)
  info$maxCor <- max(abs(paircor), na.rm=TRUE)
  return(info)
}

# q95 -------------------------------------------
#' Compute the 0.95 quantile of a vector
#'
#' @param dat vecor
#' @return quantile(dat, prob=0.95)
#' @export
#'
q95 <- function(dat){
  return(quantile(dat, prob=0.95))
}

# short -------------------------------------------
#' Shorten strings taking the first character of words
#'
#' @param dat vector of characters
#' @param nchar number of characters to keep per word (default nchar=3)
#' @return shorten characters
#' @examples
#' data(baltic)
#' short(names(baltic))
#' @export
#'
short <- function(dat, nchar = 3){
  #find the separator
  sepL <- c(length(grep("\\.", dat)), length(grep("_", dat)), length(grep(" ", dat)))
  sep <- c("\\.", "_", " ")[which.max(sepL)]
  #split
  return(unlist(lapply(strsplit(dat, sep), shorten)))
}

#' @rdname short
shorten <- function(dat, nchar = 3){
  if (length(dat)>1){
    new <- paste(substr(dat[1], 1, nchar), substr(dat[2], 1, nchar), sep = "_")
  } else {
    new <- substr(dat[1], 1, nchar*2)
  }
  return(new)
}

# negPC -------------------------------------------
#' Homogenized the sign of PC trends to negative trends
#' (internal function)
#'
#' @param mvar results of integrative trend analysis (from ita function)
#' @return same object as input, but with scores multiplied by -1 if needed to get a negative trend 
#' @export
#'
negPC <- function(mvar){
  if (!is.null(mvar$ts)){
  # check if slope of regression >0
    for (i in 1:mvar$npc){
      if(trend(mvar$ts[,i])>0){
        mvar$ts[,i] <- -mvar$ts[,i]
        if (!is.null(mvar$co)){
          mvar$co[,i] <- -mvar$co[,i]
        }
      }
    }
  } else {
    warning("ts is null")
  }
  return(mvar)
}

# panel.cor.m -------------------------------------------
#' Function for plotting pairwize correlation
#'
#' @param x,y x and y
#' @param digits precision
#' @return Pearson correlation
#' @note Inspired from https://pbil.univ-lyon1.fr/
#' @export
#'
panel.cor.m <- function(x, y, digits=2)
{
  usr <- par("usr"); on.exit(par(usr))
  graphics::par(usr = c(0, 1, 0, 1))
  r <- stats::cor(x, y,method = "pearson", use = "pairwise.complete.obs")
  r2 <- abs(cor(x, y,method = "pearson", use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  cex <- 0.8/graphics::strwidth(txt)
  test <- stats::cor.test(x,y,method = "pearson",use = "pairwise.complete.obs")
  # borrowed from printCoefmat
  #Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
  #                 cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  #                 symbols = c("***", "**", "*", ".", " "))
  graphics::text(0.5, 0.5, txt, cex = cex * r2)
  #text(.8, .8, Signif, cex=cex, col=2)
}

# getpvar -------------------------------------------
#' Function to know the percentage of variance explained
#' (internal function)
#'
#' @param ts scores of the time series on the new dimensions 
#' @param co scores of the variables on the new dimensions 
#' @param dat matrix with 
#' @param npc number of new dimensions
#' @return the cumulated percentage of variance explained by the successive new dimensions
#' @note inspired from https://ro-che.info/articles/2017-12-11-pca-explained-variance
#' @export
#'
#
getpvar <- function(ts, co, dat, npc){
  pvar <- c()
  for (i in 1:npc){
    pred <- as.matrix(ts[,1:i]) %*% as.matrix(t(co[,1:i]))
    pvar <- c(pvar, stats::cor(as.vector(pred), as.vector(unlist(dat)))**2)
  }
  return(pvar)
}
