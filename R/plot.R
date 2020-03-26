# Functions to create the different plots to compare methods for Integrative Trend Analysis

# plot_heatmap ---------------------------------------
#' Create a heatmap
#'
#' @param mvar results of multivariate analysis (from ita function)
#' @param col palette of colors (by default, traffic light colors)
#' @param n selected dimension to order the variables (default n=1) 
#' @param shortname short name of the variables, for plotting
#' @return A heatmap
#'
#' @author Romain Frelat, \email{rom.frelat@gmail.com}
#' @keywords heatmap
#' @examples
#' data(baltic)
#' # Run ita
#' mvar <- ita(baltic, npc=2)
#' # Create heatmap
#' plot_heatmap(mvar)
#' @export
#'
plot_heatmap <- function(mvar, col=c("springgreen4","chartreuse3","yellow","darkgoldenrod1", "red"), 
                         n=1, shortname=short(colnames(mvar$dat))){
  odPC <- order(mvar$co[,n])
  bk <- quantile(mvar$dat, probs = seq(0,1,length.out = length(col)+1))
  bk <- bk + c(-0.1, rep(0,length(col)-1), 0.1) #to be sure to include all values
  image(mvar$dat[,odPC], axes=FALSE, col = col, breaks=bk)
  box()
  axis(1, at = (seq(0,nrow(mvar$dat), length.out = nrow(mvar$dat)))/nrow(mvar$dat), 
       labels = row.names(mvar$dat))
  axis(2, at = (seq(0,ncol(mvar$dat), length.out = ncol(mvar$dat)))/ncol(mvar$dat), 
       labels = shortname[odPC], las=1, xpd=NA)
} 

# plot_timeseries ---------------------------------------
#' Plot time series of scores on the reduced dimensions
#'
#' @param mvar results of multivariate analysis (from ita function)
#' @param xax aestetic aspect, to plot xaxis on top of the graph 
#' @return A plot of time series
#'
#' @author Romain Frelat, \email{rom.frelat@gmail.com}
#' @keywords time series
#' @examples
#' # Load data
#' data(baltic)
#' #  Run ita
#' mvar <- ita(baltic, npc=2)
#' # Plot scores of years
#' plot_timeseries(mvar)
#' @export
#'
plot_timeseries <- function(mvar, xax=1){
  yr <- as.numeric(row.names(mvar$dat))
  limy <- range(mvar$ts, na.rm=TRUE)+c(-1,1)*(diff(range(mvar$ts, na.rm=TRUE))*0.05)
  limx <- range(yr)+c(-0.5, 0.5)
  plot(yr, mvar$ts[,1], type="l", ylab="PC score", 
       xaxt="n", ylim=limy, xlim=limx)
  if (xax==3){
    axis(3, at = yr, labels = rep("", length(yr)))
  } else {
    axis(1)
  }
  if(mvar$npc>1){
    for (i in 2:mvar$npc) lines(yr, mvar$ts[,i], lty=i)
  }
}

# plot_eig ---------------------------------------
#' Plot eigen values from multivariate analysis
#'
#' @param mvar results of multivariate analysis (from ita function)
#' @param leg whether legend is plotted (default = TRUE)
#' @param cum whether cumulative variance explained is plotted (default = FALSE) 
#' @param ... other argument for plot
#' @return A barplot of eigen vectors with percentage of variance explained in the successive components.
#' Also return an object with 
#' \itemize{
#' \item \code{barx} the x-positions of the bars in the ploting window
#' \item \code{perc} percentage of variance explained in the successive PC
#' }
#'
#' @keywords time series
#' @examples
#' # Load data
#' data(baltic)
#' # Run ita
#' mvar <- ita(baltic, npc=2)
#' # Plot scores of years
#' plot_eig(mvar)
#' @export
#' 
plot_eig <- function(mvar, leg=TRUE, cum=FALSE, ...){
  #check if eig is not null
  if (!is.null(mvar$eig)){
    coleig <- c(rep("black", mvar$npc), rep("grey", length(mvar$eig)-mvar$npc))
    perc <- mvar$eig/sum(mvar$eig)*100
    if (cum) {
      perc <- cumsum(perc)
    }
    barx <- barplot(perc, col=coleig, ylab="%", axisnames = FALSE, ...)
    if (leg){
      txtleg <- paste("PC", 1:(mvar$npc+1), ":", round(mvar$eig/sum(mvar$eig)*100)[1:(mvar$npc+1)], "%")
      legend("topright", txtleg, bty = "n", lty=c(1:mvar$npc,0))
    }
    return(data.frame("barx"=barx, "perc"=perc))
  } else {
    warnings("eig is null, no plot of eigen values can be created")
    plot.new()
  }
}

# plot_var --------------------------------------
#' Plot the scores of variables on reduced dimensions
#'
#' @param mvar results of multivariate analysis (from ita function)
#' @return A plot of scores of variables on on reduced dimensions
#'
#' @author Romain Frelat, \email{rom.frelat@gmail.com}
#' @keywords time series
#' @examples
#' # Load data
#' data(baltic)
#' #  Run ita
#' mvar <- ita(baltic, npc=2)
#' # Plot scores of variables
#' plot_var(mvar)
#' @note Use s.arrow from ade4 r-package if more than one axis, else a dotchart()
#' @export
#' 
plot_var <- function(mvar){
  lab <- colnames(mvar$dat)
  if (mvar$npc>1){
    limx <- range(mvar$co[,1])+ c(-1,1)*rep(diff(range(mvar$co[,1]))/7,2)
    ade4::s.arrow(mvar$co, label = lab, clabel = 0.8, xlim = limx)
  } else {
    odPC <- order(mvar$co[,1])
    dotchart(mvar$co[odPC], labels = lab)
  }
}

# plot_compts ---------------------------------------
#' Visualy compare the scores of the time series on reduced dimensions
#' from different methods of Integrated Trend Analysis
#'
#' @param dat original dataset
#' @param multimvar results of multiple ita (from multita function)
#' @param col pallette of colors for the different methods of Integrated Trend Analysis
#' @param pc which dimension to plot (default pc=1)
#' @param showleg whether to plot a legend (default showleg = TRUE)
#' @param mar4 additionnal margin on the right side of the plot
#' @param pccor compare best correlated PC, or fixed (default = TRUE)
#' @param ... additional arguments to plot()
#' @return The pairwize correlation of time series score compared with the first ITA of multimvar
#' @keywords time series
#' @examples
#' # Load data
#' data(baltic)
#' # Run ita
#' multi <- multita(baltic, npc=2, met=c("pca", "tsfa"))
#' # Compare scores of time series
#' plot_compts(baltic, multi)
#' @export
#'
plot_compts <- function(dat, multimvar, col=NULL, 
                        pc=1, showleg=TRUE, mar4=5,
                        pccor=TRUE,...){
  if(length(multimvar)==1){
    stop("use plot_timeseries() for single multivariate analysis")
  }
  yr <- as.numeric(row.names(dat))
  if(is.null(col)){
    if (length(multimvar)>8){
      getPalette<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
      col<-getPalette(length(multimvar))
    } else {
      col<-RColorBrewer::brewer.pal(length(multimvar), "Set2")
    }
  }
  if (showleg){
    op <- par(mar=par()$mar+c(0,0,0,mar4))
  }
  if (pccor){
    npc <- pc
    for (i in 2:length(multimvar)){
      corI <- cor(cbind(multimvar[[1]]$ts[,pc], multimvar[[i]]$ts))
      pcI <- which.max(abs(corI[-1,1]))
      npc <- c(npc, pcI)
      if(corI[-1,1][pcI]<0){
        multimvar[[i]]$ts[,pcI] <- -multimvar[[i]]$ts[,pcI]
      }
    }
  } else {
    npc <- rep(pc, length(multimvar))
  }
  
  #maxabs <- unlist(lapply(multimvar, function(x) max(abs(x$ts[,1]))))
  #limY <- c(-maxabs, maxabs)
  plot(yr, scale(multimvar[[1]]$ts[,npc[1]]), type="l", ylim=c(-2.5, 2.5), col=col[1],
       ylab="PC score", xlab="time", ...)
  if (length(multimvar)>1){
    for (i in 2:length(multimvar)){
      lines(yr, scale(multimvar[[i]]$ts[,npc[i]]), col=col[i])
    } 
  }
  
  #add a line y=0
  abline(h=0, col="grey70", lty=2)
  
  #legend
  leg <- c(paste0(names(multimvar)[1], npc[1]))
  coleg <-  col[1]
  alpha <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  stars <- c("***", "**", "*", ".", " ")
  corts <- c()
  if (length(multimvar)>1){
    for (i in 2:length(multimvar)){
      cor1i <- cor.test(multimvar[[1]]$ts[,npc[1]], multimvar[[i]]$ts[,npc[i]])
      leg <- c(leg, paste0(names(multimvar)[i], npc[i], ": ", round(cor1i$estimate,2), symnum(cor1i$p.value, alpha, stars)))
      corts <- c(corts, cor1i$estimate, cor1i$p.value)
      coleg <-  c(coleg, col[i])
    }
  }
  if (showleg){
    legend(par()$usr[2], sum(c(0.3,0.7)*par()$usr[3:4]), 
           lty = 1, col=col, cex=0.8, legend = leg, xpd=NA)
    names(corts) <- paste0(rep(c("cor", "pval"), length(multimvar)-1),
                           "_", rep(names(multimvar)[-1], each=2))
    par(op)
    return(corts)
  }
  else {
    return(list("leg"=leg, "col"=coleg))
  }
}

# plot_compvar ---------------------------------------
#' Visually compare the scores of the variables on reduced dimensions
#' from different methods of Integrated Trend Analysis
#'
#' @param dat original dataset
#' @param multimvar results of multiple ita (from multita function)
#' @param col pallette of colors for the different methods of Integrated Trend Analysis
#' @param pc which dimension to plot (default pc=1)
#' @param showleg whether to plot a legend (default showleg = TRUE)
#' @param shortname short name of the variables, use for plotting only
#'
#' @return The pairwize correlation of variable scores compared to the first ITA of \code{multimvar}.
#' @examples
#' #' # Load data
#' data(baltic)
#' # Run ita
#' multi <- multita(baltic, npc=2, met=c("pca", "tsfa"))
#' plot_compvar(baltic, multi)
#' @export
#'
plot_compvar <- function(dat, multimvar, col=NULL, 
                         pc=1, showleg=TRUE, shortname=short(colnames(dat))){
  lab <- shortname
  odPC <- order(multimvar[[1]]$co[,pc])
  
  if(is.null(col)){
    if (length(multimvar)>8){
      getPalette<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
      col<-getPalette(length(multimvar))
    } else {
      col<-RColorBrewer::brewer.pal(length(multimvar), "Set2")
    }
  }
  
  if (showleg){
    layout(matrix(1:2, ncol = 2), widths = c(3,1))
    par(mar=c(3,0,1,0))
  }
  
  dotchart(as.numeric(scale(multimvar[[1]]$co[odPC,pc])), 
           labels = lab[odPC], pch=16, xlim=c(-2.2, 2.2), 
           xlab=paste0("PC",pc))
  
  for (i in 1:length(multimvar)){
    if(!is.null(multimvar[[i]]$co)){
      points(scale(multimvar[[i]]$co[odPC,pc]), 1:nrow(multimvar[[i]]$co), pch=16, col=col[i])
    }
  }
  
  #Legend : pairwize correlation
  leg <- c(names(multimvar)[1])
  alpha <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  stars <- c("***", "**", "*", ".", " ")
  corbio <- c()
  coleg <-  col[1]
  if (length(multimvar)>1){
    for (i in 2:length(multimvar)){
      if(!is.null(multimvar[[i]]$co)){
        cor1i <- cor.test(multimvar[[1]]$co[,pc], multimvar[[i]]$co[,pc])
        leg <- c(leg, paste0(names(multimvar)[i], ": ", round(cor1i$estimate,2), symnum(cor1i$p.value, alpha, stars)))
        coleg <-  c(coleg, col[i])
        corbio <- c(corbio, cor1i$estimate, cor1i$p.value)
      }
    }
  }
  if (showleg){
    plot.new()
    legend("center", pch = 16, col=col, legend = leg, xpd=NA)
    names(corbio) <- paste0(rep(c("cor", "pval"), length(multimvar)-1),
                            "_", rep(names(multimvar)[-1], each=2))
    return(corbio)
  } else {
    return(list("leg"=leg, "col"=coleg))
  }
}

# plot_nrand ---------------------------------------
#' Compare the eigen values with randomly created matrices.
#'
#' @param mvar results of ita on original dataset (from ita function)
#' @param reig eig results from randomly created matrices (from ita.nrand function)
#' @param col color of the 95\% confidence interval of \code{reig}
#' @param cum whether cumulative variance explained is plotted (default = FALSE) 
#'
#' @return Plot of eigen values with 95\% confidence interval.
#' @examples
#' # Load data
#' data(baltic)
#' # Run ita
#' mvar <- ita(baltic, npc=2)
#' # Eigen values with multiple randomly generated data
#' mrand <- ita.nrand(baltic)
#' # Graphic comparing eigen values
#' plot_nrand(mvar, mrand)
#' 
#' # Or with multiple null model
#' randlist <- list("rand"=ita.nrand(baltic, metrand = "rand"),
#'                  "phase"=ita.nrand(baltic, metrand = "phase"))
#' plot_nrand(mvar, randlist, col=c("red", "blue"))
#' @export
#'
plot_nrand <- function(mvar, reig, col="blue", cum=FALSE){
  if(is.list(reig)){
    info <- plot_eig(mvar, leg=FALSE, cum=cum)
    if (length(reig)!= length(col)){
      warning("Number of colors doesn't match the number of null models")
      col <- rep(col, length.out=length(reig))
    }
    for (i in 1:length(reig)){
      eig95 <- apply(reig[[i]],2,q95)
      if (cum){
        cumeig <- t(apply(reig[[i]],1,cumsum))
        eig95 <- apply(cumeig,2,q95)
      }
      lines(info$barx, eig95, col=col[i], 
            type="b", pch=18, xpd=NA)
      signif <- info$perc>eig95
      star <- ifelse(signif, "*", "")
      text(x= info$barx, y= max(info$perc) + i*0.5, star, xpd=NA, col=col[i])
    }
    legend("topright", legend = names(reig), pch=18, col=col, bty="n")
  } else {
    #95% threshold
    eig95 <- apply(reig,2,q95)
    if (cum){
      cumeig <- t(apply(reig,1,cumsum))
      eig95 <- apply(cumeig,2,q95)
    }
    info <- plot_eig(mvar, leg=FALSE, cum=cum, 
                     ylim=c(0, max(c(mvar$eig/sum(mvar$eig)*100, eig95))))
    lines(info$barx, eig95, col=col, type="b", pch=18, xpd=NA)
    points(info$barx, eig95, col=col, pch=18, xpd=NA)
    signif <- info$perc>eig95
    star <- ifelse(signif, "*", "")
    text(x= info$barx, y= max(info$perc) + 1, star, xpd=NA, col=col)
  }
}

# plot_randtest ---------------------------------------
#' Compare the eigen values with randomly created matrices.
#'
#' @param obs eigen vectors from ita on original dataset (from ita function)
#' @param distri eig results from randomly created matrices (from ita.nrand function)
#' @param nclass number of bars in the histogram
#' @param coeff default =1
#' @param npc indicating which dimension to compare (default npc = 1)
#' @param ... other argument for plot
#'
#' @return Plot comparing the distribution of eigen values from random matrices, with the observed value.
#' @examples
#' # Load data
#' data(baltic)
#' # Run ita
#' mvar <- ita(baltic, npc=2)
#' #Eigen values with randomly generated data
#' mrand <- ita.nrand(baltic, nrep = 50)
#' # plot observed eigen value compared to null model 
#' plot_randtest(mvar, mrand)
#' @note inspired from ade4 r-package
#' @export
#'
plot_randtest <- function (obs, distri, nclass = 10, coeff = 1, npc=1, ...) 
{
  if ("eig"%in%names(obs)){
    obs <- obs$eig[npc]/sum(obs$eig)*100
  }
  if (ncol(distri)>1){
    distri <- distri[,npc]
  }
  r0 <- c(distri, obs)
  l0 <- max(distri) - min(distri)
  w0 <- l0/(log(length(distri), base = 2) + 1)
  w0 <- w0 * coeff
  xlim0 <- range(r0) + c(-w0, w0)
  h0 <- hist(distri, plot = FALSE, nclass = nclass)
  y0 <- max(h0$counts)
  plot(h0, xlim = xlim0, col = "grey80", ...)
  lines(c(obs, obs), c(y0/2, 0))
  points(obs, y0/2, pch = 18, cex = 2)
  invisible()
}

