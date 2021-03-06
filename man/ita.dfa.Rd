% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ita.R
\name{ita.dfa}
\alias{ita.dfa}
\title{Run a Dynamic Factor Analysis with diagonal and unequal matrix model}
\usage{
ita.dfa(dat, npc)
}
\arguments{
\item{dat}{matrix with}

\item{npc}{number of processes}
}
\value{
Result of multivariate analysis with 
\itemize{
\item \code{ts} scores of the time on the principal component
\item \code{co} scores of the variables on the principal component
\item \code{eig} eigen values
\item \code{dat} original dataset
\item \code{npc} number of selected principal component
\item \code{mita} method used for integrated trend analysis 
}
}
\description{
Run a Dynamic Factor Analysis with diagonal and unequal matrix model
}
\note{
Use MARSS pacakge, the R code is inspired from https://nwfsc-timeseries.github.io/atsa-labs/sec-dfa.html
Holmes, E. E., M. D. Scheuerell, and E. J. Ward. Applied time series analysis for fisheries and environmental data. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112.
}
\references{
Zuur AF, Fryer RJ, Jolliffe IT, Dekker R, Beukema JJ (2003) Estimating common trends in multivariate time series using dynamic factor analysis. Environmetrics 14:665–685
}
