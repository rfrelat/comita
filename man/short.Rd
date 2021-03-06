% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic.R
\name{short}
\alias{short}
\alias{shorten}
\title{Shorten strings taking the first character of words}
\usage{
short(dat, nchar = 3)

shorten(dat, nchar = 3)
}
\arguments{
\item{dat}{vector of characters}

\item{nchar}{number of characters to keep per word (default nchar=3)}
}
\value{
shorten characters
}
\description{
Shorten strings taking the first character of words
}
\examples{
data(baltic)
short(names(baltic))
}
