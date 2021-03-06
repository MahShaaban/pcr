% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_fun.R
\name{.pcr_plot_assess}
\alias{.pcr_plot_assess}
\title{Plot quality assessment graphs}
\usage{
.pcr_plot_assess(df, amount, reference_gene, method)
}
\arguments{
\item{df}{A data.frame}

\item{amount}{A numeric vector}

\item{reference_gene}{A character string}

\item{method}{A character string; 'efficiency' or 'standard_curve'}
}
\value{
A plot
}
\description{
Plot quality assessment graphs
}
\examples{
# locate and read file
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- read.csv(fl)

# make amount/dilution variable
amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

# plot the standard curves
pcr:::.pcr_plot_assess(ct3,
                amount = amount,
                reference_gene = 'GAPDH',
                method = 'standard_curve')

# plot amplification efficiency
pcr:::.pcr_plot_assess(ct3,
                       amount = amount,
                       reference_gene = 'GAPDH',
                       method = 'efficiency')

}
\keyword{internal}
