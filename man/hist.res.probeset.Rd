\name{hist.res.probeset}
\alias{hist.res.probeset}
\title{
  Histogram of probe residuals in a dataset
}
\description{
  Generate a histogram of probe residuals, for visual inspection of empirical
  distribution.
}
\usage{
hist.res.probeset(batch, res = apply.res.probeset(batch), main = "",
                  xlab = "Probe residuals", ylab = "Density", ...)
}
\arguments{
  \item{batch}{
    an 'AffyBatch' object containing the probe intensities to analyze.
  }
  \item{which}{
    indices of arrays in the batch to analyze; one plot will be generated per array.
  }
  \item{res}{
    precomputed residual matrix; if NULL, residuals are computed as usual, based on the
    definition R = Y - Y_hat, where Y is the probe intensity and Y_hat is the mean
    intensity in the associated probe set.
  }
  \item{main}{
    main title for the plot.
  }
  \item{xlab}{
    x-axis label for the plot.
  }
  \item{ylab}{
    y-axis label for the plot.
  }
  \item{\dots}{
    further arguments passed on to \code{hist}.
  }
}
\value{
  object of class '"histogram"'; see \link{hist}.
}
\references{
  Serhal, P. and Lemieux, S. (2011). Elimination of spatial autocorrelation in microarray
  data. Bioinformatics, submitted.
}
\author{
  Philippe Serhal \email{philippe.serhal@umontreal.ca} \cr
  Sébastien Lemieux \email{s.lemieux@umontreal.ca}
}
\seealso{
  \link{normalize.mnf}, \link{image.mnf.psres}
}
\examples{
  # Load the sample dataset and generate histogram of probe residuals.
  data (gse2189, package = "mnf")
  hist.res.probeset (gse2189)
}
\keyword{ ~hplot }
\keyword{ ~mnf }
