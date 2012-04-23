\name{image.pyn.psres}
\alias{image.pyn.psres}
\title{
  Spatial distribution of probe residuals
}
\description{
  Generate a heat map of probe residuals in a dataset, for visual inspection of empirical
  spatial distribution.
}
\usage{
image.pyn.psres(batch, which = 1:length(batch), transfo = log2, draw.legend = TRUE,
                shuffle = FALSE, ...)
}
\arguments{
  \item{batch}{
    an 'AffyBatch' object containing the probe intensities to analyze.
  }
  \item{which}{
    indices of arrays in the batch to analyze; one plot will be generated per array.
  }
  \item{transfo}{
    function to apply to each probe intensity before analysis.
  }
  \item{draw.legend}{
    if TRUE, draw a legend for the heat map, relating colours to values; if FALSE, don't.
  }
  \item{shuffle}{
    if TRUE, spatially shuffle the locations of the residuals before generating the
      image.
  }
  \item{\dots}{
    further arguments passed on to \code{image}.
  }
}
\value{
  NULL
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
  \link{normalize.pyn}, \link{hist.res.probeset}
}
\examples{
  # Load the sample dataset and generate heat map of probe residuals.
  data (gse2189, package = "pyn")
  image.pyn.psres (gse2189)
}
\keyword{ ~hplot }
\keyword{ ~pyn }
