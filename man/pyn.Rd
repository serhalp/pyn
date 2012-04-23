\name{normalize.pyn}
\alias{normalize.pyn}
\alias{pyn}
\title{
  Normalize microarray data
}
\description{
  Main function of the pyn package, which returns a corrected ('normalized') copy of the
  given AffyBatch.
}
\usage{
normalize.pyn(batch, interest = "probeset", bias = "spatial",
              features.i = NULL, features.b = NULL, ki = 2, kb = 20,
              summary.stat.i = "mean", summary.stat.b = "mean", res.pre = NULL,
              do.log = TRUE, do.exp = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{batch}{
    an 'AffyBatch' object containing the data to normalize.
  }
  \item{interest}{
    name of predefined interest space to use to compute residuals; only current valid
    option is "probeset".  Can also be NULL (see 'features.i' below).
  }
  \item{bias}{
    name of predefined bias space to use to estimate correction factor from residuals;
    only current valid option is "spatial".  Can also be NULL (see 'features.b' below).
  }
  \item{features.i}{
    Custom interest space coordinate map, to be used to compute residuals: probe
    intensities will mapped onto this coordinate space, and the residual of probe x will
    be computed as \code{Y_x - summary.stat.i(kNN(x))}, where \code{Y_x} is the intensity
    of x and \code{kNN(x)} returns the intensities of x's 'ki' nearest neighbours in the
    interest space.  Must be either a vector, a 1-column matrix, or a 2-column matrix;
    with the same number of rows as \code{exprs(batch)}.  Can also be NULL, if 'interest'
    is not NULL.
  }
  \item{features.b}{
    Custom bias space coordinate map, to be used to compute residuals: probe residuals
    will mapped onto this coordinate space, and the correction factor of probe x will be
    computed as \code{summary.stat.b(kNN(x))}, where \code{kNN(xx)} returns the residuals
    of probe x's 'kb' nearest neighbours in the bias space.  Must be either a vector, a
    1-column matrix, or a 2-column matrix; with the same number of rows as
    \code{exprs(batch)}.  Can also be NULL, if 'bias' is not NULL.
  }
  \item{ki}{
    Number of neighbours to pool to compute each residual.  May be ignored by some
    interest space definitions (e.g. 'interest' = "probeset").  See 'features.i' above.
  }
  \item{kb}{
    Number of neighbours to pool to compute each correction factor.  May be ignored by
    some bias space definitions.  See 'features.b' above.
  }
  \item{summary.stat.i}{
    Summary statistic function to use to pool neighbouring intensities, to compute
    residuals; either 'mean' or 'median'.  See 'features.i' above.
  }
  \item{summary.stat.b}{
    Summary statistic function to use to pool neighbouring residuals, to compute
    correction facto; either 'mean' or 'median'..  See 'features.b' above.
  }
  \item{res.pre}{
    Precomputed matrix of residuals; overrides usual computation.  If NULL, computes
    residuals as usual.
  }
  \item{do.log}{
    If TRUE, log-transform (base 2) probe intensities before doing anything else;
    otherwise, leave them unchanged.
  }
  \item{do.exp}{
    If TRUE, exponentiate (base 2) probe intensities just before returning them;
    otherwise, leave them unchanged.
  }
  \item{verbose}{
    If TRUE, prints a few, sparse messages useful for tracking progress; otherwise, prints
    nothing at all.
  }
}
\details{
}
\value{
  An 'AffyBatch' object in which the probe intensities have been corrected (or
  'normalized').
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
  \link{hist.res.probeset}, \link{image.pyn.psres}
}
\examples{
  # Load the sample dataset and apply spatial correction, returning a new AffyBatch.
  data (gse2189, package = "pyn")
  gse2189.pyn <- normalize.pyn (gse2189)

  # Run some diagnostics: residual distribution and spatial distribution, before and
  # after spatial correction.
  par (mfrow = c (1, 2))
  hist.res.probeset (gse2189)
  hist.res.probeset (gse2189.pyn)
  image.pyn.psres (gse2189)
  image.pyn.psres (gse2189.pyn)
}
\keyword{ ~manip }
\keyword{ ~AffyBatch }
\keyword{ ~normalize }
\keyword{ ~pyn }
\keyword{ ~bgcorrect }
\keyword{ ~background }
