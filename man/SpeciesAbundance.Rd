\name{SpeciesAbundance}
\alias{SpeciesAbundance}
\title{
Estimation of Species richness for individual-based (abundance) data
}
\description{
\kbd{SpeciesAbundance} is a function to provide a Species richness estimator proposed by Chao (2005) et al. \cr\cr
There are eight estimators which are Homogeneous Model(1992) estimators; Homogeneous MLE (1992) estimator; Chao1 (1984) estimator; Chao1-bc (2005) estimator; ACE (1992) estimator; Ace-1 (1992) estimator; 1st order jackknife (1978) estimators and the 2nd order jackknife (1978) estimators.}
\usage{
SpeciesAbundance(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), 
                 k = 10, conf = 0.95)
}
\arguments{
  \item{data}{
  a vector of species abundance frequency.}
  \item{method}{
  the \code{method} constructed to estimate entropy (see Details)
  }
  \item{k}{
  cut-off point. It is a value that separates frequency counts into abundant and rare groups
}
  \item{conf}{
  a positive number \eqn{\le} 1. "\code{conf}" specifies the confidence level for confidence   interval. The default is \code{0.95}.}
}
\details{
    \itemize{
      \item \code{method = "all"} : all estimators below involved. The default is "\code{all}"
      \item \code{method = "Homogeneous"} estimator, see Chao and Lee. (1992)
      \item \code{method = "Chao"} estimator, see Chao, Ann (1984 & 2005)
      \item \code{method = "CE"} estimator, see Chao and Lee. (1992)
      \item \code{method = "Jackknife"} estimator, see  Burnham, K. P. and Overton, W. S. (1978)
    }
  }
\value{
\kbd{SpeciesAbundance} returns a table of various species richness estimators, their standard error and 95 \% confidence interval which the \code{method} you choose.
}
\references{
Chao, A. and Lee, S.-M. (1992). Estimating the number of classes via sample converage. \emph{Journal of the American Statistical Association}, \bold{87}, 210-217. \cr\cr
Chao, A. (1984). Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics}, \bold{11}, 265-270. \cr\cr
Chao, A. (2005). Species estimation and applications. \emph{Encyclopedia of Statistical Sciences}, Second Edition, Vol. 12, 7907-7916 (N. Balakrishnan, C. B. Read and B. Vidakovic, Editors), Wiley, New York. \cr\cr
Burnham, K. P. and Overton, W. S. (1978). Estimaton of the size of a closed population when capture probabilities vary among naimals, \emph{Biometrika}, \bold{65}, 625-633. \cr\cr
}
\author{
K.S. Tseng ; Y.H. Lee : Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}}
\examples{
# load the individual-base (abundance) data
data(Spider_data)

# Estimation of species richness
SpeciesAbundance(Spider_data, method = "all", k=10, conf=0.95)
}