\name{ChaoSpecies}
\alias{ChaoSpecies}
\title{
Estimation of Species richness
}
\description{
\kbd{ChaoSpecies} is a function to provide a Species richness estimator proposed by Chao (2005) et al. \cr\cr
For individual-based (abundance) data, there are eight estimators which are Homogeneous Model(1992) estimators; Homogeneous MLE (1992) estimator; Chao1 (1984) estimator; Chao1-bc (2005) estimator; ACE (1992) estimator; Ace-1 (1992) estimator; 1st order jackknife (1978) estimators and the 2nd order jackknife (1978) estimators. \cr\cr
For sample-based (incidence) data, there are seven other estimators which are Homogeneous Model(1994) estimators; Chao2 (1987) estimator; Chao2-bc (2005) estimator; ICE (1994) estimator; Ace-1 estimator; 1st order jackknife (1978) estimators and the 2nd order jackknife (1978) estimators.
}
\usage{
ChaoSpecies(data, datatype = c("abundance", "incidence"), k = 10, conf = 0.95)
}
\arguments{
  \item{data}{
   a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.}
  \item{datatype}{
  the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).}
   \item{method}{
  the \code{method} constructed to estimate entropy (see Details).}
  \item{k}{
  cut-off point. It is a value that separates frequency counts into abundant and rare groups
}
  \item{conf}{
  a positive number \eqn{\le} 1. "\code{conf}" specifies the confidence level for confidence   interval. The default is \code{0.95}.}
    \item{detail}{
  the basic data information of input data.}
}
\value{
\kbd{ChaoSpecies} returns a table of various species estimators, their standard error and 95 \% confidence interval which the \code{method} you choose.
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
\references{
Chao, A. and Lee, S.-M. (1992). Estimating the number of classes via sample converage. \emph{Journal of the American Statistical Association}, \bold{87}, 210-217. \cr\cr
Chao, A. (1984). Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics}, \bold{11}, 265-270. \cr\cr
Chao, A. (2005). Species estimation and applications. \emph{Encyclopedia of Statistical Sciences}, Second Edition, Vol. 12, 7907-7916 (N. Balakrishnan, C. B. Read and B. Vidakovic, Editors), Wiley, New York. \cr\cr
Burnham, K. P. and Overton, W. S. (1978). Estimaton of the size of a closed population when capture probabilities vary among naimals, \emph{Biometrika}, \bold{65}, 625-633. \cr\cr
}
\author{
K.S. Tseng ; Y.H. Lee : Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}
\seealso{
  \code{\link[ChaoSpecies]{SpeciesAbundance}}, 
  \code{\link[ChaoSpecies]{SpeciesIncidence}}
}
\examples{
##ex1
# load the individual-base (abundance) data
data(Birds)

# Estimation of species richness
ChaoSpecies(Birds, datatype="abundance", k=10, conf=0.95)


##ex2
# load the individual-base (abundance) data
data(Coin)

# Use the "DataTransform" function to transform the type of Coin in order to estimate the species richness.
Coin <- DataTransform(Coin, "FreqCount")

# Estimation of species richness
ChaoSpecies(Coin, datatype="abundance", k=10, conf=0.95)


##ex3
# load the sample-base (incidence) data
data(Seedlings_Pre)

# Use the "DataTransform" function to transform the type of Seedlings_Pre in order to estimate the species richness.
Seedlings_Pre <- DataTransform(Seedlings_Pre, "MatrixInci")

# Estimation of species richness
ChaoSpecies(Seedlings_Pre, datatype="incidence", k=10, conf=0.95)


##ex4
# load the individual-base (abundance) data
data(Seedlings_Freq)

# Use the "DataTransform" function to transform the type of Seedlings_Freq in order to estimate the species richness.
Seedlings_Freq <- DataTransform(Seedlings_Freq, "MatrixAbun")

# Estimation of species richness
ChaoSpecies(Seedlings_Freq, datatype="abundance", k=10, conf=0.95)


##ex5
# load the sample-base (incidence) data
data(Seedlings_Counts)

# Use the "DataTransform" function to transform the type of Seedlings_Counts in order to estimate the species richness.
Seedlings_Counts <- DataTransform(Seedlings_Counts, "InciCount")

# Estimation of species richness
ChaoSpecies(Seedlings_Counts, datatype="incidence", k=10, conf=0.95)


##ex6
# load the individual-base (abundance) data
data(Spider)

# Estimation of species richness
ChaoSpecies(Spider, datatype="abundance", k=10, conf=0.95)


##ex7
# load the sample-base (incidence) data
data(Ant)

# Estimation of species richness
ChaoSpecies(Ant, datatype="incidence", k=10, conf=0.95)
}
\keyword{SPECIES}
