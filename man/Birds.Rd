\name{Birds}
\alias{Birds}
\docType{data}
\title{
Birds data
}
\description{
Birds data  is from Magurran (1988, p,152). The data set contains 25 observed species.
}
\usage{data(Birds)}
\format{
  Their observed frequencies are repectively: \cr
 c(752, 276, 194, 126, 121, ..., 2, 2, 1, 1, 1) \cr\cr
 (that is, the first species is represented by 752 individuals, the second species by 276 individuals, ..., and the last three species are each represented by onlu one indivudual.)
}
\details{
The type of Birds data is called \strong{Species Frequency} or \strong{Abundance Data}.\cr
First, the frequency data must be  entered in a column.\cr
Second, if species names are recorded in your orginal data, they must be removed to conform to the data format shown below. \cr
Third, rows with \emph{frequency 0 for an unobserved species} may be included, but that species will not have any effect in the analysis.
}
\source{
Magurran, A. E. (1988). \emph{Ecological Diversity and Its Measurement.} Princeton, Princeton University Press, New Jersey.}
\examples{
data(Birds)
}
\keyword{datasets}
