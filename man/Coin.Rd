\name{Coin}
\alias{Coin}
\docType{data}
\title{
Coin data
}
\description{
Coin data is from Holst (1982) and Chao and Lee (1992): 204 coins were found in a hoard of ancient coins.
}
\usage{data(Coin)}
\format{
Data must be arranged in the following order: (\emph{r}, \emph{m}, 1, \emph{f1}, 2, \emph{f2}, ...,r, \emph{fr}), where \emph{r} denotes the maximum frequency and \emph{m} denotes the number of values of \emph{fk}'s. \cr
The data entry for Coin data is \cr
(\emph{r}, \emph{m}, 1, \emph{f1}, 2, \emph{f2}, ...,r, \emph{fr}) = (7, 7, 1, 102, 2, 26, ..., 6, 1, 7, 1); here, each number needs to be separated by at least one blank space. Here \emph{m} also equals the number of pairs that follow the number \emph{m} in the data sequence.
}
\details{
The type of Coin data is called \strong{Frequencies of Frequencies Data}. \cr
Species frequency data are often  classified by their frequencies into a simple form of (\emph{f1}, \emph{f2},..., \emph{fr}), where \emph{r} denotes the maximum frequency and \emph{fk} denotes the number of species represented by exactly \emph{k} individuals/times in the sample.
}
\source{
Holst, L. (1981). Some asymptotic results for incomplete multinomial or Poisson samples. \emph{Scandinavian Journal of Statistics} \strong(8), 243-246.
}
\examples{
data(Coin)
}
\keyword{datasets}
