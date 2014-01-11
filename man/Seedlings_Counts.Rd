\name{Seedlings_Counts}
\alias{Seedlings_Counts}
\docType{data}
\title{
Seedlings Incidence Counts Data
}
\description{
Seedlings incidence counts data is from Colwell and Coddington (1994).
}
\usage{data(Seedlings_Counts)}
\format{
Data must be arranged in the following order: (\emph{t}, \emph{r}, \emph{m}, 1, \emph{Q1}, 2, \emph{Q2}, ...,r, \emph{Qr}), where \emph{t} denotes the samples/quadrats, \emph{r} denotes the maximum frequency and \emph{m} denotes the number of values of \emph{fk}'s. \cr
The data entry for Seedlings data is \cr
(\emph{t}, \emph{r}, \emph{m}, 1, \emph{Q1}, 2, \emph{Q2}, ...,r, \emph{Qr}) = (121, 61, 18, 1, 3, 2, 2, ..., 52, 1, 61, 1); here, each number needs to be separated by at least one blank space. Here \emph{m} also equals the number of pairs that follow the number \emph{m} in the data sequence.
}
\details{
The type of Seedlings incidence counts data is called \strong{Incidence Counts Data for Multiple Samples/Quadrats}. \cr
The presemce/absence data for \emph{t} samples/quadrats are often summarized by incidence counts (\emph{Q1}, \emph{Q2},..., \emph{Qr}), where \emph{Qk} denotes the number of species that were detected in exactly emph{k} samples/quadrats in the data; and emph{r} denotes the number of samples/quadrats in which the most frequent species are found.
}
\source{
Holst, L. (1981). Some asymptotic results for incomplete multinomial or Poisson samples. \emph{Scandinavian Journal of Statistics} \strong(8), 243-246.
}
\examples{
data(Seedlings_Counts)
}
\keyword{datasets}
