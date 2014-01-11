\name{DataTransform}
\alias{DataTransform}
\title{
Transform the data.
}
\description{
\kbd{DataTransform} is a function to transform the data so that it can be caluculated by \kbd{ChaoSpecies} function. Because \kbd{ChaoSpecies} function can only input certain form of data.
}
\usage{
DataTransform(data, type = c("FreqCount", "MatrixInci", "MatrixAbun", "InciCount"))
}
\arguments{
  \item{data}{
  a vector or a matrix which is not a vector of species abundance or incidence frequency. (\kbd{ChaoSpecies} 
  can only imput a vector of species abundance or incidence frequency.)
}
  \item{type}{
  the type of data. That is frequencies of frequencies (frequency counts) data, presence/absence data for 
  multiple samples/quadrats (multiple incidence data), frequency data for multiple samples/quadrats or 
  incidence counts data for multiple samples/quadrats.
}
}
\details{
  \kbd{DataTransform} can transform four types of data.
    \itemize{
      \item \code{type = "FreqCount"} : Frequencies of frequencies (frequency counts) data. \cr\cr
      Data must be arranged in the following order: (\emph{r}, \emph{m}, 1, \emph{f1}, 2, \emph{f2},       
      ...,r, \emph{fr}), where \emph{r} denotes the maximum frequency and \emph{m} denotes the number   
      of values of \emph{fk}'s. \cr
      
      \item \code{type = "MatrixInci"} : Presence\ Absence data for multiple samples/ quadrats. (Multiple 
      incidence data) \cr\cr
      The data is arranged in a matrix. The presence of any species in a   
      sample is denoted by 1 and its absence is denoted by 0. \cr
      
      \item \code{type = "MatrixAbun"} : Frequencies data for multiple samples/ quadrats.\cr\cr
      The data is arranged in a matrix. Each data entry in the data matrix 
      denotes the sample abundance or frequency. \cr
      
      \item \code{type = "InciCount"} : Incidence counts data for multiple samples/ quadrats.\cr\cr
      Data must be arranged in the following order: (\emph{t}, \emph{r}, \emph{m}, 1, \emph{Q1}, 2, 
      \emph{Q2}, ...,r, \emph{Qr}), where \emph{t} denotes the samples/quadrats, \emph{r} denotes the 
      maximum frequency and \emph{m} denotes the number of values of \emph{fk}'s. \cr
   }
}
\author{
K.S. Tseng ; Y.H. Lee : Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}}
}
\examples{
##ex1
# load the individual-base (abundance) data
data(Coin)

# Use the "DataTransform" function to transform the type of Coin in order to estimate the species richness.
Coin <- DataTransform(Coin, "FreqCount")


##ex2
# load the sample-base (incidence) data
data(Seedlings_Pre)

# Use the "DataTransform" function to transform the type of Seedlings_Pre.
Seedlings_Pre <- DataTransform(Seedlings_Pre, "MatrixInci")


##ex3
# load the individual-base (abundance) data
data(Seedlings_Freq)

# Use the "DataTransform" function to transform the type of Seedlings_Freq.
Seedlings_Freq <- DataTransform(Seedlings_Freq, "MatrixAbun")


##ex4
# load the sample-base (incidence) data
data(Seedlings_Counts)

# Use the "DataTransform" function to transform the type of Seedlings_Counts.
Seedlings_Counts <- DataTransform(Seedlings_Counts, "InciCount")
}
