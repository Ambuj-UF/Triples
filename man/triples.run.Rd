\name{triples}
\alias{triples}
\title{SuperMatrix Rooted Triple (SMRT)}
\description{Calculates triples resolution from the input phylogeny tree file. }
\usage{triples(input = "", auth = "", outgroup = "", output)}
\arguments{
    \item{auth}{File that contains list of taxa to be included for 
    calculating triples resolution. If \code{auth = ""} (the default) 
    then the program terminates with message}
    \item{input}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{input = ""} (the default) then the
    program terminates with message}
    \item{outgroup}{Tree outgroup taxa name. If \code{outgroup = ""} (the default) 
    then the program terminates with message}
    \item{output}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{output = ""} (the default) programs generates
    output in Triples.txt file}
}

\references{
  DeGiorgio M, Degnan JH. Fast and Consistent Estimation of Species Trees Using Supermatrix Rooted Triples. Mol Biol Evol. 2010 Mar;27(3):552-69. 
  \url{http://mbe.oxfordjournals.org/content/27/3/552.long}

}

\author{Ambuj Kumar \email{ambuj@ufl.edu}}

\examples{triples(input = "", auth = "", outgroup = "", output)}

