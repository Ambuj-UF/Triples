\name{triples}
\alias{triples}
\title{SuperMatrix Rooted Triple (SMRT)}
\description{Calculates triples resolution from the input phylogeny tree file. }
\usage{triples(inpFile = "", authorityFile = "", outgroup = "", output)}
\arguments{
    \item{authorityFile}{File that contains list of taxa to be included for 
    calculating triples resolution. If \code{authorityFile = ""} (the default) 
    then the program terminates with message}
    \item{inpFile}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{inpFile = ""} (the default) then the
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

\examples{triples(inpFile = "", authorityFile = "", outgroup = "", output)}

