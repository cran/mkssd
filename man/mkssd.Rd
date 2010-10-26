\name{mkssd}
\alias{mkssd}
\title{ Efficient multilevel k-circulant supersaturated designs}
\description{
\code{mkssd} generates an efficient multilevel k-circulant supersaturated design for given 
number of factors (m), number of levels (q), number of runs (n) and order of circulant (k).
The algorithm is designed to generate a supersaturated design with no pair of aliased columns. 
However, if for a given parametric combination a supersaturated design without 
aliased pair of columns does not exist, then it may generate a supersaturated design
with aliased pair of columns.
}
\usage{
mkssd(m,n,q,k,mef)
}
\arguments{
\item{m}{number of factors}
\item{n}{number of runs}
\item{q}{number of levels}
\item{k}{order of circulant}
\item{mef}{minimum efficiency required, should be between 0 and 1}
}
\value{
a list containing:
\item{m}{number of factors}
\item{n}{number of runs}
\item{q}{number of levels}
\item{k}{order of circulant}
\item{genv}{generator vector of the supersaturated design}
\item{Deff}{chisquare efficiency of the design}
\item{max_chisq}{maximum chisquare value of the design. This is not printed if number of level (q) = 2}
\item{t_taken}{time taken to generate the design}
\item{aliased}{number of pairs of aliased columns in the design} 
}
\references{
B. N. Mandal, V.K. Gupta and Rajender Parsad. (2010). Construction of Efficient Multi-level k-circulant Supersaturated Designs
}
\author{
B N Mandal\email{  mandal.stat@gmail.com}
}
\examples{
##
mkssd(24,9,3,3,1)
}
\keyword{multilevel k-cirulant supersaturated design}