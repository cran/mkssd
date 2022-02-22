\name{mkssd}
\alias{mkssd}
\title{Efficient multi-level k-circulant supersaturated designs}
\description{
		mkssd is a package that generates efficient balanced non-aliased multi-level k-circulant supersaturated designs by interchanging the elements of the generator vector. The package tries to generate a supersaturated design that has chisquare efficiency more than user specified efficiency level (mef). The package also displays the progress of generation of an efficient multi-level k-circulant design through a progress bar. The progress of 100 per cent means that one full round of interchange is completed. More than one full round (typically 4-5 rounds) of interchange may be required for larger designs.
	    }
\usage{mkssd(m,n,q,k,mef)}
\value{A list containing following items	
	\item{m}{number of factors}
	\item{n}{number of runs}
	\item{q}{number of levels}
	\item{k}{order of circulation}
	\item{generator.vector}{generator vector}
	\item{design}{design}
	\item{efficiency}{chi-square efficiency}
	\item{max.chisq}{maximum chi-square}
	\item{time.taken}{time taken to generate the design}
	\item{number.aliased.pairs}{number of aliased pairs}
	}
\arguments{
	\item{m}{number of factors}
	\item{n}{number of runs}
	\item{q}{number of levels}
	\item{k}{order of circulation}
	\item{mef}{minimum efficiency required, should be between 0 to 1}
	}
\references{B. N. Mandal, V. K. Gupta & Rajender Parsad (2014) Construction of Efficient Multi-Level k-Circulant Supersaturated Designs, Communications in Statistics - Theory and Methods, 43:3, 599-615}
\author{B N Mandal}
\keyword{supersaturated design}
\keyword{multi-level}
\keyword{k-circulant}
\keyword{efficiency}
\keyword{mkssd}
\examples{mkssd(10,6,3,2,1)}
