\name{d13c}
\alias{d13c}
\docType{data}
\title{
delta carbon 13 concentrations in treerings
}
\description{
carbon 13 concentrations in treerings expressed as delta (promille) deviation from a measurement standard. Data are described in McCarroll et al (2008).
}
\usage{data(d13c)}
\format{
  The format is:
 num [1:183] -22.9 -23.7 -23.7 -23.4 -23.9  from 1820 to 2002: ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
McCarroll D., Gagen M. H., Loader N. J., Robertson I., Anchukatis K. J., Los S., Young G. H. F., Jalkanen R., Kirchhefer A., and Waterhouse J. S., 2008, Correction of tree ring stable carbon isotope chronologies for changes in the carbon dioxide content of the atmosphere. GEOCHIMICA et COSMOCHIMICA ACTA, 73, 1539--1574, doi: 10.1016/j.gca.2008.11.041.
}
\examples{
data(d13c)
data(year)
plot(year, d13c)
}
\keyword{datasets}
