\name{co2.new}
\alias{co2.new}
\docType{data}
\title{
Mauna Loa atmospheric CO2 concentration extended to 2010
}
\description{
Smoothed Mauna Loa annual CO2 concentrations in ppmv
}
\usage{data(co2)}
\format{
  The format is:
 num [1:191] 285 285 285  from 1820 to 2010: ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
NOAA; data initially collected by Charles Keeling
}
\references{
Keeling, C. D., et al. 1976, Atmospheric carbon dioxide variations at Mauna Loa Observatory, Hawaii, Tellus, 28, 538-551.

McCarroll D., Gagen M. H., Loader N. J., Robertson I., Anchukatis K. J., Los S., Young G. H. F., Jalkanen R., Kirchhefer A., and Waterhouse J. S., 2008, Correction of tree ring stable carbon isotope chronologies for changes in the carbon dioxide content of the atmosphere. GEOCHIMICA et COSMOCHIMICA ACTA, 73, 1539--1574, doi: 10.1016/j.gca.2008.11.041.

Thoning, K. W. et al., 1989, Atmospheric carbon dioxide at Mauna Loa Observatory 2. Analysis of the NOAA GMCC data, 1974-1985. Journal of Geophysical Researc, 94, 8549--8565.
}
\examples{
data(co2.new)
data(year.new)
plot(year.new, co2.new)
}
\keyword{datasets}
