\name{pin}
\alias{pin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pre Industrial (PIN) correction of d13C vegetation isotope series for atmospheric pCO2 changes
}
\description{
     PIN corrects the decline in tree ring stable isotope ratios (d13C)
     over the industrial period (> 1850) resulting from increased
     atmospheric CO2 levels using estimates of pre-industrial values
     for the internal (leaf) CO2 concentration (ci).

}
\usage{
pin(year, d13c, co2, alpha = 0.9, firstyear = 1820, lastyear = 2002, plot.it = TRUE, verbose = TRUE, tol = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{
 the years corresponding to the d13C (tree ring) measurements.
}
  \item{d13c}{
 the delta 13 C values

}
  \item{co2}{
the ambient annual atmospheric co2 levels, the Mauna Loa data from Keeling are included as a proxy

}
  \item{alpha}{
the span value used by the loess routine; alpha is adjusted for the length of the time series to allow for compatibility of the algorithm between time series of different length
}
  \item{firstyear}{
  firstyear / lastyear allows selecting a subset of the data
}
  \item{lastyear}{
  firstyear / lastyear included for compatibility with MatLab
          routine
}
  \item{plot.it}{
if TRUE a plot will be generated of corrected and uncorrected data, otherwise the correction will be applied without generating the plots.

}
  \item{verbose}{
if TRUE plots showing details of the correction will be generated; the user will need to hit the return key either in the plot window (Windows platforms) or at the command line prompt (UNIX and look alikes)
}
  \item{tol}{
  tolerance interval
}
}
\details{
    The relationship between ambient co2 (ca) and ci is estimated
     using the following steps: a. (STEP 2) the d13C observations are
     split in a low frequency and high frequency component (graph 1 if
     verbose=TRUE); b. (STEP 3) low frequency ci (internal - within
     leaf - CO2 concentration) is calculated from low frequency d13C;
     c. (STEP 4 + 5) two scenarios for ci are calculated; one assumes
     constant ci/ca ratio, the other a constant difference between
     ci-ca; d. (STEP 6) the d13C values for the two scenarios are
     calculated (graph2 if verbose=TRUE);  e. (STEP 8) the increments
     are calculated for the low frequency d13C and for the two
     (constant ratio and constant difference) scenarios; f. (STEP 9)
     final increments are calculated; if an observed (low frequency)
     increment is positive it is retained else if the d13C increment is
     more negative than the ci-ca increment (also negative) the ci-ca
     increment is subtracted (resulting in a smaller negative number),
     otherwise the (low  frequency) d13C increment does not change
     (graph 3 if verbose=TRUE); g. (STEP 10 + 11) the cumulative of the
     final increment is calculated (also graph 3) and is added to the
     initial low frequency d13C value (1850 value if present, otherwise
     the first year of the  time series), this is the corrected low
     frequency variation in d13C (graph 4 shows correction if
     verbose=TRUE)); h. (STEP 12) the high frequency d13C values are
     added to the corrected low frequency d13C values.

     Details: see references below.

}
\value{
an n by 2 matrix: the first column contains the years and the second
          column the corrected d13C values

 a sequence of 4 graphs:

 a graph: a graph showing the uncorrected and corrected values of d13C
          with their low frequency variation

}
\references{
McCarroll D., Gagen M. H., Loader N. J., Robertson I., Anchukatis K. J., Los S., Young G. H. F., Jalkanen R., Kirchhefer A., and Waterhouse J. S., 2008, Correction of tree ring stable carbon isotope chronologies for changes in the carbon dioxide content of the atmosphere. GEOCHIMICA et COSMOCHIMICA ACTA, 73, 1539--1574, doi: 10.1016/j.gca.2008.11.041.
}
\author{
Sietse O. Los
}
\note{
 the algorithm uses the exact form of loess, control=loess.control(surface="direct", statistics="exact"), on some platforms the approximate version of loess leads to slightly different answers

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
'loess' and 'loess.control'
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (year, d13c, co2, alpha=0.9, firstyear=1820, lastyear=2002, plot.it=TRUE, verbose=TRUE, tol=0)
{
# pin - preindustrial correction for carbon isotope series for pCO2
#
# see McCarroll et al for details
#
# inputs year = 1D array of years
#        d13c = vector of atmospherically corrected d13 C values
#        co2  = smoothed atmospheric co2 record (e.g. from mauna loa)

#
# check for missing values and unequal length of input array
#
  if (any(is.na(c(year, d13c, co2))))
  {
    stop("Error: missing values in input data")
  }
  if ( ( max(d13c) > -15 ) | (min(d13c) < -30) )
  {
    print( paste("WARNING: Unusual values encountered in d13c data range is between", as.character(min(d13c)), "-", as.character(max(d13c)) ))
  }
  if ( !( (length(year) == length(d13c)) & (length(year)== length(co2)) ) )
  {
    stop("Error: d13c, co2 or year are of unequal length")
  }

#
# STEP 1.a adjust length arrays to common period
# provide a warning if length is adjusted
#
  if (firstyear > lastyear)
  {
    xxx <- lastyear
    lastyear <- firstyear
    firstyear <- xxx
    print("WARNING firstyear > lastyear; code reversed these")
    print("WARNING check input files and rerun if necessary")
  }
  firstyear.data <- max( min(year), firstyear)
  lastyear.data <- min( max(year), lastyear)
  if ( (firstyear.data != firstyear) | (lastyear.data != lastyear) )
  {
    print(paste("Warning: time series cover only", firstyear.data, "-", lastyear.data) )
  }
  if (year[1] > year[length(year)] )
  {
    print("Warning: year is in reversed order; all data are reordered accordingly")
    co2 <- co2[length(co2):1]
    d13c <- d13c[length(d13c):1]
    year <- year[length(year):1]
  }
  ca <- co2 <- co2[ (year >= firstyear) & (year <= lastyear) ]
  d13c <- d13c[ (year >= firstyear) & (year <= lastyear) ]
  year <- year[(year >= firstyear) & (year <= lastyear) ]
  n.yr <- length(year)

#
# STEP 1.a split d13c into a low and high frequency component, only low frequency component is corrected
#
  loess.d13c <- loess(d13c ~ year,  span=alpha, control=loess.control(statistics="exact", surface="direct") )
  d13c.low   <- predict(loess.d13c, year)
  d13c.high  <- d13c - d13c.low
  if (verbose)
  {
    old.ask <- par("ask")
    old.mfrow <- par("mfrow")
    par(ask=TRUE, mfrow=c(1,1))
    plot(year, d13c, type="l")
    lines(year, d13c.low, lwd=3, col="red")
    legend(year[1], max(d13c, na.rm=T), lwd=c(1,3), col=c("black", "red"), legend=c("d13C data", "d13C low frequency") )
  }

#
# STEP 2 convert low frequency component of d13C to ci
#
  a <- -4.4
  b <- -27
  deltac <- 6.4 + d13c.low
  cica <- (deltac - a) / (b-a)
  ci.actual <- ca * cica

#
# STEP 3 Define low frequency trends in ci as a function of ca
# choose exact form of loess because the approximate form deviates from the exact
#
  loess.ci <- loess( ci.actual ~ ca, span=alpha, control=loess.control(statistics="exact", surface="direct"))
  ci.predicted <- predict(loess.ci, ca)
  if (verbose)
  {
    plot( ca, ci.actual, type="l" )
    lines( ca, ci.predicted, col="red", lwd=3)
    legend( min(ca), max(ci.actual),  lwd=c(1,3), col=c("black", "red"), legend=c("ci actual", "ci predicted") )
  }

#
# STEP 4 Calculate annual increments in ca (1st logical constraint)
#
  max.increment <- diff(ca)
  loess.increment <- diff(ci.predicted)
#
# STEP 5 Define low-frequency trends in ci that would have occurred had ci/ca remained constant (2nd logical constraint)
#
  if (firstyear < 1850)
  {
    initial.ca <- 285
    initial.cica <- mean( cica[ year < 1850] )
  }
  else
  {
    initial.ca <- ca[1]
    initial.cica <- cica[1]
  }
  ci.constant <- initial.cica*ca
  min.increment <- diff(ci.constant)
#
# STEP 6 Incrementally correct the long-term trend in ci as a function of ca
#
  max.increment <- max.increment*(1+tol/100)
  min.increment <- min.increment*(1-tol/100)
  final.increment <- apply( cbind(loess.increment, max.increment), 1, min, na.rm=T)
  final.increment <- apply( cbind(final.increment, min.increment), 1, max, na.rm=T)
  if (verbose)
  {
    par(mfrow=c(2,1))
    plot(year, c(NA, loess.increment), ylim=range(c(min.increment, max.increment)), ylab="ci increment")
    lines(year, c(NA,  max.increment), lwd=2, col="blue", type="l")
    lines(year, c(NA, min.increment), lwd=2, col="red")
    lines(year, c(NA, final.increment), lwd=3, col="grey")
    legend(year[1], max(max.increment), lty=c(0,1,1,1), pch=c(1,0,0,0), col=c("black", "blue", "red", "grey"), lwd=c(1,2,2,3), legend=c("loess", "max increment", "min increment", "final increment") )
    plot(year, cumsum(c(ci.constant[1], loess.increment)), ylim=range(c(cumsum(max.increment)+ci.constant[1],cumsum( min.increment)+ci.constant[1])), ylab="cumulative increment")
    lines(year, cumsum(c(ci.constant[1], max.increment) ), lwd=2, col="blue")
    lines(year, cumsum(c(ci.constant[1], min.increment) ), lwd=2, col="red")
    lines(year, cumsum(c(ci.constant[1], final.increment)), lwd=3, col="grey")
  }
  cumulative.minimum <- cumsum(final.increment)
  corrected.pin <- -10.8 + (-22.6 * ((ci.actual[2:n.yr] - cumulative.minimum) / initial.ca))
  corrected.pinyear <- year[2:n.yr]
  splice <- c(1:n.yr)[year==min(corrected.pinyear)]
  corrected.d13c.low <- c( d13c.low[1:(splice-1)], corrected.pin)
  corrected.d13c     <- corrected.d13c.low+d13c.high
  corrected.year <- c( year[1:(splice-1)], corrected.pinyear)
  if (verbose)
  {
    plot(year, corrected.d13c.low, type="b", ylab="low frequency d13c", pch="*", ylim=range(c(d13c.low, corrected.d13c.low)) )
    lines(year, d13c.low, type="l")
    legend(year[1], min(corrected.d13c.low), lwd=c(1,1), lty=c(0,1), pch=c("*"," "), legend=c("after correction", "before correction") )
    title("Before and after correction")
    plot(year, d13c.high, type="l", ylab="high frequency d13c")
    title("d13c high frequency; Average d13c removed")

  }
  if (plot.it)
  {
    par(mfrow=c(1,1))
    plot(year, d13c, type="l", ylim=range(c(d13c, corrected.d13c), na.rm=T), lwd=3)
    lines(corrected.year, corrected.d13c, lwd=1, col="red")
    legend(year[1], max( c(d13c, corrected.d13c)), lwd=c(3,1), col=c("black", "red"), legend=c("d13c measured", "d13c corrected") )
    par(ask=old.ask, mfrow=c(old.mfrow))
  }
  cbind( corrected.year, corrected.d13c)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
