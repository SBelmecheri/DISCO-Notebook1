pin <-
function (year, d13c, co2, alpha=0.9, firstyear=1820, lastyear=2002, plot.it=TRUE, verbose=TRUE, tol=0)
{
# pin - preindustrial correction for carbon isotope series for pCO2
#
# see McCarroll et al 2008 GEOCHIMICA et COSMOCHIMICA ACTA, 73, 1539--1574, doi: 10.1016/j.gca.2008.11.041 for details
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
# a fraction tol/100 is added or subtracted in case upper and lower bound need
# to be adjusted
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

