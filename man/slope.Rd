% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope_and_angle.R
\name{slope}
\alias{slope}
\title{Computes slope}
\usage{
slope(time, volume, degree = TRUE)
}
\arguments{
\item{time}{A \code{vector} of time.}

\item{volume}{A \code{vector} of volume.}

\item{degree}{Default \code{TRUE} will give angle in degrees and \code{FALSE} will return in radians.}
}
\value{
Returns the slope and a \code{fit} object.
}
\description{
\code{slope} returns the slope for given time and volume data.
}
\examples{
time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
volume<- c(250.8, 320.4, 402.3, 382.6, 384, 445.9, 460.2, 546.8, 554.3, 617.9)
sl <- slope(time, volume)
par(pty="s")
xylimit <- range(c(time, volume))
plot(time, volume, type = "b", xlim = xylimit, ylim = xylimit)
abline(lm(volume~time))
}
