% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosquito.R
\name{HeFUN}
\alias{HeFUN}
\title{Calculate Heterozygosity}
\usage{
HeFUN(x, var = FALSE)
}
\arguments{
\item{x}{Microsatellite data.}

\item{var}{Should variance also be computed according to:
\deqn{2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2)/(n * (n - 1))}}
}
\value{
A list containing:
\item{H}{Point estimate of population heterozygosity}
\item{var.H}{Variance associated with the point estimate}
\item{n}{Number of samples used in the heterozygosity point estimate}
}
\description{
Calculate a heterozygosity estimate given a set of microsatellite data.
}
\details{
This function is based on the Heterozygosity estimate function within the
\pkg{pegas} package by Emmanuel Paradis. However Martin Donnelly believed
there was an error in the variance calculation so it was modified.
}
\examples{
Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02",
                "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07",
                "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")

MicroData <- formatMicrosat(HarrCraigData, Locusnames)

apply(MicroData[, 4:ncol(MicroData)], 2, HeFUN)
apply(MicroData[, 4:ncol(MicroData)], 2, HeFUN, var = TRUE)
}
\references{
Based on the Heterozygosity function in the \pkg{pegas}
package by Emmanuel Paradis.
}
\seealso{
[HeBootstrapFUN()]
}
\author{
Stu Field and Martin Donnelly
}
