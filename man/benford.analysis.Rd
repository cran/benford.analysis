\docType{package}
\name{benford.analysis}
\alias{benford.analysis}
\alias{benford.analysis-package}
\title{Benford Analysis for data validation and forensic analytics}
\description{
  The Benford Analysis package provides tools that make it
  easier to validate data using Benford's Law. The main
  purpose of the package is to identify suspicious data
  that need further verification.
}
\details{
  More information can be found on its help documentation.

  The main function is \code{\link{benford}}. It generates
  a \code{Benford} S3 object.

  The package defines S3 methods for plotting and printing
  Benford type objects.

  After running \code{benford} you can easily get the
  "suspicious" data by using the functions:
  \code{\link{suspectsTable}}, \code{\link{getSuspects}},
  \code{\link{duplicatesTable}} and
  \code{\link{getDuplicates}}. See help documentation and
  examples for further details.

  The package also includes 6 real datasets for
  illustration purposes.
}
\examples{
data(corporate.payment) #gets data
cp <- benford(corporate.payment$Amount) #generates benford object
cp #prints
plot(cp) #plots
suspectsTable(cp) #prints the digits by decreasing order of discrepancies
suspects <- getSuspects(cp, corporate.payment) #gets suspicious observations
duplicatesTable(cp) #prints the duplicates by decreasing order
duplicates <- getDuplicates(cp, corporate.payment) #gets the duplicates
}
\references{
  Alexander, J. (2009). Remarks on the use of Benford's
  Law. Working Paper, Case Western Reserve University,
  Department of Mathematics and Cognitive Science. \cr\cr
  Berger, A. and Hill, T. (2011). A basic theory of
  Benford's Law. Probability Surveys, 8, 1-126. \cr\cr
  Hill, T. (1995). A statistical derivation of the
  significant-digit law. Statistical Science, 10(4),
  354-363. \cr\cr Nigrini, M. J. (2012). Benford's Law:
  Application for Forensic Accounting, Auditing and Fraud
  Detection. Wiley and Sons: New Jersey. \cr\cr Nigrini, M.
  J. (2011). Forensic Analyticis: Methods and Techniques
  for Forensic Accounting Investigations.Wiley and Sons:
  New Jersey. \cr\cr
}

