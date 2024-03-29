% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_transtree.R
\name{plot_transtree}
\alias{plot_transtree}
\title{Plot the transmission tree}
\usage{
plot_transtree(df, max.infector = 50, lyt = "kk")
}
\arguments{
\item{df}{case_data generated by outbreak_model}

\item{max.infector}{maximum number of infectors shown on the tree}

\item{lyt}{layout of the transmission tree. can be one of "kk"(defauk=lt),"tree", or ""}
}
\value{
a visNetwork object.
}
\description{
Plot the transmission tree
}
\examples{

\dontrun{
}
}
\author{
Ta-Chou Ng
}
