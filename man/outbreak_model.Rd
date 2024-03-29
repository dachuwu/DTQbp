% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outbreak_model.R
\name{outbreak_model}
\alias{outbreak_model}
\title{Run a single instance of the branching process model}
\usage{
outbreak_model(
  num.initial.cases = NULL,
  cap_gen = NULL,
  cap_max_days = NULL,
  cap_cases = NULL,
  r0isolated = NULL,
  r0community = NULL,
  disp.iso = NULL,
  disp.com = NULL,
  prop.asym = NULL,
  relR.asym = NULL,
  prop.ascertain = NULL,
  detect_sen = NULL,
  prop.presym = NULL,
  si_omega = NULL,
  incu_shape = NULL,
  incu_scale = NULL,
  delay_shape = NULL,
  delay_scale = NULL,
  delay_off = 0,
  quarant.days = NULL,
  quarant.retro.days = NULL,
  save.case_data = F
)
}
\arguments{
\item{num.initial.cases}{number of initial cases}

\item{cap_gen}{maximum number of generations to simulate}

\item{cap_max_days}{maximum number of days to simulate}

\item{cap_cases}{maximum number of cases to simulate}

\item{r0isolated}{reproduction number of isolated cases (must be >0)}

\item{r0community}{reproduction number of active cases in the community (must be >0)}

\item{disp.iso}{dispersion parameter for isolated cases (must be >0)}

\item{disp.com}{dispersion parameter for active cases in the community  (must be >0)}

\item{prop.asym}{proportion of asymptomatic cases}

\item{relR.asym}{relative transmissibility of asymptomatic cases}

\item{prop.ascertain}{ascertainment probability (probability of being contraced)}

\item{detect_sen}{detection probability}

\item{prop.presym}{proportion of pre-symptomatic transmission}

\item{si_omega}{standard deviation of the generation interval}

\item{incu_shape}{shape parameter of incubation period distribution}

\item{incu_scale}{scale parameter of incubation period distribution}

\item{delay_shape}{shape parameter of onset-to-test delay distribution}

\item{delay_scale}{scale parameter of onset-to-test delay distribution}

\item{delay_off}{offset of onset-to-test delay distribution}

\item{quarant.days}{quarantine period (days)}

\item{quarant.retro.days}{back-tracking period (days). When a contact is traced, s/he will be isolated if s/he had recall symptoms within the back-tracking period.}

\item{save.case_data}{True/False (default) whether to export case line list data}
}
\value{
data.table of cases by week, cumulative cases, and the effective reproduction number of the outreak
}
\description{
Run a single instance of the branching process model
}
\examples{

\dontrun{

}
}
