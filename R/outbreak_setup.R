#' Set up initial cases for branching process
#' @author Joel Hellewell
#'
#' @param num.initial.cases Integer number of initial cases
#' @param incfn function that samples from incubation period Weibull distribution; generated using dist_setup
#' @param delayfn function that samples from the onset-to-hospitalisation delay Weibull distribution; generated using dist_setup
#' @param k Numeric skew parameter for sampling the serial interval from the incubation period
#' @param prop.asym Numeric proportion of cases that are sublinical (between 0 and 1)
#'
#' @return data.table of cases in outbreak so far
#' @export
#'
#' @importFrom data.table data.table
#'
#' @examples
#'
#'\dontrun{
#' # incubation period sampling function
#' incfn <- dist_setup(dist_shape = 2.322737,dist_scale = 6.492272)
#' # delay distribution sampling function
#' delayfn <- dist_setup(delay_shape, delay_scale)
#' outbreak_setup(num.initial.cases = 5,incfn,delayfn,k=1.95,prop.asym=0)
#'}
outbreak_setup <- function(num.initial.cases, incfn, delayfn,  prop.asym) {
  # Set up table of initial cases
  inc_samples <- incfn(num.initial.cases)

  case_data <- data.table(exposure = rep(0, num.initial.cases), # Exposure time of 0 for all initial cases
                          asym = purrr::rbernoulli(num.initial.cases, prop.asym),
                          caseid = 1:(num.initial.cases), # set case id
                          infector = 0,
                          missed = TRUE,
                          onset = inc_samples,
                          generation = 0,
                          cluster = 1:(num.initial.cases),
                          R0 = NA,
                          Re = NA
                          )

  # set isolation/quarantine time for cluster to minimum time of onset of symptoms + draw from delay distribution
  case_data[,`:=`(
    isolated_time = onset + delayfn(.N),
    isolated = FALSE,
    quart_time = Inf,
    quart_end = Inf

  )]
  case_data$isolated_time[case_data$asym] <- Inf
  case_data[,`:=`(
     cluster_obs = data.table::fifelse(isolated_time==Inf, NA_integer_, cluster)
  )]

  # return
  return(case_data)
}

