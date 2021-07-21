
#' Run a single instance of the branching process model
#' @author Joel Hellewell
#' @inheritParams outbreak_step
#' @param delay_shape numeric shape parameter of delay distribution
#' @param delay_scale numeric scale parameter of delay distribution
#'
#' @return data.table of cases by week, cumulative cases, and the effective reproduction number of the outreak
#' @export
#'
#' @importFrom data.table rbindlist
#'
#' @examples
#'
#'\dontrun{
#'}
pf_outbreak_model <- function(num.initial.cases = NULL,
                              cap_gen = NULL,

                              r0isolated = NULL, r0community = NULL,
                              disp.iso = NULL, disp.com = NULL,

                              prop.asym = NULL, relR.asym = NULL,

                              prop.ascertain = NULL, detect_sen=NULL,
                              prop.presym = NULL, si_omega = NULL,

                              incu_shape = NULL, incu_scale = NULL,
                              delay_shape = NULL, delay_scale = NULL, delay_off = 0,

                              quarant.days = NULL, quarant.retro.days = NULL) {

  # from weighted-average community R0 -> symptomatic R0
  r0symp <- r0community/(1 - prop.asym + prop.asym*relR.asym)

  # Set up functions to sample from distributions
  incfn <- dist_setup(dist_shape = incu_shape,
                      dist_scale = incu_scale)

  delayfn <- function(n){
    x <- dist_setup(delay_shape, delay_scale, delay_off)(n)
    x <- pmax(x, 0.001)
    x/purrr::rbernoulli(n, detect_sen)
  }

  # calculate skewness from prop. pre-symptomatic transmission
  k <- presym2k(prop.presym)
  inf_fn <- purrr::partial(inf_fn_pre, k = k, omga = si_omega)

  # Initial setup
  case_data <- outbreak_setup(num.initial.cases = num.initial.cases,
                                 incfn = incfn,
                                 prop.asym = prop.asym,
                                 delayfn = delayfn)
  case_data[,c("exposure","R0","Re"):=NULL]

  # Set initial values for loop indices
  max.gen <- 0
  extinct <- FALSE

  # Model loop
  while (max.gen < cap_gen & !extinct) {

    case_data <- pf_outbreak_step(case_data = case_data,
                         disp.iso = disp.iso,
                         disp.com = disp.com,
                         r0isolated = r0isolated,
                         r0symp = r0symp,
                         incfn = incfn,
                         delayfn = delayfn,
                         inf_fn = inf_fn,
                         prop.ascertain = prop.ascertain,
                         prop.asym = prop.asym,
                         relR.asym = relR.asym,
                         quarant.days = quarant.days,
                         quarant.retro.days = quarant.retro.days)

    max.gen <- max(case_data$generation)
    extinct <- all(case_data$isolated)
  }
  return(case_data)

}
