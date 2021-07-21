#' Run a specified number of simulations with identical parameters
#' @author Joel Hellewell
#' @param n.sim number of simulations to run
#' @param num.initial.cases Initial number of cases in each initial cluster
#' @param num.initial.clusters Number of initial clusters
#' @param prop.ascertain Probability that cases are ascertained by contact tracing
#' @param cap_max_days Maximum number of days to run process for
#' @param cap_cases Maximum number of cases to run process for
#' @param r0isolated basic reproduction number for isolated cases
#' @param r0community basic reproduction number for non-isolated cases
#' @param disp.iso dispersion parameter for negative binomial distribution for isolated cases
#' @param disp.com dispersion parameter for negative binomial distribution for non-isolated cases
#' @param delay_shape shape of distribution for delay between symptom onset and isolation
#' @param delay_scale scale of distribution for delay between symptom onset and isolation
#'
#' @importFrom purrr safely
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' res <- scenario_sim(n.sim = 5,
#' num.initial.cases = 5,
#' cap_max_days = 365,
#' cap_cases = 2000,
#' r0isolated = 0,
#' r0community = 2.5,
#' disp.iso = 1,
#' disp.com = 0.16,
#' k = 0.7,
#' delay_shape = 2.5,
#' delay_scale = 5,
#' prop.asym = 0,
#' prop.ascertain = 0)
#' #' }
#'
scenario_sim <- function(n.sim = NULL,
                         save.case_data = F,
                         num.initial.cases = NULL,
                         cap_gen = NULL, cap_max_days = NULL, cap_cases = NULL,
                         r0isolated = NULL, r0community = NULL,
                         disp.iso = NULL, disp.com = NULL,
                         prop.asym = NULL,relR.asym = NULL,

                         prop.ascertain = NULL,
                         prop.presym = NULL, si_omega = NULL,
                         incu_shape = NULL, incu_scale = NULL,
                         delay_shape = NULL, delay_scale = NULL,  delay_off = 0,
                         detect_sen=NULL,
                         quarant.days = NULL, quarant.retro.days = NULL) {

  # Run n.sim number of model runs and put them all together in a big data.frame
  res <- purrr::map(.x = 1:n.sim, ~ outbreak_model(num.initial.cases = num.initial.cases,
                                                   cap_gen = cap_gen,
                                                   cap_max_days = cap_max_days,
                                                   cap_cases = cap_cases,
                                                   r0isolated = r0isolated,
                                                   r0community = r0community,
                                                   disp.iso = disp.iso,
                                                   disp.com = disp.com,
                                                   prop.asym = prop.asym,
                                                   relR.asym = relR.asym,
                                                   prop.ascertain = prop.ascertain,
                                                   prop.presym = prop.presym,
                                                   si_omega = si_omega,
                                                   incu_shape = incu_shape,
                                                   incu_scale = incu_scale,
                                                   delay_shape = delay_shape,
                                                   delay_scale = delay_scale,
                                                   delay_off = delay_off,
                                                   detect_sen=detect_sen,
                                                   quarant.days = quarant.days,
                                                   quarant.retro.days = quarant.retro.days,
                                                   save.case_data = save.case_data
                                                   ))

  if(save.case_data){
    case_data <- purrr::map(res,1)
  } else {
    case_data <- list()
  }

  sum_data <- purrr::map(res,2) %>% data.table::rbindlist()

  return(list(case_data = case_data, sum_data = sum_data))
}
