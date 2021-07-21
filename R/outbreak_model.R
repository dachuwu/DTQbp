
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
#' incfn <- dist_setup(dist_shape = 2.322737,dist_scale = 6.492272)
#' # delay distribution sampling function
#' delayfn <- dist_setup(2, 4)
#' # generate initial cases
#' case_data <- outbreak_setup(num.initial.cases = 5,
#'                             incfn=incfn,
#'                             delayfn = delayfn,
#'                             k=1.95,
#'                             prop.asym=0)
#' # generate next generation of cases
#' case_data <- outbreak_step(case_data = case_data,
#'                            disp.iso = 1,
#'                            disp.com = 0.16,
#'                            r0isolated = 0,
#'                            r0community = 2.5,
#'                            prop.asym = 0,
#'                            incfn = incfn,
#'                            delayfn = delayfn,
#'                            prop.ascertain = 0,
#'                            k = 1.95,
#'                            quarantine = FALSE)
#'}
outbreak_model <- function(num.initial.cases = NULL,
                           cap_gen = NULL, cap_max_days = NULL, cap_cases = NULL,

                           r0isolated = NULL, r0community = NULL,
                           disp.iso = NULL, disp.com = NULL,

                           prop.asym = NULL, relR.asym = NULL,

                           prop.ascertain = NULL, detect_sen=NULL,
                           prop.presym = NULL, si_omega = NULL,

                           incu_shape = NULL, incu_scale = NULL,
                           delay_shape = NULL, delay_scale = NULL,  delay_off = 0,

                           quarant.days = NULL, quarant.retro.days = NULL,
                           save.case_data = F) {
  # from weighted-average community R0 -> symptomatic R0
  r0symp <- r0community/(1 - prop.asym + prop.asym*relR.asym)

  # Set up functions to sample from distributions
  # incubation period sampling function
  incfn <- dist_setup(dist_shape = incu_shape,
                      dist_scale = incu_scale)

  # incfn <- dist_setup(dist_shape = 3.303525,dist_scale = 6.68849) # incubation function for ECDC run
  # onset to isolation delay sampling function
  # delayfn <- dist_setup(delay_shape,
  #                       delay_scale)

  delayfn <- function(n){
    x <- dist_setup(delay_shape, delay_scale, delay_off)(n)
    x <- pmax(x, 0.001)
    x/purrr::rbernoulli(n, detect_sen)
  }

  # calculate skewness from prop. pre-symptomatic transmission
  k <- presym2k(prop.presym)
  inf_fn <- purrr::partial(inf_fn_pre,  k = k, omga = si_omega)


  # Set initial values for loop indices
  total.cases <- num.initial.cases
  latest.onset <- 0
  max.gen <- 0
  extinct <- FALSE

  # Initial setup
  case_data <- outbreak_setup(num.initial.cases = num.initial.cases,
                            incfn = incfn,
                            prop.asym = prop.asym,
                            delayfn = delayfn)

  # Preallocate
  effective_r0_vect <- c()


  # Model loop
  while (latest.onset <= cap_max_days & total.cases < cap_cases &  max.gen < cap_gen & !extinct) {

    case_data <- outbreak_step(case_data = case_data,
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

    latest.onset <- min(case_data[generation == max(generation), exposure]) #max(case_data$onset)
    max.gen <- max(case_data$generation)
    total.cases <- nrow(case_data[infector > 0])
    extinct <- all(case_data$isolated)
  }
  case_data$Re[is.na(case_data$Re)] <- 0



  ### summarized info
  suminfo <- outbreak_summary(dat=case_data)

  sum_data <- cbind(
    data.frame(total.cases = total.cases, latest.onset = latest.onset,
               total_gen = max.gen, l_extinct = extinct),
    suminfo$sum_outbreak
    )

  if(!save.case_data) case_data <- list()

  return(list(case_data = case_data,
              sum_data = sum_data,
              index_data = suminfo$sum_cluster_index, isExtinct = extinct))

}