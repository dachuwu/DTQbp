#' Create partial function to sample from gamma distributions
#' @author Joel Hellewell
#' @param dist_shape numeric shape parameter of Weibull distribution
#' @param dist_scale numeric scale parameter of Weibull distribution
#'
#' @return partial function that takes a numeric argument for number of samples
#'
#' @importFrom purrr partial
#' @export
#'
dist_setup <- function(dist_shape = NULL, dist_scale = NULL, dist_off = 0){
  out <- function(n){
    pmax(rgamma(n = n, shape = dist_shape, scale = dist_scale) + dist_off, 0)
  }
    # purrr::partial(rgamma,
    #             shape = dist_shape,
    #             scale = dist_scale)
  return(out)
}


#' Samples the serial interval for given incubation period samples
#'
#' @param inc_samp vector of samples from the incubation period distribution
#' @param k numeric skew parameter for sampling the serial interval from the incubation period
#'
#' @return
#'
#' @importFrom sn rsn
#' @export
#'
inf_fn_pre <- function(inc_samp = NULL, k = NULL, omga = NULL) {

  out <- sn::rsn(n = length(inc_samp),
                 xi = inc_samp,
                 omega = omga,
                 alpha = k)

  out <- ifelse(out < 0.001, 0.001, out)

  return(out)
}



# Calculate k from pre-symptomatic infectiousness%
presym2k <- function(ia){
  ff <- function(k) (ia - sn::psn(0, xi = 0, omega = 2, alpha = k))^2
  tmp <- optim(c(1), ff, method = "BFGS")
  round(tmp$par,3)
}

# A vectorised version of isTRUE
vect_isTRUE <- function(x) {
  sapply(x, isTRUE)
}

vect_max <- function(x, y) {
  pmax.int(x, y)
  #purrr::map2_dbl(x, y, max)
}

vect_min <- function(x, y) {
  pmin.int(x, y)
  #purrr::map2_dbl(x, y, min)
}


# summarizing cluster information
outbreak_summary_cluster <- function(x0, x1, lasym, init.src=NULL){

  g <- igraph::graph_from_data_frame(data.frame(from=x0, to=x1))

  res <- lapply(1:length(init.src) , function(m){
    tmp <- igraph::bfs(g, as.character(init.src[m]),
                       neimode="out", unreachable = F, dist=T)
    ngen <- max(tmp$dist,na.rm = T)
      #max(tmp$dist[!is.na(tmp$dist)& tmp$dist>=0])

    tmp <- tmp$order[!is.na(tmp$order)]
    ncsize_all <- length(tmp)
    ncsize_sym <- sum(!lasym[as.integer(names(tmp))])

    return(data.frame(caseid = init.src[m], ngen = ngen,
                      ncsize_all = ncsize_all, ncsize_sym = ncsize_sym
    ))
  }) %>% data.table::rbindlist()

  return(res)
}




#' ks statistics
#' @export
#'
pf_ks.test <- function(x, y) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))

  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  STATISTIC <-  max(abs(z))

  return(STATISTIC)
}
