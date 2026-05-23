#' @rdname survfitJM
#' @inheritParams survfitjmcs
#' @export
survfitJM.jmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL,
                           u = NULL, Last.time = NULL, obs.time = NULL,
                           LOCF = FALSE, LOCFcovariate = NULL,
                           clongdata = NULL,
                           method = c("Laplace", "GH"),
                           quadpoint = NULL, ...) {
  
  survfitjmcs(
    object = object,
    seed = seed,
    ynewdata = ynewdata,
    cnewdata = cnewdata,
    u = u,
    Last.time = Last.time,
    obs.time = obs.time,
    LOCF = LOCF,
    LOCFcovariate = LOCFcovariate,
    clongdata = clongdata,
    method = method,
    quadpoint = quadpoint,
    ...
  )
}