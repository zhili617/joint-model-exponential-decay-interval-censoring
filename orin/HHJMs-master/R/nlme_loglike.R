nlme_loglike_EXP_Delay <- function(nlmeObject){
  
  vars <- fmReverse(nlmeObject$fm)
  resp <- vars$resp
  
  A_names   <- nlmeObject$A_names
  ind_names <- nlmeObject$ind_names
  lag_names <- nlmeObject$lag_names
  lambda_nm <- nlmeObject$lambda_name
  ran_map   <- nlmeObject$ran_map                  # list(lambda = "b_11")
  re_disp   <- nlmeObject$re_disp                  # list(b_11 = "sigma_b11")
  lambda_scale <- if (is.null(nlmeObject$lambda_scale)) "log" else nlmeObject$lambda_scale
  lambda_div   <- if (is.null(nlmeObject$lambda_div)) 1 else nlmeObject$lambda_div
  
  # parameter + (scale * random effect)
  add_re <- function(par_nm){
    re <- ran_map[[par_nm]]
    if (is.null(re)) return(par_nm)
    sc <- if (!is.null(re_disp)) re_disp[[re]] else NULL
    if (is.null(sc)) paste0("(", par_nm, "+", re, ")")
    else             paste0("(", par_nm, "+", sc, "*", re, ")")
  }
  
  # decay rate: log scale exp(lambda + sc*b), or original scale (lambda + sc*b)/div
  lambda_i <- if (lambda_scale == "log") {
    paste0("exp(", add_re(lambda_nm), ")")
  } else {
    paste0("(", add_re(lambda_nm), ")/", lambda_div)
  }
  
  A_i <- vapply(A_names, add_re, FUN.VALUE = character(1))
  
  terms <- character(length(A_names))
  for (k in seq_along(A_names)) {
    terms[k] <- paste0(A_i[k], "*", ind_names[k],
                       "*exp(-", lambda_i, "*", lag_names[k], ")")
  }
  mu <- paste(terms, collapse = "+")
  
  sigma <- nlmeObject$sigma
  loglike <- paste0(
    "-0.5*(", resp, "-(", mu, "))^2/", sigma, "^2",
    "-log(", sigma, ")-0.5*log(2*pi)"
  )
  
  # scale parameters listed as 'fixed'; JMfit_HL moves them to dispersion via nlmeObject$disp
  scale_pars <- if (!is.null(re_disp)) unlist(re_disp, use.names = FALSE) else NULL
  fixed  <- c(A_names, lambda_nm, scale_pars)
  raneff <- if (length(ran_map) > 0) unlist(ran_map, use.names = FALSE) else NULL
  
  list(
    loglike = loglike, fixed = fixed, raneff = raneff,
    mu = mu, linear_pred = mu,
    rvX = c(ind_names, lag_names), rvZ = names(ran_map),
    resp = resp, sigma = sigma
  )
}