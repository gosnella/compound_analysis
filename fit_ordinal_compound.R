fit_ordinal_compound <- function (formula, compounds, fingerprints, data, weights, subset,
          offset, flink = c("logit", "probit", "Gosset", "Wallace",
                            "boxcox", "boxcoxmod", "GEV", "GEVmod", "GEVmodNS", "gennorm"),
          linkpar = 0, fcor = c("matern", "tanimoto", "powerexponential",
                                "exponential", "gaussian", "independent"), corpar = 0,
          startingvals = list(), n_int_pnts = 10L, nlm_fscale = 1,
          nlm_typsize = 1, nlm_gradtol = 1e-06, nlm_steptol = 1e-06,
          returnData = TRUE)
{
  cl <- match.call()
  flink <- match.arg(flink)
  ilink <- match(flink, eval(formals()$flink))
  fcor <- match.arg(fcor)
  icor <- match(fcor, eval(formals()$fcor))
  nophi <- fcor %in% c("tanimoto", "independent")
  n_int_pnts <- as.integer(n_int_pnts + n_int_pnts + 1L)
  if (n_int_pnts <= 0)
    stop("Need positive n_int_points.")
  if (missing(data))
    data <- environment(formula)
  if (length(formula) != 3)
    stop("The formula input is incomplete.")
  formula <- stats::update(formula, ~. + 1)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$formula <- formula
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  DesignFixed <- delete.intercept(model.matrix(mt, mf))
  l <- ncol(DesignFixed)
  yy <- model.response(mf)
  if (!is.factor(yy)) {
    stop("The response must be a factor.")
  }
  n <- length(yy)
  ynm <- levels(yy)
  cm1 <- nlevels(yy) - 1L
  wghts <- model.weights(mf)
  if (!is.null(wghts)) {
    if (length(wghts) != n) {
      stop(gettextf("Number of weights is %d, should equal %d (number of observations)",
                    length(wghts), n), domain = NA)
    }
    else {
      wghts <- as.double(wghts)
    }
  }
  else {
    wghts <- rep(1, n)
  }
  oofset <- as.vector(model.offset(mf))
  if (!is.null(oofset)) {
    if (length(oofset) != n) {
      stop(gettextf("Number of offsets is %d, should equal %d (number of observations)",
                    length(oofset), n), domain = NA)
    }
    else {
      oofset <- as.double(oofset)
    }
  }
  else {
    oofset <- double(n)
  }
  compounds <- stats::update(compounds, NULL ~ . - 1)
  mfatc <- mfc
  mfatc$weights <- mfatc$offset <- NULL
  mfatc$formula <- compounds
  mfatc$drop.unused.levels <- FALSE
  mfat <- eval(mfatc, parent.frame())
  mt <- attr(mfat, "terms")
  mfat <- set_contrasts(mt, mfat)
  DesignGP <- model.matrix(mt, mfat)
  is_y_finite <- is.finite(yy)
  n_finite <- sum(is_y_finite)
  y <- as.integer(yy) - 1L
  DistMat <- compounds::compound_euclidean_dist(fingerprints)
  m <- nrow(fingerprints)
  if (is.null(startingvals[["u"]])) {
    startingvals[["u"]] <- ustart <- rep(0, m)
  }
  else {
    ustart <- as.double(startingvals[["u"]])
    if (length(ustart) != m) {
      stop("The length of starting values for u is not equal to the number of fingerprints.")
    }
  }
  npar <- cm1 + l + 1 + (!nophi)
  par <- numeric(cm1 + l + 2)
  if (is.null(startingvals$alpha)) {
    startingvals$alpha <- alpha <- as.double(seq(-1, 1, length = cm1))
  }
  else {
    alpha <- as.double(startingvals$alpha)
    if (length(alpha) != cm1)
      stop("Wrong length of parameter alpha.")
  }
  par[1] <- alpha[1]
  if (cm1 > 1)
    par[2:cm1] <- log(alpha[-1] - alpha[-cm1])
  if (is.null(startingvals$beta)) {
    startingvals$beta <- beta <- double(l)
  }
  else {
    beta <- as.double(startingvals$beta)
    if (length(beta) != l)
      stop("Wrong length of parameter beta.")
  }
  par[(cm1 + 1):(cm1 + l)] <- beta
  if (is.null(startingvals$ssq)) {
    startingvals$ssq <- ssq <- 1
  }
  else {
    ssq <- as.double(startingvals$ssq)
    if (length(ssq) != 1)
      stop("Wrong length of parameter ssq.")
  }
  par[cm1 + l + 1] <- log(ssq)
  if (nophi || is.null(startingvals$phi)) {
    startingvals$phi <- phi <- 1
  }
  else {
    phi <- as.double(startingvals$phi)
    if (length(phi) != 1)
      stop("Wrong length of parameter phi.")
  }
  if (!nophi) {
    par[cm1 + l + 2] <- log(phi)
  }
  nlm_typsize <- rep(nlm_typsize, length = npar)
  ll <- 0
  ll_dalpha <- rep(0, cm1)
  ll_dbeta <- rep(0, l)
  ll_dssq <- ll_dphi <- 0
  u_grad <- rep(0, m)
  u_hess <- u_hessinv <- matrix(0, m, m)
  op_time <- system.time({
    f90 <- .Fortran("fit_model_compounds", ll, ll_dalpha,
                    ll_dbeta, ll_dssq, ll_dphi, phi, DistMat, y[is_y_finite],
                    alpha, DesignFixed[is_y_finite, ], beta, DesignGP[is_y_finite,
                    ], ssq, oofset[is_y_finite], ustart, u_grad,
                    u_hess, u_hessinv, cm1, l, m, n_finite, linkpar,
                    ilink, corpar, icor)
  })
  alpha <- f90[[9]]
  beta <- f90[[11]]
  ssq <- f90[[13]]
  phi <- f90[[6]]
  par[1] <- alpha[1]
  if (cm1 > 1)
    par[2:cm1] <- pmax(-15, log(alpha[-1] - alpha[-cm1]))
  par[(cm1 + 1):(cm1 + l)] <- beta
  par[cm1 + l + 1] <- max(-15, log(ssq))
  if (!nophi) {
    par[cm1 + l + 2] <- max(-15, log(phi))
  }
  op <- nlm(nll_compound, par[1:npar], DamageCat = yy[is_y_finite],
            DesignFixed = DesignFixed[is_y_finite, ], DesignGP = DesignGP[is_y_finite,
            ], distmat = DistMat, ustart = ustart, flink = flink,
            linkpar = linkpar, fcor = fcor, corpar = corpar, offset = oofset[is_y_finite],
            nophi = nophi, hessian = TRUE, check.analyticals = FALSE,
            fscale = nlm_fscale, typsize = nlm_typsize, gradtol = nlm_gradtol,
            steptol = nlm_steptol)
  out <- list()
  out$optimisation <- op
  out$optimisation$time <- op_time
  out$parameters <- list()
  par[1:npar] <- op$estimate
  alpha <- par[1:cm1]
  if (cm1 > 1) {
    alpha[2:cm1] <- exp(alpha[2:cm1])
    alpha <- cumsum(alpha)
  }
  names(alpha) <- paste0(ynm[1L:cm1], "|", ynm[2L:(cm1 + 1L)])
  out$parameters$alpha <- alpha
  beta <- par[(cm1 + 1):(cm1 + l)]
  names(beta) <- colnames(DesignFixed)
  out$parameters$beta <- beta
  ssq <- exp(par[cm1 + l + 1])
  out$parameters$ssq <- ssq
  if (!nophi)
    phi <- exp(par[cm1 + l + 2])
  out$parameters$phi <- if (nophi)
    NA
  else phi
  nll <- nll_compound(par[1:npar], yy[is_y_finite], DesignFixed[is_y_finite,
  ], DesignGP[is_y_finite, ], DistMat, ustart, flink, linkpar,
  fcor, corpar, oofset[is_y_finite], nophi)
  out$negloglik <- op$minimum
  par_hessian <- matrix(0, cm1 + l + 2, cm1 + l + 2)
  par_hessian[1:npar, 1:npar] <- op$hessian
  par_hessian[cm1 + l + 1, cm1 + l + 1] <- par_hessian[cm1 +
                                                         l + 1, cm1 + l + 1] * ssq
  if (!nophi)
    par_hessian[cm1 + l + 2, cm1 + l + 2] <- par_hessian[cm1 +
                                                           l + 2, cm1 + l + 2] * phi
  if (cm1 > 1) {
    for (i in 2:cm1) par_hessian[i, i] <- par_hessian[i,
                                                      i] * exp(par[i])
    Jmat <- matrix(0, cm1, cm1)
    Jmat[lower.tri(Jmat, diag = TRUE)] <- 1
    par_hessian[1:cm1, 1:cm1] <- Jmat %*% par_hessian[1:cm1,
                                                      1:cm1] %*% t(Jmat)
  }
  colnames(par_hessian) <- rownames(par_hessian) <- c(paste0("alpha_",
                                                             1:cm1), paste0("beta_", 1:l), "ssq", "phi")
  out$parameters$hessian <- par_hessian
  out$GP <- list()
  uhat <- attr(nll, "uhat")
  out$GP$u <- uhat
  out$GP$gradient <- attr(nll, "u_grad")
  out$GP$hessian <- attr(nll, "u_hess")
  uVar <- attr(nll, "u_hessinv")
  out$GP$hessian_inv <- uVar
  ww <- tt <- double(n_int_pnts)
  f90b <- .Fortran("hermite_rule", n_int_pnts, tt, ww)
  ww <- f90b[[3]]/sqrt(2 * pi)
  tt <- f90b[[2]]
  uu <- uhat + sqrt(diag(uVar)) %o% tt
  etaa <- c(oofset) + c(DesignFixed %*% beta) + DesignGP %*%
    uu
  eta <- outer(alpha, etaa, "+")
  Geta <- linkinv(eta, flink, linkpar)
  int_Geta <- matrix(matrix(Geta, ncol = n_int_pnts) %*% ww,
                     cm1)
  pred_prob <- rbind(int_Geta, 1)
  pred_prob[2:(cm1 + 1), ] <- pred_prob[2:(cm1 + 1), ] - pred_prob[1:cm1,
  ]
  out$prediction <- pred_prob
  out$AICC <- 2 * (op$minimum + n * length(op$estimate)/(n -
                                                           1 - length(op$estimate)))
  out$call <- cl
  out$Data <- list()
  if (isTRUE(returnData)) {
    out$Data$Response <- yy
    out$Data$DesignFixed <- DesignFixed
    out$Data$DesignGP <- DesignGP
    out$Data$Fingerprints <- fingerprints
    out$Data$DistMat <- DistMat
    out$Data$offset <- oofset
    out$Data$Finite <- is_y_finite
  }
  else {
    out$Data$Response <- NULL
    out$Data$DesignFixed <- NULL
    out$Data$DesignGP <- NULL
    out$Data$Fingerprints <- NULL
    out$Data$DistMat <- NULL
    out$Data$offset <- NULL
    out$Data$Finite <- NULL
  }
  out$Model <- list()
  out$Model$formula_fixed <- formula
  out$Model$formula_compounds <- compounds
  out$Model$flink <- flink
  out$Model$ilink <- ilink
  out$Model$linkpar <- linkpar
  out$Model$fcor <- fcor
  out$Model$icor <- icor
  out$Model$corpar <- corpar
  out$Model$nophi <- nophi
  out$Fit <- list()
  out$Fit$startingvals <- startingvals
  out$Fit$n_int_pnts <- n_int_pnts
  out
}