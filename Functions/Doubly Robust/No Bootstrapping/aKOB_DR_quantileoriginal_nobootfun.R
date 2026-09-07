aKOB_DR_quantileoriginal_nobootfun <-
  function(data, survey_weights, outcome, prob) {
    out <- rlang::as_string(rlang::ensym(outcome))
    wt  <- rlang::as_string(rlang::ensym(survey_weights))
    
    D <- data
    D$union <- factor(D$union, levels = c("Nonunion", "Union"))
    
    u <- D$union
    u_numeric <- as.numeric(u) - 1
    w <- D[[wt]]
    y <- D[[out]]
    n <- nrow(D)
    p_quantile <- DescTools::Quantile(x = y,
                                      weights = w,
                                      probs = prob)
    y_binary <- y > p_quantile
    total_w <- sum(w)
    norm_w <- w / total_w
    density_outcome <- density(x = y, weights = norm_w)
    density_p_quantile <- approx(
      x = density_outcome$x,
      y = density_outcome$y,
      xout = p_quantile,
      rule = 2
    )$y
    
    idxU <- which(u == "Union")
    idxN <- which(u == "Nonunion")
    
    means <- c(Nonunion = weighted.mean(y_binary[idxN], w[idxN]),
               Union    = weighted.mean(y_binary[idxU], w[idxU]))
    gap <- (1 / density_p_quantile) * diff(means)
    
    props <- c(Nonunion = sum(w[idxN]) / total_w,
               Union    = sum(w[idxU]) / total_w)
    
    make_X <- function(rhs)
      model.matrix(reformulate(rhs), data = D)
    Xs_RI <- list(
      m3_1 = make_X(c(
        "union", "nonwhite", "education", "experience"
      )),
      m3_2 = make_X(c(
        "union", "nonwhite", "education", "marr"
      )),
      m3_3 = make_X(c(
        "union", "nonwhite", "experience", "marr"
      )),
      m3_4 = make_X(c(
        "union", "education", "experience", "marr"
      )),
      m4   = make_X(c(
        "union", "nonwhite", "education", "experience", "marr"
      ))
    )
    
    Xs_W <- list(
      m3_1 = make_X(c("nonwhite", "education", "experience")),
      m3_2 = make_X(c("nonwhite", "education", "marr")),
      m3_3 = make_X(c("nonwhite", "experience", "marr")),
      m3_4 = make_X(c("education", "experience", "marr")),
      m4   = make_X(c(
        "nonwhite", "education", "experience", "marr"
      ))
    )
    
    ucol <- sapply(Xs_RI, function(X) {
      j <- grep("^union", colnames(X))
      if (length(j) != 1L)
        stop("Couldn't uniquely identify the union dummy in model matrix.")
      j
    })
    Xs_cf <- mapply(function(X, j) {
      X2 <- X
      X2[idxU, j] <- 0
      X2
    }, Xs_RI, ucol, SIMPLIFY = FALSE)
    
    coefs_RI <- lapply(Xs_RI, function(X)
      fastglm::fastglm(
        x = X,
        y = as.numeric(y_binary),
        weights = norm_w,
        family = binomial()
      )$coef)
    
    predict_RI <- mapply(function(Xcf, b) {
      plogis(as.vector(Xcf %*% b))
    }, Xs_cf, coefs_RI, SIMPLIFY = FALSE)
    
    theta_RI <- mapply(function(mu) {
      weighted.mean(mu[idxU], w[idxU])
    }, predict_RI)
    
    coefs_W <- lapply(Xs_W, function(X)
      fastglm::fastglm(
        x = X,
        y = u_numeric,
        weights = norm_w,
        family = binomial()
      )$coef)
    
    est_probs_W <- mapply(function(X, b) {
      plogis(as.vector(X %*% b))
    }, Xs_W, coefs_W, SIMPLIFY = FALSE)
    
    weights_pscore_W <- lapply(est_probs_W, function(estprob)
      (1 - u_numeric) / props["Union"] * estprob / (1 - estprob))
    
    signal_theta <- mapply(function(weights_pscore_W, predict_RI, theta_RI) {
      weights_pscore_W * (y_binary - predict_RI) + (u_numeric / props["Union"]) * (predict_RI - theta_RI) + theta_RI
    },
    weights_pscore_W,
    predict_RI,
    theta_RI,
    SIMPLIFY = FALSE)
    
    cf_full <- (1 / density_p_quantile) * sapply(signal_theta, function(signal_theta)
      weighted.mean(x = signal_theta, w = w))
    
    explained <- (cf_full - means["Nonunion"] / density_p_quantile)
    unexplained <- (means["Union"] / density_p_quantile - cf_full)
    percentage <- explained / as.numeric(gap)
    
    nm <- names(Xs_RI)
    
    names(explained) <- names(unexplained) <- names(percentage) <- nm
    
    list(
      total_gap = list(means = means, gap = as.numeric(gap)),
      explained = explained,
      unexplained = unexplained,
      percentage = percentage,
      counterfactual = cf_full
    )
  }
