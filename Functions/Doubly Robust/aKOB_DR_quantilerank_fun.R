aKOB_DR_quantilerank_fun <-
  function(data,
           survey_weights,
           outcome,
           prob,
           bootstrap_iters = 50) {
    out <- rlang::as_string(rlang::ensym(outcome))
    wt  <- rlang::as_string(rlang::ensym(survey_weights))
    
    D <- data
    D$union <- factor(D$union, levels = c("Nonunion", "Union"))
    
    u <- D$union
    u_numeric <- as.numeric(u) - 1
    w <- D[[wt]]
    y <- D[[out]]
    l <- D$lwage
    n <- nrow(D)
    total_w <- sum(w)
    norm_w <- w / total_w
    y_binary <- y > prob
    
    idxU <- which(u == "Union")
    idxN <- which(u == "Nonunion")
    
    means <- c(Nonunion = weighted.mean(y_binary[idxN], w[idxN]),
               Union    = weighted.mean(y_binary[idxU], w[idxU]))
    gap <- diff(means)
    
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
    
    cf_full <- sapply(signal_theta, function(signal_theta)
      weighted.mean(x = signal_theta, w = w))
    
    explained <- (cf_full - means["Nonunion"])
    unexplained <- (means["Union"] - cf_full)
    percentage <- explained / as.numeric(gap)
    
    one_boot <- function() {
      idx <- sample.int(n, replace = TRUE)
      u_b <- u[idx]
      u_numeric_b <- as.numeric(u_b) - 1
      w_b <- w[idx]
      l_b <- l[idx]
      total_w_b <- sum(w_b)
      norm_w_b <- w_b / total_w_b
      y_b <- wtd.rank(signif(l_b, 10), weights = w_b) / total_w_b
      y_binary_b <- y_b > prob
      
      idxU_b <- which(u_b == "Union")
      idxN_b <- which(u_b == "Nonunion")
      
      means_b <- c(Nonunion = weighted.mean(y_binary_b[idxN_b], w_b[idxN_b]),
                   Union    = weighted.mean(y_binary_b[idxU_b], w_b[idxU_b]))
      
      props_b <- c(Nonunion = sum(w_b[idxN_b]) / total_w_b,
                   Union    = sum(w_b[idxU_b]) / total_w_b)
      
      Xs_RI_idx   <- lapply(Xs_RI, function(X)
        X[idx, , drop = FALSE])
      Xs_W_idx   <- lapply(Xs_W, function(X)
        X[idx, , drop = FALSE])
      
      Xcf_b <- lapply(Xs_cf, function(Xcf)
        Xcf[idx, , drop = FALSE])
      
      coefs_RI_b <- lapply(Xs_RI_idx, function(X)
        fastglm::fastglm(
          x = X,
          y = as.numeric(y_binary_b),
          weights = norm_w_b,
          family = binomial()
        )$coef)
      
      predict_RI_b <- mapply(function(Xcf, b) {
        plogis(as.vector(Xcf %*% b))
      }, Xcf_b, coefs_RI_b, SIMPLIFY = FALSE)
      
      theta_RI_b <- mapply(function(mu)
        weighted.mean(mu[idxU_b], w_b[idxU_b]), predict_RI_b)
      
      coefs_W_b <- lapply(Xs_W_idx, function(X)
        fastglm::fastglm(
          x = X,
          y = u_numeric_b,
          weights = norm_w_b,
          family = binomial()
        )$coef)
      
      est_probs_W_b <- mapply(function(X, b) {
        plogis(as.vector(X %*% b))
      }, Xs_W_idx, coefs_W_b, SIMPLIFY = FALSE)
      
      weights_pscore_W_b <- lapply(est_probs_W_b, function(estprob)
        (1 - u_numeric_b) / props_b["Union"] * estprob / (1 - estprob))
      
      signal_theta_b <- mapply(function(weights_pscore_W,
                                        predict_RI,
                                        theta_RI) {
        weights_pscore_W * (y_binary_b - predict_RI) + (u_numeric_b / props_b["Union"]) * (predict_RI - theta_RI) + theta_RI
      },
      weights_pscore_W_b,
      predict_RI_b,
      theta_RI_b,
      SIMPLIFY = FALSE)
      
      cf_b <- sapply(signal_theta_b, function(signal_theta)
        weighted.mean(x = signal_theta, w = w_b))
      
      c(explained   = (cf_b - means_b["Nonunion"]),
        unexplained = (means_b["Union"] - cf_b))
    }
    
    boot_res <- replicate(bootstrap_iters, one_boot(), simplify = "array")
    
    ex_boot <- t(boot_res[1:5, , drop = FALSE])
    un_boot <- t(boot_res[6:10, , drop = FALSE])
    
    se_ex <- apply(ex_boot, 2, sd, na.rm = TRUE)
    se_un <- apply(un_boot, 2, sd, na.rm = TRUE)
    
    nm <- names(Xs_RI)
    
    names(explained) <- names(unexplained) <- names(percentage) <- names(se_ex) <- names(se_un) <- nm
    
    list(
      total_gap = list(means = means, gap = as.numeric(gap)),
      explained = explained,
      unexplained = unexplained,
      percentage = percentage,
      counterfactual = cf_full,
      explained_se = se_ex,
      unexplained_se = se_un,
      boot_res = boot_res
    )
  }
