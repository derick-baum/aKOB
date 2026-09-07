ABO_union_wLogit_quantilerank_fast <- function(data,
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
  y_binary <- y > prob
  total_w <- sum(w)
  norm_w <- w / total_w
  
  idxU <- which(u == "Union")
  idxN <- which(u == "Nonunion")
  
  means <- c(Nonunion = weighted.mean(y_binary[idxN], w[idxN]),
             Union    = weighted.mean(y_binary[idxU], w[idxU]))
  gap <- diff(means)
  
  props <- c(Nonunion = sum(w[idxN]) / total_w,
             Union    = sum(w[idxU]) / total_w)
  
  make_X <- function(rhs)
    model.matrix(reformulate(rhs), data = D)
  Xs <- list(
    m3_1 = make_X(c("nonwhite", "education", "experience")),
    m3_2 = make_X(c("nonwhite", "education", "marr")),
    m3_3 = make_X(c("nonwhite", "experience", "marr")),
    m3_4 = make_X(c("education", "experience", "marr")),
    m4   = make_X(c(
      "nonwhite", "education", "experience", "marr"
    ))
  )
  
  coefs <- lapply(Xs, function(X)
    fastglm::fastglm(
      x = X,
      y = u_numeric,
      weights = norm_w,
      family = binomial()
    )$coef)
  
  est_probs <- mapply(function(X, b) {
    plogis(as.vector(X %*% b))
  }, Xs, coefs, SIMPLIFY = FALSE)
  
  weights_pscore <- lapply(est_probs, function(estprob)
    (1 - u_numeric) / props["Union"] * estprob / (1 - estprob))
  
  cf_full <- sapply(weights_pscore, function(weightspscore)
    weighted.mean(x = y_binary, w = w * weightspscore))
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
    
    X_idx   <- lapply(Xs, function(X)
      X[idx, , drop = FALSE])
    
    betas_b <- lapply(X_idx, function(X)
      fastglm::fastglm(
        x = X,
        y = u_numeric_b,
        weights = norm_w_b,
        family = binomial()
      )$coef)
    
    est_probs_b <- mapply(function(X, b) {
      plogis(as.vector(X %*% b))
    }, X_idx, betas_b, SIMPLIFY = FALSE)
    
    weights_pscore_b <- lapply(est_probs_b, function(estprob)
      (1 - u_numeric_b) / props_b["Union"] * estprob / (1 - estprob))
    
    cf_b <- sapply(weights_pscore_b, function(weightspscore)
      weighted.mean(x = y_binary_b, w = w_b * weightspscore))
    
    c(explained   = (cf_b - means_b["Nonunion"]),
      unexplained = (means_b["Union"] - cf_b))
  }
  
  boot_res <- replicate(bootstrap_iters, one_boot(), simplify = "array")
  
  ex_boot <- t(boot_res[1:5, , drop = FALSE])
  un_boot <- t(boot_res[6:10, , drop = FALSE])
  
  se_ex <- apply(ex_boot, 2, sd, na.rm = TRUE)
  se_un <- apply(un_boot, 2, sd, na.rm = TRUE)
  
  names(explained) <- names(unexplained) <- names(percentage) <- names(se_ex) <- names(se_un) <- names(Xs)
  
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
