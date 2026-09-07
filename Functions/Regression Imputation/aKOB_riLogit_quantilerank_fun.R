aKOB_riLogit_quantilerank_fun <-
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
    
    make_X <- function(rhs)
      model.matrix(reformulate(rhs), data = D)
    Xs <- list(
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
    
    ucol <- sapply(Xs, function(X) {
      j <- grep("^union", colnames(X))
      if (length(j) != 1L)
        stop("Couldn't uniquely identify the union dummy in model matrix.")
      j
    })
    Xs_cf <- mapply(function(X, j) {
      X2 <- X
      X2[idxU, j] <- 0
      X2
    }, Xs, ucol, SIMPLIFY = FALSE)
    
    coefs <- lapply(Xs, function(X)
      fastglm::fastglm(
        x = X,
        y = as.numeric(y_binary),
        weights = norm_w,
        family = binomial()
      )$coef)
    
    cf_full <- mapply(function(Xcf, b) {
      mu <- plogis(as.vector(Xcf[idxU, , drop = FALSE] %*% b))
      weighted.mean(mu, w[idxU])
    }, Xs_cf, coefs)
    
    explained <- (cf_full - means["Nonunion"])
    unexplained <- (means["Union"] - cf_full)
    percentage <- explained / as.numeric(gap)
    
    one_boot <- function() {
      idx <- sample.int(n, replace = TRUE)
      u_b <- u[idx]
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
      
      X_idx   <- lapply(Xs, function(X)
        X[idx, , drop = FALSE])
      Xcf_U_b <- lapply(Xs_cf, function(Xcf)
        Xcf[idx[idxU_b], , drop = FALSE])
      
      betas_b <- lapply(X_idx, function(X)
        fastglm::fastglm(
          x = X,
          y = as.numeric(y_binary_b),
          weights = norm_w_b,
          family = binomial()
        )$coef)
      
      cf_b <- mapply(function(Xcf, b) {
        mu <- plogis(as.vector(Xcf %*% b))
        weighted.mean(mu, w_b[idxU_b])
      }, Xcf_U_b, betas_b)
      
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
