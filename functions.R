
# Method 0
homo_func <- function(n_0, n, x_0, x, s2_0, s2, direction = "upper") {
  k <- length(x)
  n_homo <- c(n_0, n)
  nu <- sum(n_homo - 1)
  s2_homo <- sum((n_homo - 1) * c(s2_0, s2)) / nu
  t <- (x - x_0) / sqrt(s2_homo) / sqrt(1 / n + 1 / n_0)
  cr <- matrix(0, nrow = k, ncol = k)
  for(i in 1:k) {
    for (j in 1:k) {
      cr[i, j] <- sqrt(n[i] * n[j] / (n[i] + n_0) / (n[j] + n_0))
    }
  }
  diag(cr) <- 1
  padj <- rep(NA, k)
  for (j in 1:k) {
    if (direction == "lower") {
      padj[j] <- 1 - pmvt(lower = rep(t[j], k), sigma = cr, df = nu)
    } else if (direction == "upper") {
      padj[j] <- 1 - pmvt(upper = rep(t[j], k), sigma = cr, df = nu)
    } else if (direction == "two.sided") {
      padj[j] <- 1 - pmvt(lower = rep(-abs(t[j]), k),
                          upper = rep(abs(t[j]), k), sigma = cr, df = nu)
    } else {
      stop("direction has to be one of 'lower', 'upper', 'two.sided'")
    }
  }
  return(padj)
}

# Method 1
ind_func <- function(n_0, n, x_0, x, s2_0, s2, direction = "upper") {
  k <- length(x)
  t <- (x - x_0) / sqrt(s2 / n + s2_0 / n_0)
  nu <- (s2 / n + s2_0 / n_0)^2 /
    (s2^2 / n^2 / (n - 1) + s2_0^2 / n_0^2 / (n_0 - 1))
  padj <- rep(NA, k)
  for (j in 1:k) {
    if (direction == "lower") {
      padj[j] <- 1 - prod(1 - pt(t[j], df = nu))
    } else if (direction == "upper") {
      padj[j] <- 1 - prod(pt(t[j], df = nu))
    } else if (direction == "two.sided") {
      padj[j] <- prod(1 - pt(abs(t[j]), df = nu) + pt(-abs(t[j]), df = nu))
    } else {
      stop("direction has to be one of 'lower', 'upper', 'two.sided'")
    }
  }
  return(padj)
}

# Method 2
PI_func <- function(n_0, n, x_0, x, s2_0, s2, direction = "upper") {
  k <- length(x)
  t <- (x - x_0) / sqrt(s2 / n + s2_0 / n_0)
  nu <- (s2 / n + s2_0 / n_0)^2 /
    (s2^2 / n^2 / (n - 1) + s2_0^2 / n_0^2 / (n_0 - 1))
  lambda <- sqrt(s2_0 / n_0 / (s2_0 / n_0 + s2 / n))
  cr <- matrix(0, nrow = k, ncol = k)
  for(i in 1:k) {
    for (j in 1:k) {
      cr[i, j] <- lambda[i] * lambda[j]
    }
  }
  diag(cr) <- 1
  padj <- rep(0, k)
  for (j in 1:k) {
    if (direction == "lower") {
      padj[j] <- 1 - pmvt(lower = rep(t[j], k), sigma = cr, df = round(nu[j]))
    } else if (direction == "upper") {
      padj[j] <- 1 - pmvt(upper = rep(t[j], k), sigma = cr, df = round(nu[j]))
    } else if (direction == "two.sided") {
      padj[j] <- 1 - pmvt(lower = rep(-abs(t[j]), k),
                          upper = rep(abs(t[j]), k), sigma = cr, df = round(nu[j]))
    } else {
      stop("direction has to be one of 'lower', 'upper', 'two.sided'")
    }
  }
  return(padj)
}

# Method 3
sim_based_func <- function(n_0, n, x_0, x, s2_0, s2, direction, nsim) {
  k <- length(x)
  t <- (x - x_0) / sqrt(s2 / n + s2_0 / n_0)
  nu <- (s2 / n + s2_0 / n_0)^2 /
    (s2^2 / n^2 / (n - 1) + s2_0^2 / n_0^2 / (n_0 - 1))
  lambda <- sqrt(s2_0 / n_0 / (s2_0 / n_0 + s2 / n))
  cr <- matrix(0, nrow = k, ncol = k)
  for(i in 1:k) {
    for (j in 1:k) {
      cr[i, j] <- lambda[i] * lambda[j]
    }
  }
  diag(cr) <- 1
  eta <- s2_0 / n_0 / sqrt(n_0 - 1) /
    sqrt((s2_0^2 / n_0^2 / (n_0 - 1) + s2^2 / n^2 / (n - 1)))
  c_val <- eta^2 * nu
  d_val <- nu - c_val
  ranking <- rank(c_val, ties.method = "first")
  v_rank <- matrix(0, nrow = nsim, ncol = k)
  v_rank[, 1] <- rchisq(nsim, df = c_val[ranking == 1])
  for (j in 2:k) {
    v_rank[, j] <- v_rank[, j - 1] +
      rchisq(nsim, df = c_val[ranking == j] - c_val[ranking == (j - 1)])
  }
  v <- v_rank
  for (j in 1:k) {
    v[, j] <- v_rank[, ranking[j]]
  }
  w <- matrix(0, nrow = nsim, ncol = k)
  for (j in 1:k) {
    w[, j] <- rchisq(nsim, df = d_val[j])
  }
  Z <- rmvnorm(nsim, mean = rep(0, k), sigma = cr)
  tt <- Z / sqrt(v + w) * matrix(rep(sqrt(nu), each = nsim), nrow = nsim)
  padj <- rep(NA, k)
  for (j in 1:k) {
    if (direction == "lower") {
      padj[j] <- mean(apply(tt, 1, min) <= t[j])
    } else if (direction == "upper") {
      padj[j] <- mean(apply(tt, 1, max) >= t[j])
    } else if (direction == "two.sided") {
      padj[j] <- mean(apply(tt, 1, max) >= abs(t[j])
                      | apply(tt, 1, min) <= -abs(t[j]))
    } else {
      stop("direction has to be one of 'lower', 'upper', 'two.sided'")
    }
  }
  return(padj)
}

