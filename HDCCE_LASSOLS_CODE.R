# HDCCE_HAC_bandwidth.R
# HDCCE Results and robustness checks using Lasso/Least-squares estimators
# Author: Minsoo Lee  (supervised by Dr. Jessica Leung & Prof. Param Silvapulle)

# Import data
data <- openxlsx::read.xlsx("Final_Cabon_Emm_disaggregate_zscored.xlsx", colNames = TRUE, detectDates = TRUE)
p <- 22 # drivers
T_dim <- 30 #1990-2019

z_cols        <- grep("_z$", names(data), value = TRUE)
driver_cols_z <- setdiff(z_cols, "carbon_emissions_pc_z")
stopifnot(length(driver_cols_z) == p)

year_vec  <- data$Year
uniq_yrs  <- sort(unique(year_vec))
T_window  <- uniq_yrs[(length(uniq_yrs) - (T_dim - 1)):length(uniq_yrs)]
data_9019 <- data[data$Year %in% T_window, ]
data_9019 <- data_9019[order(data_9019$Country.Code), ]
data_9019 <- data_9019[(data_9019$Country.Code != "CRI"), ]

data_all    <- data_9019
data_tax    <- data_9019[data_9019$Carbon.Tax == 1, ]
data_no_tax <- data_9019[data_9019$Carbon.Tax == 0, ]
groups_list <- list("All" = data_all, "Carbon tax" = data_tax, "No carbon tax" = data_no_tax)
GROUP_ORDER <- names(groups_list)

## ------------------------------------------------------------------------
## Settings
## ------------------------------------------------------------------------
K_ALPHA_THRESHOLD <- 0.01
DK_LAG      <- 2                     # Bartlett lag window for cross-country pairs
KAPPA_MULT  <- 1.25                  # node-wise penalty rule constant (unchanged)
HAC_KERNELS <- c("uniform", "bartlett")
T_eff       <- T_dim - 1             # usable periods per country (29)
# Newey-West rounded down to whole lags (3 at T=29); T^(1/3) left unrounded (3.07)
HT_ROUND_NW  <- TRUE
HT_ROUND_T13 <- FALSE
rnd_nw       <- if(HT_ROUND_NW)  floor else identity
rnd_t13      <- if(HT_ROUND_T13) floor else identity
EXCEL_FILE  <- "./tables/HDCCE_LASSOLS_RESULTS.xlsx"
dir.create("./tables", showWarnings = FALSE)

# Bandwidth rules -> unique h_T values with a label naming the rule(s)
ht_rules <- data.frame(
  h_T  = c(T_eff,
           rnd_nw(4 * (T_eff / 100)^(2/9)),
           rnd_t13(T_eff^(1/3))),
  rule = c("h_T = T (Linton baseline)",
           if(HT_ROUND_NW)  "Newey-West (1994): floor(4(T/100)^(2/9))" else "Newey-West (1994): 4(T/100)^(2/9), unrounded",
           if(HT_ROUND_T13) "T^(1/3) rate: floor(T^(1/3))"              else "T^(1/3) rate: T^(1/3), unrounded"),
  stringsAsFactors = FALSE
)

ht_grid <- aggregate(rule ~ h_T, data = ht_rules,
                     FUN = function(x) paste(x, collapse = "; "))
print(ht_grid)

# One row per variance variant: SE type x window shape x bandwidth
variants <- expand.grid(se_type = c("cluster", "DK_cluster"),
                        kernel  = HAC_KERNELS,
                        h_T     = ht_grid$h_T,
                        stringsAsFactors = FALSE)
variants$rule <- ht_grid$rule[match(variants$h_T, ht_grid$h_T)]
variants$id   <- paste(variants$se_type, variants$kernel, variants$h_T, sep = "|")
n_var <- nrow(variants)

## ------------------------------------------------------------------------
## Variance helpers
## ------------------------------------------------------------------------
## score = x_it*e_it (LS) or node-residual*lasso-residual (lasso).
## cluster: within-country only. DK_cluster: adds a cross-country Bartlett_m term.

# Weight at lag j >= 1 for window h
kern_w <- function(j, h, kernel){
  if(kernel == "uniform") as.numeric(j <= h) else pmax(1 - j / (h + 1), 0) * (j <= h)
}

# Long-run variance of a Tn x R matrix H (row t = score at date t)
lrv_kernel <- function(H, h, kernel){
  Tn <- nrow(H)
  S  <- crossprod(H)
  # Bartlett weight at h_T = T is 1 - j/(T+1)
  if(h > 0) for(j in 1:min(h, Tn - 1)){
    w  <- kern_w(j, h, kernel)
    Gj <- crossprod(H[(j+1):Tn, , drop = FALSE], H[1:(Tn-j), , drop = FALSE])
    S  <- S + w * (Gj + t(Gj))
  }
  S
}

# Scalar version: sum of the columns' long-run variances
lrv_scalar_kernel <- function(G, h, kernel){
  Tn  <- nrow(G)
  out <- sum(G^2)
  if(h > 0) for(j in 1:min(h, Tn - 1)){
    out <- out + 2 * kern_w(j, h, kernel) *
      sum(G[(j+1):Tn, , drop = FALSE] * G[1:(Tn-j), , drop = FALSE])
  }
  out
}

sum_over_countries <- function(score, n_sub, Tn){
  R <- ncol(score)
  matrix(apply(array(score, dim = c(Tn, n_sub, R)), c(1, 3), sum),
         nrow = Tn, ncol = R)
}

# Clip negative eigenvalues (combined windows aren't guaranteed PSD)
psd_truncate <- function(M){
  M  <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
}

# R x R score covariance for every variant (LS). Returns a list named by variants$id.
score_cov_variants <- function(score, n_sub, Tn, m, variants){
  blocks <- lapply(seq_len(n_sub),
                   function(i) score[((i-1)*Tn + 1):(i*Tn), , drop = FALSE])
  S_dk       <- lrv_kernel(sum_over_countries(score, n_sub, Tn), m, "bartlett")
  S_within_m <- Reduce(`+`, lapply(blocks, lrv_kernel, h = m, kernel = "bartlett"))
  cache <- list()
  out   <- setNames(vector("list", nrow(variants)), variants$id)
  for(r in seq_len(nrow(variants))){
    key <- paste(variants$kernel[r], variants$h_T[r])
    if(is.null(cache[[key]])){
      cache[[key]] <- Reduce(`+`, lapply(blocks, lrv_kernel,
                                         h = variants$h_T[r], kernel = variants$kernel[r]))
    }
    S_within_h <- cache[[key]]
    out[[r]] <- if(variants$se_type[r] == "cluster") psd_truncate(S_within_h)
    else psd_truncate(S_dk + S_within_h - S_within_m)
  }
  out
}

# Scalar score variance for every variant (lasso); G is node-residual * lasso-residual
score_var_variants <- function(G, m, variants){
  S_dk       <- lrv_scalar_kernel(matrix(rowSums(G), ncol = 1), m, "bartlett")
  S_within_m <- lrv_scalar_kernel(G, m, "bartlett")
  cache <- list()
  vals  <- numeric(nrow(variants))
  for(r in seq_len(nrow(variants))){
    key <- paste(variants$kernel[r], variants$h_T[r])
    if(is.null(cache[[key]])){
      cache[[key]] <- lrv_scalar_kernel(G, variants$h_T[r], variants$kernel[r])
    }
    S_within_h <- cache[[key]]
    vals[r] <- if(variants$se_type[r] == "cluster") max(S_within_h, 0)
    else max(S_dk + S_within_h - S_within_m, 0)
  }
  vals
}

# Significance stars from p x length(alpha) CI matrices, alpha = c(0.01, 0.05, 0.1)
make_stars <- function(conf_min, conf_max){
  n_alpha <- ncol(conf_min)
  codes <- numeric(nrow(conf_min))
  for(j in seq_len(nrow(conf_min))){
    for(a in seq_len(n_alpha)){
      if(isTRUE(conf_min[j,a] < 0 && 0 < conf_max[j,a])) codes[j] <- a
    }
  }
  vapply(codes, function(cd) paste(rep("*", n_alpha - cd), collapse = ""), character(1))
}

## ------------------------------------------------------------------------
## HD-CCE projection (identical construction to the main script)
## ------------------------------------------------------------------------
prep_hdcce <- function(data_sub, alpha_threshold = K_ALPHA_THRESHOLD){
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  Y_sub <- data_sub$carbon_emissions_pc
  X_sub_var <- apply(X_sub, MARGIN = 2, FUN = var)
  
  X_bar <- matrix(NA, ncol = p, nrow = T_dim)
  for(t in 1:T_dim){
    indices <- seq(t, n_sub * T_dim, by = T_dim)
    X_bar[t,] <- colMeans(X_sub[indices,])
  }
  Cov_X_bar <- (1/T_dim) * t(X_bar) %*% X_bar
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric = TRUE)
  eigen_values <- Cov_X_bar_eigen$values / Cov_X_bar_eigen$values[1]
  K_hat <- sum(alpha_threshold < eigen_values)
  
  W_hat_tmp <- X_bar %*% Cov_X_bar_eigen$vectors[,1:K_hat]
  W_hat <- cbind((rep(1,(T_dim-1) )), W_hat_tmp[-T_dim,], W_hat_tmp[-1,])
  Pi_hat <- diag((T_dim-1)) - W_hat %*% solve(t(W_hat) %*% W_hat)  %*% t(W_hat)
  
  Y_hat <- rep(NA, (T_dim-1) * n_sub)
  X_hat <- matrix(NA, nrow = n_sub * (T_dim-1), ncol = p)
  for(i in 1:n_sub){
    index1 <- ((i-1) * (T_dim-1) + 1):(i * (T_dim-1))
    index2 <- ((i-1) * T_dim + 1):(i * T_dim -1)
    Y_hat[index1] <- Pi_hat %*% t(t(Y_sub[index2]))
    X_hat[index1,] <- Pi_hat %*% t(t(X_sub[index2,]))
  }
  list(n_sub = n_sub, X_sub = X_sub, Y_sub = Y_sub, X_sub_var = X_sub_var,
       X_bar = X_bar, K_hat = K_hat, X_hat = X_hat, Y_hat = Y_hat)
}

## ------------------------------------------------------------------------
## Least squares, all variants
## ------------------------------------------------------------------------
run_LS_variants <- function(data_sub, variants, alpha = c(0.01, 0.05, 0.1)){
  pr <- prep_hdcce(data_sub)
  X_hat <- pr$X_hat
  Y_hat <- pr$Y_hat
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  res_LS <- coef(fit_LS)
  resid_LS <- stats::residuals(fit_LS)
  XtX_inv <- summary(fit_LS)$cov.unscaled
  score <- X_hat * resid_LS
  score_covs <- score_cov_variants(score, pr$n_sub, T_dim - 1, DK_LAG, variants)
  sd_j <- sqrt(pr$X_sub_var)
  coef_scaled <- unname(res_LS) * sd_j
  
  out <- lapply(seq_len(nrow(variants)), function(r){
    Vcov <- XtX_inv %*% score_covs[[r]] %*% XtX_inv
    se_scaled <- sqrt(diag(Vcov)) * sd_j
    conf_min <- sweep(outer(se_scaled, qnorm(alpha/2)),     1, coef_scaled, "+")
    conf_max <- sweep(outer(se_scaled, qnorm(1 - alpha/2)), 1, coef_scaled, "+")
    data.frame(
      Variable      = colnames(pr$X_sub),
      se_type       = variants$se_type[r],
      kernel        = variants$kernel[r],
      h_T           = variants$h_T[r],
      rule          = variants$rule[r],
      n             = pr$n_sub,
      K_hat         = pr$K_hat,
      LS_coef       = round(coef_scaled, 4),
      LS_SE         = round(se_scaled, 4),
      LS_p_value    = round(2 * pnorm(-abs(coef_scaled / se_scaled)), 4),
      LS_Significance = make_stars(conf_min, conf_max),
      row.names     = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ------------------------------------------------------------------------
## Desparsified lasso, all variants (cv.glmnet fits shared; variance/kappa vary)
## ------------------------------------------------------------------------
run_lasso_variants <- function(data_sub, variants, alpha = c(0.01, 0.05, 0.1)){
  pr <- prep_hdcce(data_sub)
  n_sub <- pr$n_sub; X_sub <- pr$X_sub; Y_sub <- pr$Y_sub
  X_sub_var <- pr$X_sub_var; X_bar <- pr$X_bar; K_hat <- pr$K_hat
  X_hat <- pr$X_hat; Y_hat <- pr$Y_hat
  foldid <- rep(1:n_sub, each = (T_dim-1))
  
  res <- coef(
    glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                      standardize = FALSE, intercept = FALSE),
    s = "lambda.min"
  )[-1]
  res_scaled <- res * sqrt(X_sub_var)
  
  n_var   <- nrow(variants)
  n_alpha <- length(alpha)
  despar_beta_all <- matrix(NA_real_, nrow = p, ncol = n_var)
  Avar_all        <- matrix(NA_real_, nrow = p, ncol = n_var)
  conf_min_all    <- array(NA_real_, dim = c(p, n_alpha, n_var))
  conf_max_all    <- array(NA_real_, dim = c(p, n_alpha, n_var))
  
  for(COEF_INDEX in 1:p){
    
    Cov_X_bar_tilde <- (1/T_dim) * t(X_bar[,-COEF_INDEX]) %*% X_bar[,-COEF_INDEX]
    Cov_X_bar_tilde_eigen <- eigen(Cov_X_bar_tilde, symmetric = TRUE)
    W_tilde_tmp <- X_bar[,-COEF_INDEX] %*% Cov_X_bar_tilde_eigen$vectors[,1:K_hat]
    W_tilde <- cbind((rep(1,(T_dim-1) )), W_tilde_tmp[-1,], W_tilde_tmp[-T_dim,])
    Pi_tilde <- diag((T_dim-1)) -  W_tilde %*% solve(t(W_tilde) %*% W_tilde)  %*% t(W_tilde)
    
    Y_tilde <- rep(NA, (T_dim-1) * n_sub)
    X_tilde <- matrix(NA, nrow = n_sub * (T_dim-1), ncol = p)
    for(i in 1:n_sub){
      index1 <- ((i-1) * (T_dim-1) + 1):(i * (T_dim-1))
      index2 <- ((i-1) * T_dim + 1):(i * T_dim -1)
      Y_tilde[index1] <- Pi_tilde %*% t(t(Y_sub[index2]))
      X_tilde[index1,] <- Pi_tilde %*% t(t(X_sub[index2,]))
    }
    
    fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde,
                                   foldid = foldid, family = "gaussian",
                                   alpha = 1, intercept = FALSE, standardize = FALSE)
    coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
    yhat_Lasso <- stats::predict(fit_Lasso, newx = X_tilde,
                                 type = "response", s = "lambda.min")
    resid_Lasso <- Y_tilde - yhat_Lasso
    eps_mat <- matrix(resid_Lasso, nrow = (T_dim-1), ncol = n_sub)
    
    fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX],
                                        foldid = foldid, family = "gaussian",
                                        intercept = FALSE, standardize = FALSE)
    kappa_grid <- fit_node_Lasso$lambda
    kappa_cv_idx <- fit_node_Lasso$index[1]
    kappa_grid_len <- length(kappa_grid)
    
    var_scaled <- matrix(0, nrow = kappa_grid_len, ncol = n_var)
    for(k in 1:kappa_grid_len){
      yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                        type = "response", s = kappa_grid[k])
      resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso
      Delta_mat <- matrix(resid_node_Lasso, nrow = (T_dim-1), ncol = n_sub)
      G <- Delta_mat * eps_mat
      denom <- as.numeric(t(X_tilde[,COEF_INDEX]) %*% resid_node_Lasso)^2
      var_scaled[k, ] <- score_var_variants(G, DK_LAG, variants) / denom
    }
    
    for(r in seq_len(n_var)){
      v <- var_scaled[, r]
      V_TRUNC <- KAPPA_MULT * v[kappa_cv_idx]
      kappa_idx <- 1
      for(l in 1:kappa_grid_len){
        if(v[l] <= V_TRUNC){ kappa_idx <- l }
        if(v[l] > V_TRUNC){ break }
      }
      yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                        type = "response", s = kappa_grid[kappa_idx])
      resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso
      despar_beta <- as.numeric(
        coefs_Lasso[COEF_INDEX] +
          (t(resid_node_Lasso) %*% resid_Lasso) / (t(resid_node_Lasso) %*% X_tilde[, COEF_INDEX])
      )
      Avar <- sqrt(v[kappa_idx])
      conf_min_all[COEF_INDEX, , r] <- sqrt(X_sub_var[COEF_INDEX]) * (despar_beta + Avar * qnorm(alpha/2))
      conf_max_all[COEF_INDEX, , r] <- sqrt(X_sub_var[COEF_INDEX]) * (despar_beta + Avar * qnorm(1 - alpha/2))
      despar_beta_all[COEF_INDEX, r] <- despar_beta
      Avar_all[COEF_INDEX, r] <- Avar
    }
  }
  
  out <- lapply(seq_len(n_var), function(r){
    coef_d <- despar_beta_all[, r] * sqrt(X_sub_var)
    se_d   <- Avar_all[, r] * sqrt(X_sub_var)
    data.frame(
      Variable      = colnames(X_sub),
      se_type       = variants$se_type[r],
      kernel        = variants$kernel[r],
      h_T           = variants$h_T[r],
      rule          = variants$rule[r],
      n             = n_sub,
      K_hat         = K_hat,
      Lasso_coef_raw      = round(res_scaled, 4),
      Lasso_debiased_coef = round(coef_d, 4),
      Lasso_SE            = round(se_d, 4),
      Lasso_p_value       = round(2 * pnorm(-abs(coef_d / se_d)), 4),
      Lasso_Significance  = make_stars(conf_min_all[, , r], conf_max_all[, , r]),
      row.names     = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ------------------------------------------------------------------------
## Run for the three groups
## ------------------------------------------------------------------------
cat("\nRunning LS variants...\n")
ls_long <- do.call(rbind, lapply(GROUP_ORDER, function(g){
  cat("  group:", g, "\n")
  tab <- run_LS_variants(groups_list[[g]], variants)
  tab$Group <- g
  tab
}))

cat("\nRunning desparsified-lasso variants \n")
lasso_long <- do.call(rbind, lapply(GROUP_ORDER, function(g){
  cat("  group:", g, "\n")
  tab <- run_lasso_variants(groups_list[[g]], variants)
  tab$Group <- g
  tab
}))

order_cols <- function(tab){
  lead <- c("Group", "se_type", "kernel", "h_T", "rule", "Variable")
  tab <- tab[order(match(tab$Group, GROUP_ORDER), tab$se_type, tab$kernel, tab$h_T, tab$Variable), ]
  tab <- tab[, c(lead, setdiff(names(tab), lead))]
  rownames(tab) <- NULL
  tab
}
ls_long    <- order_cols(ls_long)
lasso_long <- order_cols(lasso_long)

## ------------------------------------------------------------------------
## Summaries
## ------------------------------------------------------------------------
## Compared to the DK_cluster/uniform/h_T=T baseline. SE_ratio<1 = smaller SE
## than baseline. LS coefficients don't vary; the lasso's debiased coef can.
summarise_bandwidth <- function(tab, coef_col, se_col, p_col){
  key  <- c("Group", "Variable", "se_type")
  base <- tab[tab$kernel == "uniform" & tab$h_T == T_eff, c(key, coef_col, se_col, p_col)]
  names(base)[4:6] <- c("coef_base", "se_base", "p_base")
  m <- merge(tab, base, by = key, sort = FALSE)
  m$SE_ratio    <- m[[se_col]] / m$se_base
  m$coef_change <- m[[coef_col]] - m$coef_base
  m$sig5        <- m[[p_col]] < 0.05
  m$sig5_base   <- m$p_base < 0.05
  
  combos <- unique(m[, c("Group", "se_type", "kernel", "h_T", "rule")])
  rows <- lapply(seq_len(nrow(combos)), function(r){
    s <- m[m$Group == combos$Group[r] & m$se_type == combos$se_type[r] &
             m$kernel == combos$kernel[r] & m$h_T == combos$h_T[r], ]
    data.frame(
      Group  = combos$Group[r], se_type = combos$se_type[r], kernel = combos$kernel[r],
      h_T    = combos$h_T[r],   rule    = combos$rule[r],
      n_sig_10pct = sum(s[[p_col]] < 0.10, na.rm = TRUE),
      n_sig_5pct  = sum(s[[p_col]] < 0.05, na.rm = TRUE),
      n_sig_1pct  = sum(s[[p_col]] < 0.01, na.rm = TRUE),
      n_sig5_changed_vs_baseline = sum(s$sig5 != s$sig5_base, na.rm = TRUE),
      median_SE_ratio = round(median(s$SE_ratio, na.rm = TRUE), 3),
      min_SE_ratio    = round(min(s$SE_ratio, na.rm = TRUE), 3),
      max_SE_ratio    = round(max(s$SE_ratio, na.rm = TRUE), 3),
      median_abs_coef_change = round(median(abs(s$coef_change), na.rm = TRUE), 5),
      max_abs_coef_change    = round(max(abs(s$coef_change), na.rm = TRUE), 5),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(match(out$Group, GROUP_ORDER), out$se_type, out$kernel, out$h_T), ]
  rownames(out) <- NULL
  out
}
ls_summary    <- summarise_bandwidth(ls_long,    "LS_coef",             "LS_SE",     "LS_p_value")
lasso_summary <- summarise_bandwidth(lasso_long, "Lasso_debiased_coef", "Lasso_SE",  "Lasso_p_value")
print(ls_summary)
print(lasso_summary)

## ------------------------------------------------------------------------
## Model fit (RMSE/AIC/BIC) for the two HD-CCE estimators
## ------------------------------------------------------------------------
## Same construction as compare_LS_lasso_fit() in the main replication script,
## reusing prep_hdcce() here instead of repeating Steps 1-2 inline.
compute_fit_stats_hdcce <- function(data_sub){
  pr <- prep_hdcce(data_sub)
  X_hat <- pr$X_hat
  Y_hat <- pr$Y_hat
  n_obs <- length(Y_hat)
  foldid <- rep(1:pr$n_sub, each = (T_dim-1))
  
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  RSS_LS <- sum(stats::residuals(fit_LS)^2)
  k_LS <- p
  
  fit_lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                                 standardize = FALSE, intercept = FALSE)
  beta_lasso <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  RSS_lasso <- sum((Y_hat - X_hat %*% beta_lasso)^2)
  k_lasso <- sum(beta_lasso != 0)
  
  fit_stats <- function(model_name, RSS, k){
    data.frame(
      Model = model_name,
      n     = n_obs,
      K_hat = pr$K_hat,
      k     = k,
      RMSE  = round(sqrt(RSS / n_obs), 4),
      AIC   = round(n_obs * log(RSS / n_obs) + 2 * k, 2),
      BIC   = round(n_obs * log(RSS / n_obs) + k * log(n_obs), 2)
    )
  }
  
  rbind(fit_stats("LS", RSS_LS, k_LS), fit_stats("Lasso", RSS_lasso, k_lasso))
}

cat("\nRunning model fit comparison...\n")
fit_stats_table <- do.call(rbind, lapply(GROUP_ORDER, function(g){
  tab <- compute_fit_stats_hdcce(groups_list[[g]])
  tab$Group <- g
  tab
}))
fit_stats_table <- fit_stats_table[, c("Group", setdiff(names(fit_stats_table), "Group"))]
rownames(fit_stats_table) <- NULL
print(fit_stats_table)

## ------------------------------------------------------------------------
## Residual diagnostics for the two HD-CCE estimators
## ------------------------------------------------------------------------
## Residuals of the main LS/lasso fits on the Pi_hat-projected data; residual-based,
## so independent of h_T and the SE variant. BP = pooled heteroskedasticity;
## Bartlett = variance equal across countries; Ljung-Box = per-country serial
## correlation; Pesaran CD = cross-sectional dependence.
run_diagnostics_hdcce <- function(data_sub, lb_lag = 4){
  pr <- prep_hdcce(data_sub)
  n_sub <- pr$n_sub
  X_hat <- pr$X_hat
  Y_hat <- pr$Y_hat
  n_obs <- length(Y_hat)
  foldid <- rep(1:n_sub, each = (T_dim-1))
  country_factor <- factor(rep(1:n_sub, each = (T_dim-1)))
  
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  resid_LS <- as.numeric(stats::residuals(fit_LS))
  
  fit_lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                                 standardize = FALSE, intercept = FALSE)
  beta_lasso <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  resid_lasso <- as.numeric(Y_hat - X_hat %*% beta_lasso)
  
  diagnose_one <- function(model_name, resid){
    ## Breusch-Pagan (pooled)
    bp_fit  <- lm(resid^2 ~ X_hat)
    bp_stat <- n_obs * summary(bp_fit)$r.squared
    bp_p    <- pchisq(bp_stat, df = p, lower.tail = FALSE)
    
    ## Bartlett's test across countries
    bart <- stats::bartlett.test(resid, country_factor)
    
    ## Per-country Ljung-Box
    lb_pvals <- numeric(n_sub)
    for(i in 1:n_sub){
      block <- ((i-1)*(T_dim-1)+1):(i*(T_dim-1))
      lb_pvals[i] <- stats::Box.test(resid[block], lag = lb_lag,
                                     type = "Ljung-Box")$p.value
    }
    
    ## Pesaran CD (column-major reshape puts country i's series in column i)
    resid_mat <- matrix(resid, nrow = T_dim - 1, ncol = n_sub)
    R_resid   <- cor(resid_mat)
    upper_sum <- sum(R_resid[upper.tri(R_resid)])
    CD_stat   <- sqrt(2 * (T_dim - 1) / (n_sub * (n_sub - 1))) * upper_sum
    CD_p      <- 2 * pnorm(-abs(CD_stat))
    
    data.frame(
      Model                = model_name,
      BP_stat              = round(bp_stat, 3),
      BP_p                 = round(bp_p, 4),
      Bartlett_stat        = round(unname(bart$statistic), 3),
      Bartlett_p           = round(bart$p.value, 4),
      LjungBox_pct_reject5 = round(100 * mean(lb_pvals < 0.05), 1),
      LjungBox_mean_p      = round(mean(lb_pvals), 4),
      CD_stat              = round(CD_stat, 3),
      CD_p                 = round(CD_p, 4),
      stringsAsFactors     = FALSE
    )
  }
  
  rbind(diagnose_one("LS", resid_LS), diagnose_one("Lasso", resid_lasso))
}

cat("\nRunning residual diagnostics...\n")
diagnostics_table <- do.call(rbind, lapply(GROUP_ORDER, function(g){
  tab <- run_diagnostics_hdcce(groups_list[[g]])
  tab$Group <- g
  tab
}))
diagnostics_table <- diagnostics_table[, c("Group", setdiff(names(diagnostics_table), "Group"))]
rownames(diagnostics_table) <- NULL
print(diagnostics_table)

## ------------------------------------------------------------------------
## Excel workbook
## ------------------------------------------------------------------------
XLSX_HEADER_STYLE <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                                           halign = "center", border = "Bottom")
compute_col_widths <- function(df, min_w = 8, max_w = 45){
  w <- vapply(seq_len(ncol(df)), function(j){
    header_w <- nchar(names(df)[j]) + 4
    cell_w   <- max(nchar(as.character(df[[j]])), 0, na.rm = TRUE) + 2
    max(header_w, cell_w)
  }, numeric(1))
  pmin(pmax(w, min_w), max_w)
}
add_results_sheet <- function(wb, sheet, df){
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df, withFilter = TRUE, headerStyle = XLSX_HEADER_STYLE)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = compute_col_widths(df))
  invisible(wb)
}

readme <- data.frame(
  Item = c("Purpose", "Bandwidth h_T", "Window shape (kernel)", "SE types", "Baseline",
           "Summary columns", "Estimates vs. SEs", "Not varied", "Caveats", "Fit_stats",
           "Residual_diagnostics"),
  Detail = c(
    "Sensitivity of the HD-CCE LS and desparsified-lasso standard errors to the HAC bandwidth h_T, for All / Carbon tax / No carbon tax.",
    paste0("h_T is the width of the within-country window (T = ", T_eff, " usable periods). Values: Linton h_T = T; Newey-West (1994) ",
           if(HT_ROUND_NW) "floor(4(T/100)^(2/9)) (rounded down, as in the formula)" else "4(T/100)^(2/9) (unrounded)",
           "; T^(1/3) rate ", if(HT_ROUND_T13) "floor(T^(1/3)) (rounded down)" else "T^(1/3) with constant 1 (unrounded)",
           ". See Bandwidth_rules. Only lags j <= h_T enter and Bartlett weights use the (possibly fractional) h_T, weight 1 - j/(h_T+1); a uniform window is unaffected by the fractional part."),
    "uniform: weight 1 inside the window (as in the original band_mask code; h_T = T reproduces the earlier clustered HAC). bartlett: weight 1 - lag/(h_T+1), the kernel the Newey-West rule is designed for; always positive semidefinite.",
    paste0("cluster: within-country HAC with window h_T, independence across countries (Linton's N^HAC). DK_cluster: same within-country window plus Driscoll-Kraay Bartlett weights (lag ", DK_LAG, ") for cross-country pairs, as in the thesis tables."),
    "DK_cluster, uniform window, h_T = T: reproduces the existing thesis results and is the reference for all ratios. Check that these rows match the Estimator_comparison sheet of the main workbook.",
    "n_sig_*: number of the 22 drivers significant at 10/5/1%. n_sig5_changed_vs_baseline: drivers whose 5% significance differs from the baseline. SE_ratio: SE divided by the baseline SE (same SE type). abs_coef_change: absolute change in the coefficient relative to the baseline.",
    "LS coefficients do not depend on h_T. Debiased-lasso coefficients can change slightly because the node-wise kappa selection uses the variance.",
    "The tuning constants (KAPPA_MULT = 1.25 in the node-wise rule, cross-validated lambda and kappa) and the number of factors are held fixed; their effect on both estimates and SEs is not examined here.",
    "Smaller h_T lowers the SE only by ignoring longer within-country serial dependence; a smaller SE at a smaller bandwidth is not evidence of greater precision. Truncated (uniform) windows are not guaranteed positive semidefinite; negative variances are floored at 0 (lasso) or eigenvalues clipped (LS). Baseline choice (h_T = T) should not be changed in light of the results.",
    "RMSE, AIC and BIC for the LS and lasso fits on the Pi_hat-projected data (k = 22 for LS, nonzero coefficients at lambda.min for lasso). Residual-based, so independent of h_T and of the SE variant.",
    "Breusch-Pagan (pooled heteroskedasticity), Bartlett (residual variance equal across countries), Ljung-Box (lag 4, per country: % of countries rejecting no autocorrelation at 5% and mean p-value) and Pesaran CD (cross-sectional dependence) for the LS and lasso residuals of the projected regression. Residual-based, so independent of h_T and of the SE variant."
  ),
  stringsAsFactors = FALSE
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "README")
openxlsx::writeData(wb, "README", readme, headerStyle = XLSX_HEADER_STYLE)
openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(26, 110))
openxlsx::addStyle(wb, "README", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                   rows = 2:(nrow(readme) + 1), cols = 1:2, gridExpand = TRUE)
add_results_sheet(wb, "Bandwidth_rules", ht_grid)
add_results_sheet(wb, "Lasso_summary",   lasso_summary)
add_results_sheet(wb, "LS_summary",      ls_summary)
add_results_sheet(wb, "Lasso_long",      lasso_long)
add_results_sheet(wb, "LS_long",         ls_long)
add_results_sheet(wb, "Fit_stats",       fit_stats_table)
add_results_sheet(wb, "Residual_diagnostics", diagnostics_table)
openxlsx::saveWorkbook(wb, EXCEL_FILE, overwrite = TRUE)
cat("\nResults saved to", EXCEL_FILE, "\n")