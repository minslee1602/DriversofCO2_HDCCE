# Replication File for: Drivers of Carbon Emissions in OECD and BRICS Countries - A High-Dimensional CCE Panel Approach
# Author: Minsoo Lee
# Supervised by: Dr. Jessica Leung & Prof. Param Silvapulle
# This file obtains HD-CCE least-squares and lasso estimates for a panel dataset of 42 countries divided into carbon-tax
# and non-carbon tax groups. RMSE, AIC/BIC are obtained for model comparison, and residual diagnostic checks using Breusch Pagan,
# Bartlett, Ljung-box and Pesaran (2021)'s CD test are included.


## Data
data <- openxlsx::read.xlsx("Final_Cabon_Emm_disaggregate_zscored.xlsx", colNames = TRUE, detectDates = TRUE)
p <- 22 # drivers
T_dim <- 30 #1990-2019

# Driver columns = every "*_z" column except the outcome's own z-score
# (carbon_emissions_pc_z), which we don't use. Selecting by name rather than
# position because the z-scored workbook interleaves each raw column with its
# _z column, so the drivers are no longer in one contiguous block.
z_cols       <- grep("_z$", names(data), value = TRUE)
driver_cols_z <- setdiff(z_cols, "carbon_emissions_pc_z")
stopifnot(length(driver_cols_z) == p)

year_vec  <- data$Year
uniq_yrs  <- sort(unique(year_vec))
T_window  <- uniq_yrs[(length(uniq_yrs) - (T_dim - 1)):length(uniq_yrs)]
data_9019 <- data[data$Year %in% T_window, ]
# Country.Code
data_9019 <- data_9019[order(data_9019$Country.Code), ]
# Exclude Costa Rica
data_9019 <- data_9019[(data_9019$Country.Code != "CRI"), ]

## Carbon.Tax subgroups
# Carbon.Tax is time-invariant per country in this file, so filtering the
# already-windowed/sorted/CRI-excluded panel by it keeps balanced n*T_dim
# blocks for each subgroup
data_all    <- data_9019
data_tax    <- data_9019[data_9019$Carbon.Tax == 1, ]
data_no_tax <- data_9019[data_9019$Carbon.Tax == 0, ]

cat("All countries:  n =", length(unique(data_all$Country.Code)), "\n")
cat("Carbon tax = 1: n =", length(unique(data_tax$Country.Code)), "\n")
cat("Carbon tax = 0: n =", length(unique(data_no_tax$Country.Code)), "\n")

## Adjusting alpha for tau threshold
K_ALPHA_THRESHOLD <- 0.01

## Standard-error settings
# DK_LAG: Bartlett lag window m for the time dimension. Driscoll & Kraay
# require m(T) = O(T^(1/4)) (about 2 for T_dim - 1 = 29) and use m = 2 in
# their empirical examples. See the lag-sensitivity block near the end.
DK_LAG   <- 2
SE_TYPES <- "DK_cluster"
stopifnot(DK_LAG >= 0, DK_LAG < T_dim - 1)

GROUP_ORDER <- c("All", "Carbon tax", "No carbon tax")
EXCEL_FILE  <- "./tables/LassoandLS_HDCCE_results.xlsx"


lrv_bartlett <- function(H, m){
  Tn <- nrow(H)
  stopifnot(m >= 0, m < Tn)
  S <- crossprod(H)
  if(m > 0) for(j in 1:m){
    w  <- 1 - j/(m + 1)
    Gj <- crossprod(H[(j+1):Tn, , drop = FALSE],
                    H[1:(Tn-j), , drop = FALSE])
    S  <- S + w * (Gj + t(Gj))
  }
  S
}

# Scalar version for a Tn x units matrix G: sum over columns of each column's
# Bartlett long-run variance (used by the desparsified-lasso variance)
lrv_scalar <- function(G, m){
  Tn <- nrow(G)
  stopifnot(m >= 0, m < Tn)
  out <- sum(G^2)
  if(m > 0) for(j in 1:m){
    out <- out + 2 * (1 - j/(m + 1)) *
      sum(G[(j+1):Tn, , drop = FALSE] * G[1:(Tn-j), , drop = FALSE])
  }
  out
}

# Sum a country-major stacked (n_sub*Tn) x R score matrix across countries at
# each date -> Tn x R
sum_over_countries <- function(score, n_sub, Tn){
  R <- ncol(score)
  matrix(apply(array(score, dim = c(Tn, n_sub, R)), c(1, 3), sum),
         nrow = Tn, ncol = R)
}

# Force symmetry and clip negative eigenvalues (the two-way combination
# S_cluster + S_DK - S_within is not guaranteed positive semidefinite)
psd_truncate <- function(M){
  M  <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
}

# R x R score-covariance matrices for all three SE types (LS estimator).
# Returns a named list: cluster, DK, DK_cluster.
score_cov_matrices <- function(score, n_sub, Tn, m){
  blocks <- lapply(seq_len(n_sub),
                   function(i) score[((i-1)*Tn + 1):(i*Tn), , drop = FALSE])
  S_cl     <- Reduce(`+`, lapply(blocks, function(B) tcrossprod(colSums(B))))
  S_dk     <- lrv_bartlett(sum_over_countries(score, n_sub, Tn), m)
  S_within <- Reduce(`+`, lapply(blocks, lrv_bartlett, m = m))
  list(cluster    = S_cl,
       DK         = S_dk,
       DK_cluster = psd_truncate(S_cl + S_dk - S_within))
}

# Scalar score variances for all three SE types (desparsified lasso).
# G is the (T_dim-1) x n_sub matrix of node-residual * lasso-residual.
# Returns a named vector: cluster, DK, DK_cluster.
score_var_all <- function(G, m){
  S_cl     <- sum(colSums(G)^2)
  S_dk     <- lrv_scalar(matrix(rowSums(G), ncol = 1), m)
  S_within <- lrv_scalar(G, m)
  c(cluster    = S_cl,
    DK         = S_dk,
    DK_cluster = max(S_cl + S_dk - S_within, 0))
}

# Significance stars from a p x length(alpha) pair of CI matrices, with alpha
# ordered c(0.01, 0.05, 0.1): "***" p<0.01, "**" p<0.05, "*" p<0.10
make_stars <- function(conf_min, conf_max){
  n_alpha <- ncol(conf_min)
  codes <- numeric(nrow(conf_min))
  for(j in seq_len(nrow(conf_min))){
    for(a in seq_len(n_alpha)){
      if(isTRUE(conf_min[j,a] < 0 && 0 < conf_max[j,a])) codes[j] <- a
    }
  }
  vapply(codes, function(cd) paste(rep("*", n_alpha - cd), collapse = ""),
         character(1))
}

## Scree plot of the K_hat eigenvalues 
make_scree_plot <- function(data_sub, group_name,
                            alpha_threshold = K_ALPHA_THRESHOLD,
                            file = NULL){
  
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  
  X_bar <- matrix(NA, ncol = p, nrow = T_dim)
  for(t in 1:T_dim){
    indices <- seq(t, n_sub * T_dim, by = T_dim)
    X_bar[t,] <- colMeans(X_sub[indices,])
  }
  Cov_X_bar <- (1/T_dim) * t(X_bar) %*% X_bar
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric = TRUE)
  psi_hat <- Cov_X_bar_eigen$values            # psi_hat_1 >= ... >= psi_hat_p
  tau <- alpha_threshold * psi_hat[1]          # tau = alpha * psi_hat_1
  K_hat <- sum(psi_hat > tau)
  
  if(!is.null(file)){
    png(file, width = 1600, height = 1200, res = 200)
    on.exit(dev.off(), add = TRUE)
  }
  
  plot(seq_len(p), psi_hat, type = "b", pch = 16,
       xlab = "j", ylab = expression(hat(psi)[j]),
       main = paste0("Scree plot -- ", group_name, " (K_hat = ", K_hat, ")"))
  abline(h = tau, lty = 2, col = "red")
  if(K_hat > 0){
    points(seq_len(K_hat), psi_hat[seq_len(K_hat)], pch = 16, col = "blue")
  }
  legend("topright", bty = "n",
         legend = c(paste0("tau = ", signif(tau, 4), "  (alpha = ", alpha_threshold, ")"),
                    paste0("K_hat = ", K_hat, " retained")),
         col = c("red", "blue"), lty = c(2, NA), pch = c(NA, 16))
  
  invisible(list(psi_hat = psi_hat, tau = tau, K_hat = K_hat))
}

dir.create("./figures", showWarnings = FALSE)
make_scree_plot(data_all,    "All countries",  file = "./figures/ScreePlot_All.png")
make_scree_plot(data_tax,    "Carbon tax = 1", file = "./figures/ScreePlot_CarbonTax.png")
make_scree_plot(data_no_tax, "Carbon tax = 0", file = "./figures/ScreePlot_NoCarbonTax.png")

## Least squares HD-CCE estimator 
## Returns one row per (Variable, SE_type). Point estimates do not depend on
## SE_type; only LS_SE, the CIs and the stars do.
run_LS_hdcce <- function(data_sub, alpha = c(0.01, 0.05, 0.1),
                         alpha_threshold = K_ALPHA_THRESHOLD,
                         drop_vars = NULL, 
                         se_types = SE_TYPES,
                         dk_lag = DK_LAG){
  
  # drop_vars: character vector of driver_cols_z names to exclude from this
  # fit (e.g. drop_vars = "gas_production_pc_z"), for robustness-checking
  # whether a variable's coefficient/sign is stable to a correlated driver
  # being removed.
  driver_cols_z_use <- setdiff(driver_cols_z, drop_vars)
  p_local <- length(driver_cols_z_use)
  
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z_use]
  X_sub <- X_sub[, order(names(X_sub))]
  Y_sub <- data_sub$carbon_emissions_pc
  X_sub_var <- apply(X_sub, MARGIN = 2, FUN = var)
  
  # csa of regressors
  X_bar <- matrix(NA, ncol = p_local, nrow = T_dim)
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
  X_hat <- matrix(NA, nrow = n_sub * (T_dim-1), ncol = p_local)
  for(i in 1:n_sub){
    index1 <- ((i-1) * (T_dim-1) + 1):(i * (T_dim-1))
    index2 <- ((i-1) * T_dim + 1):(i * T_dim -1)
    Y_hat[index1] <- Pi_hat %*% t(t(Y_sub[index2]))
    X_hat[index1,] <- Pi_hat %*% t(t(X_sub[index2,]))
  }
  
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  res_LS <- coef(fit_LS)
  names(res_LS) <- colnames(X_sub)
  resid_LS <- stats::residuals(fit_LS)
  
  # Sandwich variance: (X'X)^-1 S (X'X)^-1, where S is the score covariance
  # for each SE type (cluster / DK / DK_cluster, see helpers above). The
  # "cluster" S reproduces the earlier country-clustered Liang-Zeger loop.
  XtX_inv <- summary(fit_LS)$cov.unscaled
  score <- X_hat * resid_LS                       # n_obs x p_local
  score_covs <- score_cov_matrices(score, n_sub, T_dim - 1, dk_lag)
  
  res_LS_scaled <- res_LS * sqrt(X_sub_var)
  
  out <- lapply(se_types, function(s){
    Vcov_LS <- XtX_inv %*% score_covs[[s]] %*% XtX_inv
    SE_LS <- sqrt(diag(Vcov_LS))
    names(SE_LS) <- colnames(X_sub)
    SE_LS_scaled <- SE_LS * sqrt(X_sub_var)
    
    conf_LS_min <- sweep(outer(SE_LS_scaled, qnorm(alpha/2)),     1, res_LS_scaled, "+")
    conf_LS_max <- sweep(outer(SE_LS_scaled, qnorm(1 - alpha/2)), 1, res_LS_scaled, "+")
    
    data.frame(
      Variable        = colnames(X_sub),
      SE_type         = s,
      n               = n_sub,
      K_hat           = K_hat,
      LS_coef         = round(res_LS_scaled, 4),
      LS_SE           = round(SE_LS_scaled, 4),
      LS_CI95_lower   = round(conf_LS_min[,2], 4),
      LS_CI95_upper   = round(conf_LS_max[,2], 4),
      LS_Significance = make_stars(conf_LS_min, conf_LS_max),
      row.names       = NULL
    )
  })
  do.call(rbind, out)
}

## Desparsified lasso HD-CCE estimator. The variance normalization N^HAC is
## computed under all three SE types inside the same node-wise loop, so the
## expensive cv.glmnet fits are shared; only the variance calculation and the
## resulting kappa selection are repeated per SE type.
##   cluster    : sum over countries of (sum over t of g_it)^2, i.e. h_T = T_dim - 1
##                with no truncation (same as the earlier version)
##   DK         : Bartlett long-run variance (lag window dk_lag) of the
##                cross-country sums g_t = sum_i g_it
##   DK_cluster : combination described in the header
## where g_it = node-residual_it * lasso-residual_it. Steps 1-2 (K_hat, Pi_hat)
## are re-estimated separately within each subgroup.
##
## Regressors are already standardized (the *_z columns), so the manual
## variance-1 normalize_cols() step used in the raw-data version has been
## removed: X_hat / X_tilde are fed into glmnet directly, and coefficients
## are no longer divided back by a normalize_cols() scale.
run_desparsified_lasso_hdcce <- function(data_sub, alpha = c(0.01, 0.05, 0.1),
                                         alpha_threshold = K_ALPHA_THRESHOLD,
                                         se_types = SE_TYPES,
                                         dk_lag = DK_LAG){
  
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  Y_sub <- data_sub$carbon_emissions_pc
  X_sub_var <- apply(X_sub, MARGIN = 2, FUN = var)
  
  # csa of regressors
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
  
  foldid <- rep(1:n_sub, each = (T_dim-1))
  
  ## Main lasso, on the already-standardized (z-scored) regressors directly
  res <- coef(
    glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                      standardize = FALSE, intercept = FALSE),
    s = "lambda.min"
  )[-1]
  res_scaled <- res * sqrt(X_sub_var)
  names(res_scaled) <- colnames(X_sub)
  
  n_se <- length(se_types)
  n_alpha <- length(alpha)
  despar_beta_all <- matrix(NA_real_, nrow = p, ncol = n_se, dimnames = list(NULL, se_types))
  Avar_all        <- matrix(NA_real_, nrow = p, ncol = n_se, dimnames = list(NULL, se_types))
  conf_band_min_all <- setNames(replicate(n_se, matrix(NA_real_, nrow = p, ncol = n_alpha),
                                          simplify = FALSE), se_types)
  conf_band_max_all <- setNames(replicate(n_se, matrix(NA_real_, nrow = p, ncol = n_alpha),
                                          simplify = FALSE), se_types)
  
  for(COEF_INDEX in 1:p){
    
    # Empirical covariance matrix and eigenstructure, leave-one-out
    Cov_X_bar_tilde <- (1/T_dim) * t(X_bar[,-COEF_INDEX]) %*% X_bar[,-COEF_INDEX]
    Cov_X_bar_tilde_eigen <- eigen(Cov_X_bar_tilde, symmetric = TRUE)
    W_tilde_tmp <- X_bar[,-COEF_INDEX] %*% Cov_X_bar_tilde_eigen$vectors[,1:K_hat]
    W_tilde <- cbind((rep(1,(T_dim-1) )), W_tilde_tmp[-1,], W_tilde_tmp[-T_dim,])
    Pi_tilde <- diag((T_dim-1)) -  W_tilde %*% solve(t(W_tilde) %*% W_tilde)  %*% t(W_tilde)
    
    # Project the data
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
    
    # Variance along the kappa grid, one column per SE type
    var_scaled <- matrix(0, nrow = kappa_grid_len, ncol = n_se,
                         dimnames = list(NULL, se_types))
    for(k in 1:kappa_grid_len){
      yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                        type = "response", s = kappa_grid[k])
      resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso
      
      Delta_mat <- matrix(resid_node_Lasso, nrow = (T_dim-1), ncol = n_sub)
      G <- Delta_mat * eps_mat
      denom <- as.numeric(t(X_tilde[,COEF_INDEX]) %*% resid_node_Lasso)^2
      var_scaled[k, ] <- score_var_all(G, dk_lag)[se_types] / denom
    }
    
    # Kappa selection, debiased coefficient and CIs, separately per SE type
    for(s in se_types){
      v <- var_scaled[, s]
      V_TRUNC <- 1.25 * v[kappa_cv_idx]
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
      
      conf_band_min_all[[s]][COEF_INDEX, ] <- sqrt(X_sub_var[COEF_INDEX]) * (despar_beta + Avar * qnorm(alpha/2))
      conf_band_max_all[[s]][COEF_INDEX, ] <- sqrt(X_sub_var[COEF_INDEX]) * (despar_beta + Avar * qnorm(1 - alpha/2))
      despar_beta_all[COEF_INDEX, s] <- despar_beta
      Avar_all[COEF_INDEX, s] <- Avar
    }
  }
  
  # Same column names as the main script's Table 8.1, no prefix needed since
  # they don't collide with the LS_* columns from run_LS_hdcce()
  out <- lapply(se_types, function(s){
    data.frame(
      Variable      = colnames(X_sub),
      SE_type       = s,
      Lasso_coef    = round(res_scaled, 4),
      Debiased_coef = round(despar_beta_all[, s] * sqrt(X_sub_var), 4),
      Std_error     = round(Avar_all[, s] * sqrt(X_sub_var), 4),
      CI95_lower    = round(conf_band_min_all[[s]][, 2], 4),
      CI95_upper    = round(conf_band_max_all[[s]][, 2], 4),
      Significance  = make_stars(conf_band_min_all[[s]], conf_band_max_all[[s]]),
      row.names     = NULL
    )
  })
  do.call(rbind, out)
}

## Run both estimators for the three subgroups (each returns all SE types)
table_all    <- run_LS_hdcce(data_all)
table_tax    <- run_LS_hdcce(data_tax)
table_no_tax <- run_LS_hdcce(data_no_tax)

lasso_all    <- run_desparsified_lasso_hdcce(data_all)
lasso_tax    <- run_desparsified_lasso_hdcce(data_tax)
lasso_no_tax <- run_desparsified_lasso_hdcce(data_no_tax)

# Merge the desparsified lasso columns onto the LS table for each subgroup,
# matching on both Variable and SE_type
table_all    <- merge(table_all,    lasso_all,    by = c("Variable", "SE_type"), sort = FALSE)
table_tax    <- merge(table_tax,    lasso_tax,    by = c("Variable", "SE_type"), sort = FALSE)
table_no_tax <- merge(table_no_tax, lasso_no_tax, by = c("Variable", "SE_type"), sort = FALSE)

table_all$Group    <- "All"
table_tax$Group    <- "Carbon tax"
table_no_tax$Group <- "No carbon tax"

# *** p < 0.01, ** p < 0.05, * p < 0.10 (based on the 99%/95%/90% CIs above)
subgroup_table <- rbind(table_all, table_tax, table_no_tax)
subgroup_table <- subgroup_table[order(match(subgroup_table$Group, GROUP_ORDER),
                                       match(subgroup_table$SE_type, SE_TYPES),
                                       subgroup_table$Variable), ]
lead_cols <- c("Group", "SE_type", "Variable")
subgroup_table <- subgroup_table[, c(lead_cols, setdiff(names(subgroup_table), lead_cols))]
rownames(subgroup_table) <- NULL
print(subgroup_table)

dir.create("./tables", showWarnings = FALSE)

## ------------------------------------------------------------------------
## LS vs. lasso comparison tables
## ------------------------------------------------------------------------

build_estimator_comparison <- function(tbl, se_type){
  sub <- tbl[tbl$SE_type == se_type, ]
  sub <- sub[order(match(sub$Group, GROUP_ORDER), sub$Variable), ]
  data.frame(
    Group                = sub$Group,
    Variable             = sub$Variable,
    n                    = sub$n,
    K_hat                = sub$K_hat,
    LS_coef              = sub$LS_coef,
    LS_SE                = sub$LS_SE,
    LS_CI95_lower        = sub$LS_CI95_lower,
    LS_CI95_upper        = sub$LS_CI95_upper,
    LS_t                 = round(sub$LS_coef / sub$LS_SE, 2),
    LS_p_value           = round(2 * pnorm(-abs(sub$LS_coef / sub$LS_SE)), 4),
    LS_Significance      = sub$LS_Significance,
    Lasso_coef_raw       = sub$Lasso_coef,
    Lasso_debiased_coef  = sub$Debiased_coef,
    Lasso_SE             = sub$Std_error,
    Lasso_CI95_lower     = sub$CI95_lower,
    Lasso_CI95_upper     = sub$CI95_upper,
    Lasso_t              = round(sub$Debiased_coef / sub$Std_error, 2),
    Lasso_p_value        = round(2 * pnorm(-abs(sub$Debiased_coef / sub$Std_error)), 4),
    Lasso_Significance   = sub$Significance,
    SE_ratio_lasso_to_LS = round(sub$Std_error / sub$LS_SE, 3),
    Sign_agree           = ifelse(sign(sub$LS_coef) == sign(sub$Debiased_coef), "Yes", "No"),
    row.names            = NULL
  )
}
estimator_comparisons <- setNames(lapply(SE_TYPES, function(s) build_estimator_comparison(subgroup_table, s)),
                                  SE_TYPES)
estimator_comparison <- estimator_comparisons[["DK_cluster"]]
print(estimator_comparison)

## Excel workbook
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

# Write a data frame to a new sheet: styled header, autofilter, frozen top
# row, explicit column widths
add_results_sheet <- function(wb, sheet, df){
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df, withFilter = TRUE, headerStyle = XLSX_HEADER_STYLE)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = compute_col_widths(df))
  invisible(wb)
}

readme <- data.frame(
  Item = c("Purpose",
           "Estimators",
           "Standard errors",
           "DK lag window (m)",
           "Estimator_comparison sheet",
           "t-statistics and p-values",
           "Significance stars",
           "Coefficient scale",
           "Lasso columns",
           "Other sheets",
           "Caveats"),
  Detail = c(
    "HD-CCE least-squares (LS) and desparsified-lasso estimates of the drivers of per-capita carbon emissions, for All countries, Carbon tax and No carbon tax subgroups (T = 30 years, p = 22 z-scored drivers).",
    "LS = HD-CCE least squares on the Pi_hat-projected data. Lasso = desparsified (debiased) lasso HD-CCE with node-wise lasso; CV folds are countries.",
    "DK_cluster for both estimators: two-way combination of country clustering and Driscoll-Kraay (1998) SEs (Thompson 2011 / Cameron-Gelbach-Miller 2011 style, Bartlett lags in time). Same-country pairs weighted 1 at all lags, cross-country pairs get DK weights. Negative eigenvalues clipped (LS) / variance floored at 0 (lasso).",
    paste0("m = ", DK_LAG, " (Bartlett weights 1 - j/(m+1)). Sheet DK_lag_sensitivity_LS reruns LS with m = 1..4."),
    "LS and lasso side by side, all groups stacked. The same SE type and lag window are used for both estimators, so LS_SE and Lasso_SE are directly comparable. SE_ratio_lasso_to_LS > 1: debiased lasso less precise. Sign_agree: LS and debiased-lasso coefficients have the same sign.",
    "t = coefficient / SE; p-values are two-sided from the standard normal, matching the normal-quantile CIs used for the stars.",
    "*** p < 0.01, ** p < 0.05, * p < 0.10, from the 99% / 95% / 90% normal CIs. No stars = not significant at 10%.",
    "Coefficients and SEs are per one-standard-deviation change in the driver (multiplied by the driver's SD). Outcome is in projected carbon_emissions_pc units.",
    "Lasso_coef_raw is the un-debiased lasso coefficient (no valid SE). Inference uses Lasso_debiased_coef and Lasso_SE. The debiased coefficient can differ slightly across SE types because the node-wise kappa selection uses the variance.",
    "All_results (every column, long format), DK_lag_sensitivity_LS, Fit_stats (RMSE/AIC/BIC), VIF, Diagnostics (BP, Bartlett, Ljung-Box, Pesaran CD). Panel unit-root tests are run separately in Stata (see xtcips_crosscheck.do), not in this script.",
    "DK asymptotics are for T -> infinity with fixed N; here N is about 20 per subgroup and T-1 = 29, and the moment series come from Pi_hat-projected data (not covered by the original DK theory). Borderline stars deserve caution."
  ),
  stringsAsFactors = FALSE
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "README")
openxlsx::writeData(wb, "README", readme, headerStyle = XLSX_HEADER_STYLE)
openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(26, 110))
openxlsx::addStyle(wb, "README", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                   rows = 2:(nrow(readme) + 1), cols = 1:2, gridExpand = TRUE)

add_results_sheet(wb, "Estimator_comparison", estimator_comparison)
add_results_sheet(wb, "All_results",          subgroup_table)

# Save now so the main comparison is on disk even if a later step fails
openxlsx::saveWorkbook(wb, EXCEL_FILE, overwrite = TRUE)
cat("\nMain comparison sheets saved to", EXCEL_FILE, "\n")

## DK lag-window sensitivity for lasso
groups_list <- list("All" = data_all, "Carbon tax" = data_tax, "No carbon tax" = data_no_tax)

dk_lag_sens <- do.call(rbind, lapply(1:4, function(m){
  do.call(rbind, lapply(names(groups_list), function(g){
    tab <- run_LS_hdcce(groups_list[[g]], se_types = "DK_cluster", dk_lag = m)
    tab$Group  <- g
    tab$dk_lag <- m
    tab
  }))
}))
dk_lag_sens <- dk_lag_sens[order(match(dk_lag_sens$Group, GROUP_ORDER), dk_lag_sens$SE_type,
                                 dk_lag_sens$dk_lag, dk_lag_sens$Variable), ]
lead_cols_sens <- c("Group", "SE_type", "dk_lag", "Variable")
dk_lag_sens <- dk_lag_sens[, c(lead_cols_sens, setdiff(names(dk_lag_sens), lead_cols_sens))]
rownames(dk_lag_sens) <- NULL
print(dk_lag_sens)

## Model fit comparison
compare_LS_lasso_fit <- function(data_sub, alpha_threshold = K_ALPHA_THRESHOLD){
  
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  Y_sub <- data_sub$carbon_emissions_pc
  
  # Steps 1-2: K_hat, Pi_hat -- identical construction to run_LS_hdcce()/
  # run_desparsified_lasso_hdcce() above
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
  n_obs <- length(Y_hat)
  foldid <- rep(1:n_sub, each = (T_dim-1))
  
  ## LS fit -- same lm() call as run_LS_hdcce()
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  RSS_LS <- sum(stats::residuals(fit_LS)^2)
  k_LS <- p
  
  ## Main lasso fit -- same cv.glmnet() call as the Lasso_coef column in
  ## run_desparsified_lasso_hdcce()
  fit_lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                                 standardize = FALSE, intercept = FALSE)
  beta_lasso <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  fitted_lasso <- X_hat %*% beta_lasso
  RSS_lasso <- sum((Y_hat - fitted_lasso)^2)
  k_lasso <- sum(beta_lasso != 0)
  
  fit_stats <- function(model_name, RSS, k){
    data.frame(
      Model = model_name,
      n     = n_obs,
      K_hat = K_hat,
      k     = k,
      RMSE  = round(sqrt(RSS / n_obs), 4),
      AIC   = round(n_obs * log(RSS / n_obs) + 2 * k, 2),
      BIC   = round(n_obs * log(RSS / n_obs) + k * log(n_obs), 2)
    )
  }
  
  rbind(
    fit_stats("LS",    RSS_LS,    k_LS),
    fit_stats("Lasso", RSS_lasso, k_lasso)
  )
}

fitstats_all    <- compare_LS_lasso_fit(data_all)
fitstats_tax    <- compare_LS_lasso_fit(data_tax)
fitstats_no_tax <- compare_LS_lasso_fit(data_no_tax)

fitstats_all$Group    <- "All"
fitstats_tax$Group    <- "Carbon tax"
fitstats_no_tax$Group <- "No carbon tax"

fit_stats_table <- rbind(fitstats_all, fitstats_tax, fitstats_no_tax)
fit_stats_table <- fit_stats_table[, c("Group", setdiff(names(fit_stats_table), "Group"))]
print(fit_stats_table)

## VIF test
compute_VIF <- function(data_sub){
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  R <- cor(as.matrix(X_sub))
  vif_vals <- tryCatch(diag(solve(R)),
                       error = function(e) rep(NA_real_, ncol(X_sub)))
  data.frame(Variable = colnames(X_sub), VIF = round(vif_vals, 3))
}

vif_all    <- compute_VIF(data_all)
vif_tax    <- compute_VIF(data_tax)
vif_no_tax <- compute_VIF(data_no_tax)

vif_all$Group    <- "All"
vif_tax$Group    <- "Carbon tax"
vif_no_tax$Group <- "No carbon tax"

vif_table <- rbind(vif_all, vif_tax, vif_no_tax)
vif_table <- vif_table[, c("Group", "Variable", "VIF")]
print(vif_table)

## Residual diagnostics
run_diagnostics <- function(data_sub, alpha_threshold = K_ALPHA_THRESHOLD, lb_lag = 4){
  
  n_sub <- length(unique(data_sub$Country.Code))
  X_sub <- data_sub[, driver_cols_z]
  X_sub <- X_sub[, order(names(X_sub))]
  Y_sub <- data_sub$carbon_emissions_pc
  
  # Steps 1-2: K_hat, Pi_hat -- identical construction to run_LS_hdcce()/
  # run_desparsified_lasso_hdcce()/compare_LS_lasso_fit() above
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
  n_obs <- length(Y_hat)
  foldid <- rep(1:n_sub, each = (T_dim-1))
  country_factor <- factor(rep(1:n_sub, each = (T_dim-1)))
  
  ## LS residuals -- same lm() call as run_LS_hdcce()/compare_LS_lasso_fit()
  fit_LS <- lm(Y_hat ~ X_hat - 1)
  resid_LS <- stats::residuals(fit_LS)
  
  ## Main lasso residuals -- same cv.glmnet() call as compare_LS_lasso_fit()
  fit_lasso <- glmnet::cv.glmnet(x = X_hat, y = Y_hat, foldid = foldid,
                                 standardize = FALSE, intercept = FALSE)
  beta_lasso <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  resid_lasso <- as.numeric(Y_hat - X_hat %*% beta_lasso)
  
  diagnose_one <- function(model_name, resid){
    
    ## Breusch-Pagan (pooled): BP = n_obs * R^2_aux ~ chi-sq(p)
    bp_fit  <- lm(resid^2 ~ X_hat)
    bp_stat <- n_obs * summary(bp_fit)$r.squared
    bp_p    <- pchisq(bp_stat, df = p, lower.tail = FALSE)
    
    ## Bartlett's test: is residual variance equal across countries?
    bart <- stats::bartlett.test(resid, country_factor)
    
    ## Per-country Ljung-Box, summarized across countries
    lb_pvals <- numeric(n_sub)
    for(i in 1:n_sub){
      block <- ((i-1)*(T_dim-1)+1):(i*(T_dim-1))
      lb_pvals[i] <- stats::Box.test(resid[block], lag = lb_lag,
                                     type = "Ljung-Box")$p.value
    }
    
    ## Pesaran CD test. resid is stacked country-major in blocks of
    ## T_dim-1 (block i = ((i-1)*(T_dim-1)+1):(i*(T_dim-1))), so reshaping
    ## column-major into a (T_dim-1) x n_sub matrix puts country i's whole
    ## residual series into column i -- exactly matching that block layout.
    resid_mat <- matrix(resid, nrow = T_dim - 1, ncol = n_sub)
    R_resid <- cor(resid_mat)
    upper_sum <- sum(R_resid[upper.tri(R_resid)])
    CD_stat <- sqrt(2 * (T_dim - 1) / (n_sub * (n_sub - 1))) * upper_sum
    CD_p <- 2 * pnorm(-abs(CD_stat))
    
    data.frame(
      Model                = model_name,
      BP_stat              = round(bp_stat, 3),
      BP_p                 = round(bp_p, 4),
      Bartlett_stat        = round(unname(bart$statistic), 3),
      Bartlett_p           = round(bart$p.value, 4),
      LjungBox_pct_reject5 = round(100 * mean(lb_pvals < 0.05), 1),
      LjungBox_mean_p      = round(mean(lb_pvals), 4),
      CD_stat              = round(CD_stat, 3),
      CD_p                 = round(CD_p, 4)
    )
  }
  
  rbind(
    diagnose_one("LS",    resid_LS),
    diagnose_one("Lasso", resid_lasso)
  )
}

diag_all    <- run_diagnostics(data_all)
diag_tax    <- run_diagnostics(data_tax)
diag_no_tax <- run_diagnostics(data_no_tax)

diag_all$Group    <- "All"
diag_tax$Group    <- "Carbon tax"
diag_no_tax$Group <- "No carbon tax"


diagnostics_table <- rbind(diag_all, diag_tax, diag_no_tax)
diagnostics_table <- diagnostics_table[, c("Group", setdiff(names(diagnostics_table), "Group"))]
print(diagnostics_table)

## ------------------------------------------------------------------------
## Excel workbook: add the remaining sheets and save the final version
## ------------------------------------------------------------------------
add_results_sheet(wb, "DK_lag_sensitivity_LS", dk_lag_sens)
add_results_sheet(wb, "Fit_stats",             fit_stats_table)
add_results_sheet(wb, "VIF",                   vif_table)
add_results_sheet(wb, "Diagnostics",           diagnostics_table)

openxlsx::saveWorkbook(wb, EXCEL_FILE, overwrite = TRUE)
cat("\nAll results saved to", EXCEL_FILE, "\n")