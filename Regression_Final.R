library(MuMIn)

RegressionDiagnosicsPlots <- function(model) {
  # setting up the plot layout
  layout_matrix <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE)
  layout(layout_matrix)
  
  # 1) histogram of residuals
  hist(
    model$residuals,
    xlab = "Residuals",
    main = paste("Model =", deparse(substitute(model)))
  )
  
  # 2) QQ plot
  car::qqPlot(
    model$residuals,
    pch = 16,
    col = adjustcolor("black", 0.7),
    xlab = "Theoretical Quantiles (Normal)",
    ylab = "Sample Quantiles (r.hat)",
    main = "Normal Q-Q Plot"
  )
  
  # 3) fitted vs residuals
  plot(
    model$fitted.values,
    model$residuals,
    pch = 16,
    col = adjustcolor("black", 0.5),
    xlab = "Fitted Values",
    ylab = "Residuals"
  )
  abline(h = 0, lty = 2, col = 'red')
  
  # 4) residuals vs time (index)
  plot(
    model$residuals,
    pch = 16,
    col = adjustcolor("black", 0.5),
    xlab = "Time",
    ylab = "Residuals"
  )
  abline(h = 0, lty = 2, col = 'red')
  
  # 5) ACF of residuals
  acf(model$residuals)
  
  par(mfrow = c(1, 1)) # reset
}

RegressionDiagnosicsTests <- function(
    model
) {
  # Testing normality and the homogeneity of the variance
  residuals.model = residuals(model) #extracting the residuals
  print(shapiro.test(residuals.model))
  n <- length(residuals.model)  # should be 183432
  
  # create 6 groups of size 30572 (last one may be slightly shorter)
  seg <- factor(rep(1:6, each = ceiling(n / 6))[1:n])
  
  fligner.test(residuals.model, seg)
  print(fligner.test(residuals.model, seg)) # Fligner's test for H0:sigma1=sigma2=...=sigma6 (note that we made 6 chunks)
  
  # Testing randomness of the residuals
  par(mfrow = c(1, 1))
  print(randtests::runs.test(residuals.model, plot = TRUE))
}

# Fitting a linear regression model with time and month as explanatory variates
tim <- time(price_ts) # Extracting time as the explanatory variate from the time series framework of data
model <- lm(price_ts ~ tim + Energy_monthly$demand_monthly_avg)
RegressionDiagnosicsPlots(model)
RegressionDiagnosicsTests(model)

#Check Seasonality
Decomposed_AirPassengers = stl(
  price_ts,
  s.window = "periodic"
)
plot(Decomposed_AirPassengers)
var_season <- var(Decomposed_AirPassengers$time.series[, "seasonal"], na.rm = TRUE)
var_trend  <- var(Decomposed_AirPassengers$time.series[, "trend"], na.rm = TRUE)
var_resid  <- var(Decomposed_AirPassengers$time.series[, "remainder"], na.rm = TRUE)

prop_season <- var_season / (var_season + var_trend + var_resid)
prop_season

#boxcox
boxcox.model = MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = TRUE)
opt.lambda = boxcox.model$x[which.max(boxcox.model$y)]

Energy_monthly$price_monthly_avg_b <- (Energy_monthly$price_monthly_avg^opt.lambda - 1) / opt.lambda
par(mfrow=c(1,2))

#difference
diff_price2 <- diff(Energy_monthly$price_monthly_avg_b, lag=12)
acf(diff_price2, lag.max = 36)

# number of seasonal lags lost
n_lost <- 12


t_aln = tim[(n_lost + 1):length(tim)]
d_aln = Energy_monthly$demand_monthly_avg[(n_lost + 1):length(Energy_monthly$demand_monthly_avg)]
y = diff_price2



#Create a data Frame for regression
#t_aln <- tim
#d_aln <- Energy_monthly$demand_monthly_avg
#y <- Energy_monthly$price_monthly_avg_bc

df <- data.frame(y, t_aln, d_aln)



## =======================
## Variable selection
## =======================
ToPlot_t <- data.frame(
  t_aln = seq(min(df$t_aln), max(df$t_aln), length.out = 500)
)


## Choose t0 values
t0_vec <- c(2014+6/12, 2018+6/12)

output <- data.frame(p = NA, Pred.SE.2014 = NA, Pred.SE.2018 = NA,
                     AdjR2 = NA, AICc = NA, BIC = NA)

for (p in 1:6) {
  model <- lm(y ~ poly(t_aln, p) + d_aln, data=df)
  
  ## Prediction SE at t0 
  Pred <- predict.lm(
    model,
    newdata = data.frame(
      t_aln = t0_vec,
      d_aln = mean(d_aln, na.rm = TRUE)
    ),
    interval = "prediction", level = 0.95
  )
  
  ## Convert 95% PI width to Pred.SE: width = 2 * t_{.975,df} * SE_pred
  df_eff <- length(y) - (p + 2)  # p poly terms + intercept? (poly adds p columns) + d_aln
  tcrit  <- qt(0.975, df = max(df_eff, 1))
  pred_se <- (Pred[, "upr"] - Pred[, "lwr"]) / (2 * tcrit)  # length 2 (500, 750)
  
  output[p, ] <- c(
    p,
    pred_se[1],                
    pred_se[2],                
    summary(model)$adj.r.squared,
    MuMIn::AICc(model),
    BIC(model)
  )
}
par(mfrow=c(1,1))
round(output, 4)

## --- APSE  ---
set.seed(123)
n <- length(df$y)
pred.indx <- sample(1:n, floor(0.2*n), replace = FALSE)  # ~20% holdout

y_tr  <- df$y[-pred.indx]
t_tr  <- t_aln[-pred.indx]
d_tr  <- d_aln[-pred.indx]

APSE <- c()
for (p in 1:6) {
  model <- lm(y_tr ~ poly(t_tr, p) + d_tr)
  pred  <- predict(model, newdata = data.frame(
    t_tr = t_aln[pred.indx],
    d_tr = d_aln[pred.indx]
  ))
  APSE[p] <- mean((pred - y[pred.indx])^2)
}
round(APSE, 2)
which.min(APSE)
plot(APSE, pch = 16, type = "b")
points(which.min(APSE), APSE[which.min(APSE)],
       col = adjustcolor("red", 0.5), cex = 1.5, pch = 16)

## --- k-fold CV  ---
set.seed(1)
Shuffle <- sample(1:n, n, replace = FALSE)
y.shuf  <- df$y[Shuffle]
t.shuf  <- t_aln[Shuffle]
d.shuf  <- d_aln[Shuffle]

number.of.folds <- 5
Size.Per.Fold   <- round(n / number.of.folds)
Folds <- vector("list", number.of.folds)
for (i in 1:number.of.folds) {
  START <- Size.Per.Fold * (i - 1) + 1
  END   <- min(Size.Per.Fold * i, n)
  Folds[[i]] <- START:END
}

cross.valid <- matrix(0, nrow = number.of.folds, ncol = 6)
for (k in 1:number.of.folds) {
  x <- t.shuf[-Folds[[k]]]
  z <- d.shuf[-Folds[[k]]]
  yy <- y.shuf[-Folds[[k]]]
  for (i in 1:6) {
    fitmod <- lm(yy ~ poly(x, i) + z)
    Pred   <- predict(fitmod, newdata = data.frame(
      x = t.shuf[Folds[[k]]],
      z = d.shuf[Folds[[k]]]
    ))
    cross.valid[k, i] <- mean((y.shuf[Folds[[k]]] - Pred)^2)
  }
}
CV.Err <- apply(cross.valid, 2, mean)
plot(CV.Err, pch = 16, type = "b",
     ylab = "Estimated C.V. error (mean APE)",
     xlab = "p : polynomial degree in time",
     main = paste0(number.of.folds, "-fold Cross-validation"))
indx <- which.min(CV.Err)
points(indx, CV.Err[indx], pch = 16, col = adjustcolor("red", 0.5), cex = 1.5)

round(data.frame(output, APSE, CV.Err = CV.Err), 4)

Best_model_1=lm(y ~ poly(t_aln, 2) + d_aln)
RegressionDiagnosicsPlots(Best_model_1)
RegressionDiagnosicsTests(Best_model_1)


## =======================
## Stepwise AIC/BIC
## =======================
# Stepwise by AIC (both directions)
library(MASS)

X.ortho<-poly(t_aln, 6, raw = FALSE)
X.ortho.df <- as.data.frame(X.ortho)
colnames(X.ortho.df) <- c("t1","t2","t3","t4","t5","t6")

y=df$y
df2 <- cbind( X.ortho.df, y, d_aln)

model2 <- lm(y ~ t1 + t2 + t3 + t4 + t5 + t6 +d_aln, data =df2)


best_AIC<-stepAIC(model2, direction = "both")
summary(best_AIC)

# Stepwise by BIC (k = log(n))
n <- nrow(df2)
best_BIC <- step(model2, k = log(n), direction = "both")
summary(best_BIC)

RegressionDiagnosicsPlots(best_AIC)
RegressionDiagnosicsTests(best_AIC)

Best_model_2 <- lm(y ~ t1 + t2  + t6 +d_aln, data =df2)
RegressionDiagnosicsPlots(Best_model_2)
RegressionDiagnosicsTests(Best_model_2)

## =======================
## Elastic Net
## =======================
library(glmnet) # calling the glmnet library

length(df$y)
length(df$t_aln)
# --- Build the design matrix from your objects y, t_aln, d_aln ---
stopifnot(length(y) == length(t_aln), length(y) == length(d_aln))

X.ortho <- poly(t_aln, 6, raw = FALSE)
X.ortho.df <- as.data.frame(X.ortho)
colnames(X.ortho.df) <- paste0("t", 1:6)

df3 <- cbind(X.ortho.df, d_aln = as.numeric(d_aln))
y_vec <- as.numeric(df$y)

# glmnet needs a matrix for x
X <- as.matrix(df3)
Y <- y_vec

# --- Try a grid of alphas and pick the best by CV error ---
set.seed(123)

alphas <- c(0, 0.25, 0.5, 0.75, 1)   # ridge -> lasso
cv_list <- list()
summary_tbl <- data.frame(
  alpha = alphas,
  cvm_min = NA_real_,           # mean CV error at lambda.min
  cvm_1se = NA_real_,           # mean CV error at lambda.1se
  lambda_min = NA_real_,
  lambda_1se = NA_real_,
  nonzero_min = NA_integer_,
  nonzero_1se = NA_integer_
)

for (i in seq_along(alphas)) {
  a <- alphas[i]
  cvfit <- cv.glmnet(
    x = X, y = Y,
    family = "gaussian",
    alpha = a,
    nfolds = 10,            # Since there is a high number of observations
    standardize = TRUE,
    intercept = TRUE        # include intercept; safer than forcing 0
  )
  cv_list[[as.character(a)]] <- cvfit
  
  # Extract metrics
  summary_tbl$cvm_min[i]     <- min(cvfit$cvm)                 # at lambda.min
  summary_tbl$lambda_min[i]  <- cvfit$lambda.min
  summary_tbl$cvm_1se[i]     <- cvfit$cvm[cvfit$lambda == cvfit$lambda.1se]
  summary_tbl$lambda_1se[i]  <- cvfit$lambda.1se
  
  # Non-zero counts
  beta_min <- coef(cvfit, s = "lambda.min")
  beta_1se <- coef(cvfit, s = "lambda.1se")
  summary_tbl$nonzero_min[i] <- sum(beta_min != 0)
  summary_tbl$nonzero_1se[i] <- sum(beta_1se != 0)
}

print(summary_tbl[order(summary_tbl$cvm_min), ], row.names = FALSE)

# --- Choose a winner ---
# pick by cvm_1se
best_idx_1se <- which.min(summary_tbl$cvm_1se)
best_alpha_1se <- summary_tbl$alpha[best_idx_1se]
best_lambda_1se <- summary_tbl$lambda_1se[best_idx_1se]

cat("Best by lambda.1se:  alpha =", best_alpha_1se, " lambda =", best_lambda_1se, "\n\n")

# If you prefer the 1SE rule (simpler model), use:
best_cv_1se <- cv_list[[as.character(best_alpha_1se)]]
beta_hat_1se <- coef(best_cv_1se, s = "lambda.1se")
nz1 <- which(beta_hat_1se != 0)
cat("\nNon-zero coefficients (lambda.1se):\n")
print(beta_hat_1se[nz1, , drop = FALSE])

# --- packages you use in the functions ---
library(car)        # for qqPlot
library(randtests)  # for runs.test

# --- assume you have X (matrix or data.frame of predictors) and y (response) used in cv.glmnet ---
# 1) Identify the active set from your 1SE solution
beta_1se <- as.matrix(coef(best_cv_1se, s = "lambda.1se"))
active_idx <- which(beta_1se[-1, , drop = FALSE] != 0)          # drop intercept row
active_vars <- rownames(beta_1se)[-1][active_idx]               # names of selected predictors

# 2) Build a data.frame for lm and refit with the chosen variables
dat <- data.frame(y = y, X)
fit_1se_lm <- lm(reformulate(active_vars, response = "y"), data = dat)

# 3) Run your diagnostics
RegressionDiagnosicsPlots(fit_1se_lm)
RegressionDiagnosicsTests(fit_1se_lm)

#===========
#Test best models to choose 1 
#============
## -------- NEW helper: build newdata rows that match a model's terms ----------
.make_newdata_for_t0 <- function(model, t0_vec, base_df, P_basis = NULL) {
  nd <- data.frame(row = seq_along(t0_vec))
  
  # Which variables does the model require?
  fml_terms <- attr(terms(model), "term.labels")
  all_needed <- all.vars(formula(model))
  
  # Handle precomputed orthogonal poly columns (t1..t6)
  tcols_all <- paste0("t", 1:6)
  tcols_used <- intersect(tcols_all, all_needed)
  
  if (length(tcols_used) > 0) {
    if (is.null(P_basis)) {
      stop("Model uses t1..t6 but P_basis is NULL. Pass P_basis = poly(t_aln, 6, raw = FALSE).")
    }
    # Predict orthogonal basis at t0_vec, then subset to the tcols the model actually uses
    X0 <- as.data.frame(predict(P_basis, newdata = t0_vec))
    colnames(X0) <- tcols_all[seq_len(ncol(X0))]
    nd <- cbind(nd, X0[, tcols_used, drop = FALSE])
  }
  
  # If the model uses t_aln directly (e.g., poly(t_aln, p) inside formula), just provide t_aln
  if ("t_aln" %in% all_needed) {
    nd$t_aln <- t0_vec
  }
  
  # Provide d_aln if requested (use mean as a neutral value)
  if ("d_aln" %in% all_needed) {
    nd$d_aln <- mean(base_df$d_aln, na.rm = TRUE)
  }
  
  # Provide month if present (reference level)
  if ("month" %in% all_needed && "month" %in% names(base_df)) {
    nd$month <- factor(rep(levels(base_df$month)[1], length(t0_vec)),
                       levels = levels(base_df$month))
  }
  
  # Any other columns that appear in the model: pull from base_df as neutral values
  # - numeric → mean; factor → reference level; dummy (0/1) → column mean (works fine)
  skip <- c("(Intercept)", "y", "t_aln", "d_aln", "month", tcols_all)
  other_needed <- setdiff(all_needed, skip)
  for (v in other_needed) {
    if (!v %in% names(base_df)) next
    if (is.factor(base_df[[v]])) {
      nd[[v]] <- factor(rep(levels(base_df[[v]])[1], length(t0_vec)),
                        levels = levels(base_df[[v]]))
    } else {
      nd[[v]] <- mean(base_df[[v]], na.rm = TRUE)
    }
  }
  
  nd$row <- NULL
  nd
}

evaluate_model <- function(model, df_used, t0_vec, P_basis = NULL, kfold = 5, seed = 123) {
  # 1) Pred.SE at t0 values
  nd_t0 <- .make_newdata_for_t0(model, t0_vec, base_df = df_used, P_basis = P_basis)
  Pred <- predict(model, newdata = nd_t0, interval = "prediction", level = 0.95)
  tcrit <- qt(0.975, df = df.residual(model))
  pred_se <- (Pred[, "upr"] - Pred[, "lwr"]) / (2 * tcrit)  # length = length(t0_vec)
  
  # 2) Information criteria & AdjR2
  adjr2 <- summary(model)$adj.r.squared
  aicc  <- MuMIn::AICc(model)
  bicv  <- BIC(model)
  
  # 3) APSE (20% holdout, same formula)
  set.seed(seed)
  n <- nrow(df_used)
  hold <- sample(seq_len(n), floor(0.2 * n), replace = FALSE)
  form <- formula(model)
  fit_tr <- lm(form, data = df_used[-hold, , drop = FALSE])
  pred_te <- predict(fit_tr, newdata = df_used[hold, , drop = FALSE])
  APSE <- mean((df_used$y[hold] - pred_te)^2, na.rm = TRUE)
  
  # 4) k-fold CV error (refit same formula per fold)
  set.seed(seed + 1)
  idx <- sample(seq_len(n), n, replace = FALSE)
  folds <- split(idx, cut(seq_along(idx), breaks = kfold, labels = FALSE))
  err <- numeric(kfold)
  for (i in seq_along(folds)) {
    te <- folds[[i]]
    tr <- setdiff(idx, te)
    fit_i <- lm(form, data = df_used[tr, , drop = FALSE])
    pred_i <- predict(fit_i, newdata = df_used[te, , drop = FALSE])
    err[i] <- mean((df_used$y[te] - pred_i)^2, na.rm = TRUE)
  }
  CV.Err <- mean(err)
  
  list(
    Pred.SE = as.numeric(pred_se),
    AdjR2 = adjr2, AICc = aicc, BIC = bicv, APSE = APSE, CV.Err = CV.Err
  )
}


# keep this
P_basis <- poly(t_aln, 6, raw = FALSE)

res_best <- evaluate_model(
  model   = Best_model_2,
  df_used = df2,
  t0_vec  = t0_vec,
  P_basis = P_basis
)

# PASS P_basis HERE TOO (if fit_1se_lm uses any of t1..t6)
res_lm1se <- evaluate_model(
  model   = fit_1se_lm,
  df_used = dat,
  t0_vec  = t0_vec,
  P_basis = P_basis   # <-- was NULL before; provide the same basis
)


comp <- data.frame(
  Model = c("Best_model_2", "ENet_1SE_refit"),
  Pred.SE.2014 = c(res_best$Pred.SE[1],  res_lm1se$Pred.SE[1]),
  Pred.SE.2018 = c(res_best$Pred.SE[2],  res_lm1se$Pred.SE[2]),
  AdjR2        = c(res_best$AdjR2,       res_lm1se$AdjR2),
  AICc         = c(res_best$AICc,        res_lm1se$AICc),
  BIC          = c(res_best$BIC,         res_lm1se$BIC),
  APSE         = c(res_best$APSE,        res_lm1se$APSE),
  CV.Err       = c(res_best$CV.Err,      res_lm1se$CV.Err)
)


print(cbind(
  Model = comp$Model,
  round(comp[ , setdiff(names(comp), "Model")], 4)
))

#==========
#Forecast
#==========
## ---- Forecast to 2026 with 95% PI ----

# 0) prerequisites we already have:
#   P_basis, t_aln, d_aln, df2, Best_model_1_bc

# Step size (monthly if your series is monthly)
step <- if (exists("price_ts")) 1/frequency(price_ts) else median(diff(t_aln))

# Build a monthly grid from the period AFTER your last observation up to Dec-2026
last_t   <- max(t_aln)
end_t    <- 2026 + 11/12   # 2026-Dec in decimal year
t_future <- seq(from = last_t + step, to = end_t, by = step)

# Create the orthogonal polynomial columns at future times
P_fut <- as.data.frame(predict(P_basis, newdata = t_future))
colnames(P_fut) <- paste0("t", 1:ncol(P_fut))

# If you have your own future demand forecast, put it in d_future (same length as t_future)
# d_future <- <your vector here>
d_future <- rep(mean(d_aln, na.rm = TRUE), length(t_future))   # default: mean demand

# Assemble the newdata for THIS model (only the columns it uses)
newdata_bc <- data.frame(
  t1 = P_fut$t1,
  t2 = P_fut$t2,
  t6 = P_fut$t6,
  d_aln = d_future
)

# Forecast with 95% prediction intervals
fc_bc <- predict(Best_model_1_bc, newdata = newdata_bc, interval = "prediction", level = 0.95)

# ---- Plot forecast (fitted history + future with PI) ----
op <- par(mar = c(4,4,2,1))
# history on original scale y vs t_aln
plot(t_aln, df2$y, type = "l", xlab = "Time", ylab = "y", main = "Forecast to 2026 with 95% PI")
# fitted line over history (optional but nice)
lines(t_aln, fitted(Best_model_1_bc), col = "blue")

# future forecast & intervals
lines(t_future, fc_bc[, "fit"], lwd = 2,  col = "red")
lines(t_future, fc_bc[, "lwr"], lty = 2,  col = "red")
lines(t_future, fc_bc[, "upr"], lty = 2, col = "red")
legend("topleft", bty = "n",
       legend = c("Observed", "In-sample fit", "Forecast", "95% PI"),
       lty = c(1,1,1,2), col = c("black","blue","red","black"), lwd = c(1,1,2,1))
par(op)

# ---- Show and save forecast numbers ----
Forecast_tbl <- data.frame(
  Time = round(t_future, 4),      # or format as year-month
  Predicted = fc_bc[, "fit"],
  Lower_95 = fc_bc[, "lwr"],
  Upper_95 = fc_bc[, "upr"]
)

# Print first few lines in console
print(head(Forecast_tbl, 20))

