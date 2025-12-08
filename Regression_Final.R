library(glmnet)
library(MASS)
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


# Fitting a linear regression model with time
tim <- time(price_ts) # Extracting time as the explanatory variate from the time series framework of data
model <- lm(price_ts ~ tim + Energy_monthly$demand_monthly_avg)
RegressionDiagnosicsPlots(model)
RegressionDiagnosicsTests(model)


par(mfrow=c(1,2))
price_ts<- ts(Energy_monthly$price_monthly_avg,
               start = c(start_y, start_m),
               frequency = 12)
# Plot the ts
plot(price_ts, ylab = "Avg HOEP (¢/kWh)", xlab = "Time",
     main = "Monthly Average Price (ts)")
acf(Energy_monthly$price_monthly_avg)
par(mfrow=c(1,1))

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


#Create a data Frame for regression
t_aln <- tim
d_aln <- Energy_monthly$demand_monthly_avg
y <- Energy_monthly$price_monthly_avg

df <- data.frame(y, t_aln, d_aln)
## =======================
## Variable selection
## =======================
## Grid for plotting 
ToPlot_t <- data.frame(
  t_aln = seq(min(df$t_aln), max(df$t_aln), length.out = 500)
)

## Choose t0 values
t0_vec <- c(2014 + 7/12, 2018 + 7/12)

output <- data.frame(t = NA, Pred.SE.2014 = NA, Pred.SE.2018 = NA,
                     AdjR2 = NA, AICc = NA, BIC = NA)

## Model selection
for (p in 1:6) {
  model <- lm(y ~ poly(t_aln, p) + d_aln +
                sin(4*pi*t_aln) + cos(4*pi*t_aln),
              data = df)
  
  ## Prediction SE at t0 (supply t_aln & d_aln; sin/cos are computed from t_aln)
  Pred <- predict.lm(
    model,
    newdata = data.frame(
      t_aln = t0_vec,
      d_aln = mean(df$d_aln, na.rm = TRUE)
    ),
    interval = "prediction", level = 0.95
  )
  
  ## Convert 95% PI width to Pred.SE
  tcrit   <- qt(0.975, df = max(df.residual(model), 1))
  pred_se <- (Pred[, "upr"] - Pred[, "lwr"]) / (2 * tcrit)
  
  output[p, ] <- c(
    p,
    pred_se[1],
    pred_se[2],
    summary(model)$adj.r.squared,
    MuMIn::AICc(model),
    BIC(model)
  )
}

par(mfrow = c(1,1))
round(output, 4)

## APSE (20% holdout) 
set.seed(123)
n <- nrow(df)
pred.indx <- sample(1:n, floor(0.2*n), replace = FALSE)

y_tr <- df$y[-pred.indx]
t_tr <- df$t_aln[-pred.indx]
d_tr <- df$d_aln[-pred.indx]

APSE <- numeric(6)
for (p in 1:6) {
  model <- lm(y_tr ~ poly(t_tr, p) + d_tr +
                sin(4*pi*t_tr) + cos(4*pi*t_tr))
  
  pred <- predict(model, newdata = data.frame(
    t_tr = df$t_aln[pred.indx],
    d_tr = df$d_aln[pred.indx]
  ))
  APSE[p] <- mean((pred - df$y[pred.indx])^2)
}

round(APSE, 2)
which.min(APSE)
plot(APSE, pch = 16, type = "b")
points(which.min(APSE), APSE[which.min(APSE)],
       col = adjustcolor("red", 0.5), cex = 1.5, pch = 16)

## k-fold CV with harmonics
set.seed(1)
Shuffle <- sample(1:n, n, replace = FALSE)
y.shuf  <- df$y[Shuffle]
t.shuf  <- df$t_aln[Shuffle]
d.shuf  <- df$d_aln[Shuffle]

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
  x  <- t.shuf[-Folds[[k]]]
  z  <- d.shuf[-Folds[[k]]]
  yy <- y.shuf[-Folds[[k]]]
  
  for (i in 1:6) {
    fitmod <- lm(yy ~ poly(x, i) + z +
                   sin(4*pi*x) + cos(4*pi*x))
    
    Pred <- predict(fitmod, newdata = data.frame(
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

## Final best model example
Best_model_1 <- lm(y ~ poly(t_aln, 3) + d_aln +
                     sin(4*pi*t_aln) + cos(4*pi*t_aln),
                   data = df)
RegressionDiagnosicsPlots(Best_model_1)
RegressionDiagnosicsTests(Best_model_1)


#boxcox
# make sure df is a plain data.frame
df <- as.data.frame(df)

bc <- boxcox(
  y ~ poly(t_aln, 3) + d_aln +
    sin(4*pi*t_aln) + cos(4*pi*t_aln),
  data   = df,
  lambda = seq(-2, 2, 0.1),
  plotit = TRUE
)


opt.lambda = bc$x[which.max(bc$y)]

y_bc1 <- (y^opt.lambda - 1) / opt.lambda
Best_model_1_bc <- lm(y_bc1 ~ poly(t_aln, 3) + d_aln +
                     sin(4*pi*t_aln) + cos(4*pi*t_aln),
                   data = df)
summary(Best_model_1_bc)
RegressionDiagnosicsPlots(Best_model_1_bc)
RegressionDiagnosicsTests(Best_model_1_bc)

## =======================
## Stepwise AIC/BIC
## =======================
# Stepwise by AIC (both directions)

X.ortho <- poly(t_aln, 6, raw = FALSE)
X.ortho.df <- as.data.frame(X.ortho)
names(X.ortho.df) <- paste0("t", 1:6)

df2 <- data.frame(
  y      = as.numeric(y),
  t_aln  = as.numeric(t_aln),
  d_aln  = as.numeric(d_aln)
)
df2 <- cbind(df2, X.ortho.df)

## precompute seasonal harmonics (period = 6 here)
df2$sin6 <- sin(4*pi*df2$t_aln)
df2$cos6 <- cos(4*pi*df2$t_aln)


model2 <- lm(y ~ t1 + t2 + t3 + t4 + t5 + t6 + d_aln + sin6 + cos6, data = df2)

best_AIC <- stepAIC(model2, direction = "both")
summary(best_AIC)

n <- nrow(df2)
best_BIC <- step(model2, k = log(n), direction = "both")
summary(best_BIC)

Best_model_2 <- lm(y ~ t1 + t2 + t3 + d_aln + cos6, data = df2)
RegressionDiagnosicsPlots(Best_model_2)
RegressionDiagnosicsTests(Best_model_2)

#boxcox
bc <- boxcox(Best_model_2, lambda = seq(-2, 2, 0.1), plotit = TRUE)
opt.lambda <- bc$x[which.max(bc$y)]

df2$y_bc2 <- if (abs(opt.lambda) < 1e-6) log(df2$y) else (df2$y^opt.lambda - 1)/opt.lambda

Best_model_2_bc <- lm(y_bc2 ~ t1 + t2 + t3 + d_aln + cos6, data = df2)
RegressionDiagnosicsPlots(Best_model_2_bc)
RegressionDiagnosicsTests(Best_model_2_bc)

## =======================
## Elastic Net
## =======================
stopifnot(length(y) == length(t_aln), length(y) == length(d_aln))

# orthogonal polys from t_aln
P_basis <- poly(as.numeric(t_aln), 6, raw = FALSE)
X.ortho.df <- as.data.frame(P_basis)
names(X.ortho.df) <- paste0("t", 1:6)

# seasonal harmonics from t_aln 
sin6 <- sin(4*pi*as.numeric(t_aln))
cos6 <- cos(4*pi*as.numeric(t_aln))

# master df
D <- data.frame(
  y      = as.numeric(y),
  t_aln  = as.numeric(t_aln),
  d_aln  = as.numeric(d_aln),
  X.ortho.df,
  sin6 = sin6,
  cos6 = cos6
)

# drop rows with NA in any modeling columns
D <- na.omit(D)

# Elastic Net 
# predictors 
 predictor_cols <- c(paste0("t",1:6), "d_aln", "sin6", "cos6")


X <- as.matrix(D[, predictor_cols, drop = FALSE])
Y <- D$y

set.seed(123)
alphas <- c(0, 0.25, 0.5, 0.75, 1)
cv_list <- list()
summary_tbl <- data.frame(
  alpha       = alphas,
  cvm_min     = NA_real_,
  cvm_1se     = NA_real_,
  lambda_min  = NA_real_,
  lambda_1se  = NA_real_,
  nonzero_min = NA_integer_,
  nonzero_1se = NA_integer_
)

for (i in seq_along(alphas)) {
  a <- alphas[i]
  cvfit <- cv.glmnet(
    x = X, y = Y,
    family = "gaussian",
    alpha = a,
    nfolds = 10,
    standardize = TRUE,
    intercept = TRUE
  )
  cv_list[[as.character(a)]] <- cvfit
  
  summary_tbl$cvm_min[i]     <- min(cvfit$cvm)
  summary_tbl$lambda_min[i]  <- cvfit$lambda.min
  summary_tbl$cvm_1se[i]     <- cvfit$cvm[cvfit$lambda == cvfit$lambda.1se]
  summary_tbl$lambda_1se[i]  <- cvfit$lambda.1se
  
  beta_min <- coef(cvfit, s = "lambda.min")
  beta_1se <- coef(cvfit, s = "lambda.1se")
  summary_tbl$nonzero_min[i] <- sum(beta_min != 0)
  summary_tbl$nonzero_1se[i] <- sum(beta_1se != 0)
}

print(summary_tbl[order(summary_tbl$cvm_min), ], row.names = FALSE)

## choose by lambda.min instead of lambda.1se
best_idx_min    <- which.min(summary_tbl$cvm_min)
best_alpha_min  <- summary_tbl$alpha[best_idx_min]
best_cv_min     <- cv_list[[as.character(best_alpha_min)]]
best_lambda_min <- summary_tbl$lambda_min[best_idx_min]
cat("Best by lambda.min:  alpha =", best_alpha_min, " lambda =", best_lambda_min, "\n\n")

## coefficients at lambda.min
beta_min <- as.matrix(coef(best_cv_min, s = "lambda.min"))

## drop intercept row -> active vars
nz_idx <- which(beta_min[-1, , drop = FALSE] != 0)
vars_from_coeffs <- rownames(beta_min)[-1][nz_idx]
active_vars <- vars_from_coeffs



# Refit OLS on selected vars 
Best_model_3 <- lm(reformulate(active_vars, response = "y"), data = D)
RegressionDiagnosicsPlots(Best_model_3)
RegressionDiagnosicsTests(Best_model_3)

##Box–Cox on OLS model, refit on transformed y 
bc <- boxcox(Best_model_3, lambda = seq(-2, 2, 0.1), plotit = TRUE)
opt.lambda <- bc$x[which.max(bc$y)]

D$y_bc <- (D$y^opt.lambda - 1)/opt.lambda

Best_model_3_bc <- lm(reformulate(active_vars, response = "y_bc"), data = D)
summary(Best_model_3_bc)
RegressionDiagnosicsPlots(Best_model_3_bc)
RegressionDiagnosicsTests(Best_model_3_bc)

#===========
#Test best models to choose 1 
#============
.make_newdata_for_t0 <- function(model, t0_vec, base_df, P_basis = NULL) {
  nd <- data.frame(row = seq_along(t0_vec))
  
  fml_terms  <- attr(terms(model), "term.labels")
  all_needed <- all.vars(formula(model))
  
  tcols_all  <- paste0("t", 1:6)
  tcols_used <- intersect(tcols_all, all_needed)
  
  if (length(tcols_used) > 0) {
    if (is.null(P_basis)) stop("Model uses t1..t6 but P_basis is NULL.")
    X0 <- as.data.frame(predict(P_basis, newdata = t0_vec))
    colnames(X0) <- tcols_all[seq_len(ncol(X0))]
    nd <- cbind(nd, X0[, tcols_used, drop = FALSE])
  }
  
  if ("t_aln" %in% names(base_df) && "t_aln" %in% all_needed)
    nd$t_aln <- t0_vec
  
  if ("d_aln" %in% all_needed) {
    m <- match(t0_vec, base_df$t_aln)
    d_at_t0 <- base_df$d_aln[m]
    d_at_t0[is.na(d_at_t0)] <- mean(base_df$d_aln, na.rm = TRUE)
    nd$d_aln <- d_at_t0
  }
  
  ## seasonal columns
  if ("sin6" %in% all_needed) nd$sin6 <- sin(4*pi*as.numeric(t0_vec))
  if ("cos6" %in% all_needed) nd$cos6 <- cos(4*pi*as.numeric(t0_vec))
  
  # fill any other needed terms with means/baselines
  skip <- c("(Intercept)","y","y_bc2","t_aln","d_aln",tcols_all,"sin6","cos6") 
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
  form <- formula(model)
  resp <- all.vars(form)[1]  # actual response used by this model
  
  # 1) Pred.SE at t0
  nd_t0 <- .make_newdata_for_t0(model, t0_vec, base_df = df_used, P_basis = P_basis)
  Pred <- predict(model, newdata = nd_t0, interval = "prediction", level = 0.95)
  tcrit <- qt(0.975, df = df.residual(model))
  pred_se <- (Pred[, "upr"] - Pred[, "lwr"]) / (2 * tcrit)
  
  # 2) Info criteria & AdjR2
  adjr2 <- summary(model)$adj.r.squared
  aicc  <- MuMIn::AICc(model)
  bicv  <- BIC(model)
  
  # 3) APSE (20% holdout) 
  set.seed(seed)
  n <- nrow(df_used)
  hold <- sample(seq_len(n), floor(0.2 * n), replace = FALSE)
  fit_tr <- lm(form, data = df_used[-hold, , drop = FALSE])
  pred_te <- predict(fit_tr, newdata = df_used[hold, , drop = FALSE])
  APSE <- mean((df_used[[resp]][hold] - pred_te)^2, na.rm = TRUE)
  
  # 4) k-fold CV error
  set.seed(seed + 1)
  idx <- sample(seq_len(n), n, replace = FALSE)
  folds <- split(idx, cut(seq_along(idx), breaks = kfold, labels = FALSE))
  err <- numeric(kfold)
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(idx, te)
    fit_i <- lm(form, data = df_used[tr, , drop = FALSE])
    pred_i <- predict(fit_i, newdata = df_used[te, , drop = FALSE])
    err[i] <- mean((df_used[[resp]][te] - pred_i)^2, na.rm = TRUE)
  }
  CV.Err <- mean(err)
  
  list(Pred.SE = as.numeric(pred_se), AdjR2 = adjr2,
       AICc = aicc, BIC = bicv, APSE = APSE, CV.Err = CV.Err)
}


P_basis <- poly(df2$t_aln, 6, raw = FALSE)

res_best <- evaluate_model(
  model   = Best_model_2_bc,
  df_used = df2,
  t0_vec  = t0_vec,
  P_basis = P_basis
)



res_lm1se <- evaluate_model(
  model   = Best_model_3_bc,
  df_used = D,
  t0_vec  = t0_vec,
  P_basis = P_basis
)
res_lm1se


comp <- data.frame(
  Model = c("Best_model_2", "Best_model_3"),
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

## =========================
## Forecast to Dec-2026 (+ PI) from Best_model_2_bc
## =============================
# Step size from your ts() frequency
step   <- 1 / frequency(price_ts)
last_t <- max(t_aln, na.rm = TRUE)
end_t  <- 2026 + 11/12
t_future <- seq(from = last_t + step, to = end_t, by = step)

# FUTURE orthogonal-polynomial columns from SAME basis
P_fut <- as.data.frame(predict(P_basis, newdata = t_future))
colnames(P_fut) <- paste0("t", 1:ncol(P_fut))

# FUTURE d_aln scenario = mean training d_aln
mf <- model.frame(Best_model_2_bc)
fit_d <- lm(d_aln ~ sin(4*pi*t_aln) + cos(4*pi*t_aln), data = D)
d_future <- predict(fit_d, newdata = data.frame(t_aln = t_future))


# Assemble NEWDATA matching the exact terms the model uses
need       <- attr(terms(Best_model_2_bc), "term.labels")
tcols_used <- intersect(paste0("t", 1:6), need)

newdata_bc <- data.frame(row = seq_along(t_future))
if (length(tcols_used) > 0) newdata_bc <- cbind(newdata_bc, P_fut[, tcols_used, drop = FALSE])
if ("d_aln" %in% need) newdata_bc$d_aln <- d_future
if ("sin6" %in% need) newdata_bc$sin6 <- sin(4*pi*as.numeric(t_future))
if ("cos6" %in% need) newdata_bc$cos6 <- cos(4*pi*as.numeric(t_future))
newdata_bc$row <- NULL

# Sanity checks
stopifnot(all(setdiff(need, names(newdata_bc)) == character(0)))
stopifnot(nrow(newdata_bc) == length(t_future))

# Forecast on Box–Cox scale
fc_bc <- predict(Best_model_2_bc, newdata = newdata_bc, interval = "prediction", level = 0.95)

# Inverse Box–Cox helper
inv_boxcox <- function(z, lambda) {
  if (abs(lambda) < 1e-6) exp(z) else pmax(0, (lambda*z + 1))^(1/lambda)
}

# Back-transform
fit_fut <- inv_boxcox(fc_bc[, "fit"], opt.lambda)
lwr_fut <- inv_boxcox(fc_bc[, "lwr"], opt.lambda)
upr_fut <- inv_boxcox(fc_bc[, "upr"], opt.lambda)

# Historical fitted on training frame 
pred_hist_bc <- predict(Best_model_2_bc, newdata = mf)
fit_hist <- inv_boxcox(pred_hist_bc, opt.lambda)

## =============================
## Plot observed, fitted, forecast
## =============================
# Avoid 'type' name collision; coerce numeric
if (exists("type", inherits = TRUE)) rm(type)

y_obs_num   <- as.numeric(model.response(mf))         # this is y_bc scale
y_obs_bt    <- inv_boxcox(y_obs_num, opt.lambda)      # back-transformed for original scale
fit_hist_n  <- as.numeric(fit_hist)
fit_fut_n   <- as.numeric(fit_fut)
lwr_fut_n   <- as.numeric(lwr_fut)
upr_fut_n   <- as.numeric(upr_fut)

op <- par(mar = c(4,4,2,1))

# plot as before, let it draw the regular 5-year ticks
x_all <- c(as.numeric(t_aln), as.numeric(t_future))
y_all <- c(y_obs_bt, fit_hist_n, lwr_fut_n, upr_fut_n)
plot(as.numeric(t_aln), y_obs_bt, type = "l",
     xlab = "Time", ylab = "y",
     main = "Forecast to 2026 with 95% Prediction Intervals",
     xlim = range(x_all, na.rm = TRUE),
     ylim = range(y_all, na.rm = TRUE))

# fitted + forecast + PI
graphics::lines(as.numeric(t_aln),    fit_hist_n,  lwd = 2, col = "blue")
graphics::lines(as.numeric(t_future), fit_fut_n,   lwd = 2, col = "red")
polygon(c(as.numeric(t_future), rev(as.numeric(t_future))),
        c(lwr_fut_n, rev(upr_fut_n)),
        border = NA, col = adjustcolor("red", 0.2))

# add the three dotted year markers
abline(v = 2024, lty = 2, col = "darkgray", lwd = 1.5)
abline(v = 2025, lty = 2, col = "darkgray", lwd = 1.5)
abline(v = 2026, lty = 2, col = "darkgray", lwd = 1.5)


legend("topleft",
       legend = c("Observed", "Fitted (hist)", "Forecast", "95% PI"),
       lty = c(1,1,1,NA), lwd = c(1,2,2,NA),
       pch = c(NA, NA, NA, 15),
       col = c("black","blue","red", adjustcolor("red", 0.2)),
       pt.cex = 2, bty = "n")

par(op)
