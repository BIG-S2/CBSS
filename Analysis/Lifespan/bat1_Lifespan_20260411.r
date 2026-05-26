#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/Lifespan
#module load r/4.4.0;R
iii=1
nets <- c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Sub","Vis1","Vis2","VMN")
paraT=expand.grid(length(nets):1,1:length(nets))
r00=paraT[iii,1]
c00=paraT[iii,2]
library(data.table)
source('gamm_hacks.r')
library(itsadug)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")
source("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R")
require(mgcv)	
S=fread('NetSC_combined.csv',data.table=F)
datalist=c("hcp","hcpa","pnc","ping","ukb","adnigo2")
indd=which(S$study%in%datalist)
indd1=which(S[,7]=='adnigo2')
#Remove MCI, AD subjects from ADNI
S1=S[indd1,];S1[,'ID']=as.numeric(unlist(lapply(strsplit(S1[,'ID'],'_'),'[',3)))
Diag=fread('/proj/tengfei/ADNI/clinical/adnimerge_20200312.csv',data.table=F)
Diag=Diag[,c('RID','EXAMDATE','DX')];Diag=Diag[!is.na(Diag[,3]),];Diag=Diag[!is.na(Diag[,2]),]
Judge=rep(NA,dim(S1)[1])
for(ii in 1:length(Judge))
{
indtemp=which(Diag[,'RID']==S1[ii,'ID']);
if(length(indtemp)>0){
Days=abs(as.numeric(as.Date(Diag[indtemp,'EXAMDATE'],format='%Y-%m-%d')-as.Date(S1[ii,'ID0'],format='%Y%m%d')))
indtemp=indtemp[which.min(Days)]; Judge[ii]=Diag[indtemp,'DX']}
}
Judge[is.na(Judge)]='Dementia'
S=S[-indd1[Judge!='CN'],]
S=S[!duplicated(S[,'ID']),]
ind00=which(S$age<90) #15 subjects excluded
S=S[ind00,];
colors <- c("purple", "orange","pink","brown","black","grey","red", "blue")
S0=S
############Best####################################################
pair_names <- colnames(S)[8:73]   # length 66
ii_vec     <- 1:66                # your ii index corresponds to pair_names[ii]
idx_mat <- matrix(NA_integer_, nrow=length(nets), ncol=length(nets),dimnames=list(nets, nets))
for (ii in ii_vec) {
  ab <- strsplit(pair_names[ii], "_", fixed=TRUE)[[1]]
  a <- ab[1]; b <- ab[2]
  ra <- match(a, nets); cb <- match(b, nets)
  if (!is.na(ra) && !is.na(cb)) {
    idx_mat[ra, cb] <- ii
    idx_mat[cb, ra] <- ii   # make symmetric so both triangles are filled
  }
}
colors <- c("brown", "grey","pink","purple", "orange","red","blue","black")
##############
library(gamlss)
library(gamlss.dist)
library(MASS)
# -----------------------------
# 1. extract selected powers
# -----------------------------
extract_fp_powers <- function(fit, param = c("mu", "sigma")) {
  param <- match.arg(param)
  obj <- switch(
    param,
    mu = fit$mu.coefSmo,
    sigma = fit$sigma.coefSmo
  )

  if (is.null(obj) || length(obj) == 0 || is.null(obj[[1]]$power)) {
    stop(sprintf("No fp powers found for %s", param))
  }

  obj[[1]]$power
}

# -----------------------------
# 2. build explicit FP basis
#    Handles repeated powers:
#    p, p, p -> x^p, x^p log(x), x^p log(x)^2
#    for p=0 -> log(x), log(x)^2, log(x)^3
# -----------------------------
make_fp_basis <- function(x, powers, prefix = "fp") {
  if (any(x <= 0, na.rm = TRUE)) {
    stop("Age must be > 0 for fractional polynomial terms")
  }

  out <- vector("list", length(powers))
  names(out) <- paste0(prefix, seq_along(powers))

  seen <- numeric(0)

  for (j in seq_along(powers)) {
    p <- powers[j]

    occ <- sum(seen == p) + 1
    seen <- c(seen, p)

    if (p == 0) {
      # repeated 0 powers: log(x), log(x)^2, log(x)^3, ...
      term <- log(x)^occ
    } else {
      # repeated nonzero powers: x^p, x^p log(x), x^p log(x)^2, ...
      term <- x^p * (log(x)^(occ - 1))
    }

    out[[j]] <- term
  }

  as.data.frame(out)
}
make_formula_with_interactions <- function(response, sex_var, study_var, basis_names) {
  rhs <- paste0(
    sex_var,
    " * (", paste(basis_names, collapse = " + "), ") + ",
    study_var
  )
  as.formula(paste(response, "~", rhs))
}

refit_jsu_explicit_fp <- function(data,
                                  response = "y",
                                  age_var = "age",
                                  sex_var = "sex",
                                  study_var = "study1",
                                  fit_fp,
                                  family_obj = JSU(),
                                  trace = FALSE,
                                  data_name = NULL) {

  dat <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.null(rownames(data))) {
    rownames(dat) <- rownames(data)
  }

  stopifnot(response %in% names(dat))
  stopifnot(age_var %in% names(dat))
  stopifnot(sex_var %in% names(dat))
  stopifnot(study_var %in% names(dat))

  if (!is.factor(dat[[sex_var]]))   dat[[sex_var]]   <- as.factor(dat[[sex_var]])
  if (!is.factor(dat[[study_var]])) dat[[study_var]] <- as.factor(dat[[study_var]])

  # capture the actual family expression from the function call, e.g. JSU()
  family_call <- substitute(family_obj)

  mu_powers <- extract_fp_powers(fit_fp, "mu")
  sigma_powers <- extract_fp_powers(fit_fp, "sigma")

  mu_basis <- make_fp_basis(dat[[age_var]], mu_powers, prefix = "mu_fp")
  sigma_basis <- make_fp_basis(dat[[age_var]], sigma_powers, prefix = "sg_fp")

  rownames(mu_basis) <- rownames(dat)
  rownames(sigma_basis) <- rownames(dat)

  for (nm in names(mu_basis)) dat[[nm]] <- mu_basis[[nm]]
  for (nm in names(sigma_basis)) dat[[nm]] <- sigma_basis[[nm]]

  mu_names <- names(mu_basis)
  sg_names <- names(sigma_basis)

  mu_formula <- as.formula(
    paste(response, "~", sex_var, "* (", paste(mu_names, collapse = " + "), ") +", study_var)
  )

  sigma_formula <- as.formula(
    paste("~", sex_var, "* (", paste(sg_names, collapse = " + "), ") +", study_var)
  )

  if (is.null(data_name)) {
    data_name <- paste0("S_augmented_fit_", sample.int(1e6, 1))
  }

  assign(data_name, dat, envir = .GlobalEnv)

  fit_explicit <- gamlss(
    formula = mu_formula,
    sigma.formula = sigma_formula,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    family = family_obj,
    data = get(data_name, envir = .GlobalEnv),
    trace = trace
  )

  # overwrite stored call so predict() can find them later
  fit_explicit$call$data <- as.name(data_name)
  fit_explicit$call$family <- family_call

  # attach metadata for augmenting newdata later
  fit_explicit$fp_meta <- list(
    age_var = age_var,
    sex_var = sex_var,
    study_var = study_var,
    mu_powers = mu_powers,
    sigma_powers = sigma_powers,
    mu_basis_names = mu_names,
    sigma_basis_names = sg_names,
    data_name = data_name
  )

  fit_explicit
}
augment_fp_newdata <- function(newdata, fit) {
  nd <- as.data.frame(newdata, stringsAsFactors = FALSE)

  meta <- fit$fp_meta
  if (is.null(meta)) stop("fit$fp_meta is missing")

  mu_basis <- make_fp_basis(nd[[meta$age_var]], meta$mu_powers, prefix = "mu_fp")
  sigma_basis <- make_fp_basis(nd[[meta$age_var]], meta$sigma_powers, prefix = "sg_fp")

  for (nm in names(mu_basis)) nd[[nm]] <- mu_basis[[nm]]
  for (nm in names(sigma_basis)) nd[[nm]] <- sigma_basis[[nm]]

  nd
}

predict_jsu_quantiles_explicit_fp <- function(fit, newdata, probs = c(0.025, 0.5, 0.975)) {
  nd <- augment_fp_newdata(newdata, fit)

  mu    <- as.numeric(predict(fit, newdata = nd, what = "mu",    type = "response"))
  sigma <- as.numeric(predict(fit, newdata = nd, what = "sigma", type = "response"))
  nu    <- as.numeric(predict(fit, newdata = nd, what = "nu",    type = "response"))
  tau   <- as.numeric(predict(fit, newdata = nd, what = "tau",   type = "response"))

  out <- sapply(probs, function(p) qJSU(p, mu = mu, sigma = sigma, nu = nu, tau = tau))
  out <- as.data.frame(out)
  colnames(out) <- paste0("q", probs)
  out
}

######################
#Bootstrap version
######################
refit_same_explicit_jsu <- function(data, fit, trace = FALSE, data_name = NULL) {
  dat <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.factor(dat[[fit$fp_meta$sex_var]])) {
    dat[[fit$fp_meta$sex_var]] <- factor(dat[[fit$fp_meta$sex_var]])
  }
  if (!is.factor(dat[[fit$fp_meta$study_var]])) {
    dat[[fit$fp_meta$study_var]] <- factor(dat[[fit$fp_meta$study_var]])
  }

  # keep factor levels aligned with original training data
  dat_train <- get(fit$fp_meta$data_name, envir = .GlobalEnv)
  dat[[fit$fp_meta$sex_var]] <- factor(
    dat[[fit$fp_meta$sex_var]],
    levels = levels(dat_train[[fit$fp_meta$sex_var]])
  )
  dat[[fit$fp_meta$study_var]] <- factor(
    dat[[fit$fp_meta$study_var]],
    levels = levels(dat_train[[fit$fp_meta$study_var]])
  )

  # rebuild the SAME explicit basis using stored powers
  mu_basis <- make_fp_basis(dat[[fit$fp_meta$age_var]], fit$fp_meta$mu_powers, prefix = "mu_fp")
  sigma_basis <- make_fp_basis(dat[[fit$fp_meta$age_var]], fit$fp_meta$sigma_powers, prefix = "sg_fp")

  for (nm in names(mu_basis)) dat[[nm]] <- mu_basis[[nm]]
  for (nm in names(sigma_basis)) dat[[nm]] <- sigma_basis[[nm]]

  if (is.null(data_name)) {
    data_name <- paste0("S_boot_aug_", sample.int(1e6, 1))
  }
  assign(data_name, dat, envir = .GlobalEnv)

  fitb <- gamlss(
    formula = fit$mu.formula,
    sigma.formula = fit$sigma.formula,
    nu.formula = fit$nu.formula,
    tau.formula = fit$tau.formula,
    family = JSU(),
    data = get(data_name, envir = .GlobalEnv),
    trace = trace
  )

  fitb$call$data <- as.name(data_name)
  fitb$call$family <- quote(JSU())
  fitb$fp_meta <- fit$fp_meta
  fitb$fp_meta$data_name <- data_name
  fitb
}
peak_jsu<-function(fit,study="hcpa",p=0.5,age_grid=NULL,nsim=500,seed=1,hessian.fun = c("PB", "R"),conf.level = 0.95,data_boot = NULL,B_boot=200)
{
  hessian.fun <- match.arg(hessian.fun);alpha <- 1 - conf.level;set.seed(seed)
  if (is.null(fit$fp_meta)) stop("fit$fp_meta is missing")
  dat_train <- get(fit$fp_meta$data_name, envir = .GlobalEnv)
  if (is.null(age_grid)) {
    age_vals <- dat_train[[fit$fp_meta$age_var]]
    age_min <- min(age_vals[age_vals > 0], na.rm = TRUE)
    age_max <- max(age_vals, na.rm = TRUE)
    age_grid <- seq(age_min, age_max, length.out = 1001)
  }
  # prediction data
  newF <- data.frame(age = age_grid,sex = factor("Female", levels = levels(dat_train[[fit$fp_meta$sex_var]])),
    study1 = factor(study, levels = levels(dat_train[[fit$fp_meta$study_var]])))
  newM <- data.frame(age = age_grid,sex = factor("Male", levels = levels(dat_train[[fit$fp_meta$sex_var]])),
    study1 = factor(study, levels = levels(dat_train[[fit$fp_meta$study_var]])))
  curveF_hat <- predict_jsu_quantiles_explicit_fp(fit, newF, probs = p)[, 1]
  curveM_hat <- predict_jsu_quantiles_explicit_fp(fit, newM, probs = p)[, 1]
  peakF_hat <- age_grid[which.max(curveF_hat)];peakM_hat <- age_grid[which.max(curveM_hat)];diff_hat <- peakM_hat - peakF_hat
  vcov_warn <- character(0)
  V_try <- withCallingHandlers(try(vcov(fit, hessian.fun = hessian.fun), silent = TRUE),
	  warning = function(w) {vcov_warn <<- c(vcov_warn, conditionMessage(w));invokeRestart("muffleWarning")})
  bad_warn <- any(grepl("NaNs produced|diag\\(V\\) had non-positive or NA entries|non-finite result may be dubious",vcov_warn))
  use_boot <- inherits(V_try, "try-error") || bad_warn
  if (!use_boot) {V <- V_try;beta_hat <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients,fit$tau.coefficients)
    if (!is.matrix(V)) stop("vcov(fit) did not return a matrix")
    if (nrow(V) != ncol(V)) stop("vcov(fit) is not square")
    if (length(beta_hat) != nrow(V)) {
      stop(sprintf("length(beta_hat)=%d but nrow(vcov)=%d",
                   length(beta_hat), nrow(V)))
    }
    rownames(V) <- names(beta_hat)
    colnames(V) <- names(beta_hat)
    ndF <- augment_fp_newdata(newF, fit)
    ndM <- augment_fp_newdata(newM, fit)
    Xmu_F <- model.matrix(delete.response(terms(fit$mu.terms)), data = ndF, contrasts.arg = fit$contrasts)
    Xmu_M <- model.matrix(delete.response(terms(fit$mu.terms)), data = ndM, contrasts.arg = fit$contrasts)
    Xsigma_F <- model.matrix(delete.response(terms(fit$sigma.terms)), data = ndF, contrasts.arg = fit$contrasts)
    Xsigma_M <- model.matrix(delete.response(terms(fit$sigma.terms)), data = ndM, contrasts.arg = fit$contrasts)
    Xnu_F <- model.matrix(delete.response(terms(fit$nu.terms)), data = ndF, contrasts.arg = fit$contrasts)
    Xnu_M <- model.matrix(delete.response(terms(fit$nu.terms)), data = ndM, contrasts.arg = fit$contrasts)
    Xtau_F <- model.matrix(delete.response(terms(fit$tau.terms)), data = ndF, contrasts.arg = fit$contrasts)
    Xtau_M <- model.matrix(delete.response(terms(fit$tau.terms)), data = ndM, contrasts.arg = fit$contrasts)
    mu_n    <- names(fit$mu.coefficients)
    sigma_n <- names(fit$sigma.coefficients)
    nu_n    <- names(fit$nu.coefficients)
    tau_n   <- names(fit$tau.coefficients)
    mu_linkinv    <- make.link.gamlss(fit$mu.link)$linkinv
    sigma_linkinv <- make.link.gamlss(fit$sigma.link)$linkinv
    nu_linkinv    <- make.link.gamlss(fit$nu.link)$linkinv
    tau_linkinv   <- make.link.gamlss(fit$tau.link)$linkinv
    sim_beta <- MASS::mvrnorm(nsim, mu = beta_hat, Sigma = V)
    peakF_sim <- rep(NA_real_, nsim)
    peakM_sim <- rep(NA_real_, nsim)
    for (b in seq_len(nsim)) {
      bb <- sim_beta[b, ]
      eta_mu_F    <- as.numeric(Xmu_F[, mu_n, drop = FALSE] %*% bb[mu_n])
      eta_sigma_F <- as.numeric(Xsigma_F[, sigma_n, drop = FALSE] %*% bb[sigma_n])
      eta_nu_F    <- as.numeric(Xnu_F[, nu_n, drop = FALSE] %*% bb[nu_n])
      eta_tau_F   <- as.numeric(Xtau_F[, tau_n, drop = FALSE] %*% bb[tau_n])
      mu_F    <- mu_linkinv(eta_mu_F)
      sigma_F <- sigma_linkinv(eta_sigma_F)
      nu_F    <- nu_linkinv(eta_nu_F)
      tau_F   <- tau_linkinv(eta_tau_F)
      curveF_b <- qJSU(p, mu = mu_F, sigma = sigma_F, nu = nu_F, tau = tau_F)
      peakF_sim[b] <- age_grid[which.max(curveF_b)]
      eta_mu_M    <- as.numeric(Xmu_M[, mu_n, drop = FALSE] %*% bb[mu_n])
      eta_sigma_M <- as.numeric(Xsigma_M[, sigma_n, drop = FALSE] %*% bb[sigma_n])
      eta_nu_M    <- as.numeric(Xnu_M[, nu_n, drop = FALSE] %*% bb[nu_n])
      eta_tau_M   <- as.numeric(Xtau_M[, tau_n, drop = FALSE] %*% bb[tau_n])
      mu_M    <- mu_linkinv(eta_mu_M)
      sigma_M <- sigma_linkinv(eta_sigma_M)
      nu_M    <- nu_linkinv(eta_nu_M)
      tau_M   <- tau_linkinv(eta_tau_M)
      curveM_b <- qJSU(p, mu = mu_M, sigma = sigma_M, nu = nu_M, tau = tau_M)
      peakM_sim[b] <- age_grid[which.max(curveM_b)]
    }
    ok <- is.finite(peakF_sim) & is.finite(peakM_sim)
    peakF_sim <- peakF_sim[ok]
    peakM_sim <- peakM_sim[ok]
    diff_sim  <- peakM_sim - peakF_sim
  }
# -------- bootstrap fallback --------
	if (use_boot) {
	  if (is.null(data_boot)) {
		stop("vcov() failed. Please provide data_boot for bootstrap fallback.")}
	  peakF_sim <- rep(NA_real_, B_boot)
	  peakM_sim <- rep(NA_real_, B_boot)
	  print('Boots: ');
	  for (b in seq_len(B_boot)) {
	    print(B_boot-b)
		idx <- sample(seq_len(nrow(data_boot)), replace = TRUE)
		db <- data_boot[idx, , drop = FALSE]
		fitb <- try(
		  refit_same_explicit_jsu(
			data = db,
			fit = fit,
			trace = FALSE,
			data_name = paste0("S_boot_aug_", b, "_", sample.int(1e6, 1))
		  ),
		  silent = TRUE
		)
		if (inherits(fitb, "try-error")) next
		predF_b <- try(
		  predict_jsu_quantiles_explicit_fp(fitb, newF, probs = p)[, 1],
		  silent = TRUE
		)
		predM_b <- try(
		  predict_jsu_quantiles_explicit_fp(fitb, newM, probs = p)[, 1],
		  silent = TRUE
		)
		if (inherits(predF_b, "try-error") || inherits(predM_b, "try-error")) next
		peakF_sim[b] <- age_grid[which.max(predF_b)]
		peakM_sim[b] <- age_grid[which.max(predM_b)]
	  }
	  ok <- is.finite(peakF_sim) & is.finite(peakM_sim)
	  peakF_sim <- peakF_sim[ok]
	  peakM_sim <- peakM_sim[ok]
	  diff_sim  <- peakM_sim - peakF_sim
	}
  if (length(diff_sim) < 10) {
    warning("Very few valid simulated peaks; CI/p-value may be unstable")
  }
  p_diff <- 2 * min(
    (sum(diff_sim >= 0) + 1) / (length(diff_sim) + 1),
    (sum(diff_sim <= 0) + 1) / (length(diff_sim) + 1)
  )
  p_diff <- min(p_diff, 1)
  list(
    conf.level = conf.level,
    p_quantile = p,
    female = c(
      peak = peakF_hat,
      lwr = unname(quantile(peakF_sim, alpha / 2, na.rm = TRUE)),
      upr = unname(quantile(peakF_sim, 1 - alpha / 2, na.rm = TRUE))
    ),
    male = c(
      peak = peakM_hat,
      lwr = unname(quantile(peakM_sim, alpha / 2, na.rm = TRUE)),
      upr = unname(quantile(peakM_sim, 1 - alpha / 2, na.rm = TRUE))
    ),
    diff_male_minus_female = c(
      estimate = diff_hat,
      lwr = unname(quantile(diff_sim, alpha / 2, na.rm = TRUE)),
      upr = unname(quantile(diff_sim, 1 - alpha / 2, na.rm = TRUE))
    ),
    p_value_diff = p_diff,
    nsim_used = length(diff_sim),
    method = if (use_boot) "bootstrap_fallback" else "vcov_simulation"
  )
}
##############
# ---------- main ----------
age_min <- max(min(S$age[S$age > 0], na.rm = TRUE), 1e-6);age_grid <- seq(age_min, 90, length.out = 181)
plot_one_pair <- function(ii) {
  S <- S0
  S$y <- S[, ii + 7]
  eps <- 1e-4
  S$y <- pmin(pmax(S$y, eps), 1 - eps)
  S$study1 <- relevel(factor(S$study), ref = "hcpa")
  S$sex    <- factor(S$male, levels = c(0, 1), labels = c("Female", "Male"))
  fit_jsu0 <- gamlss(y ~ fp(age, npoly = 3)+sex + study1,sigma.formula = ~ fp(age, npoly = 3)+sex + study1,
    nu.formula = ~ 1,tau.formula = ~ 1,family = JSU(),data = S,trace = FALSE)
  fit_jsu1 <- gamlss(y ~ fp(age, npoly = 2)+sex + study1,sigma.formula = ~ fp(age, npoly = 2)+sex + study1,
    nu.formula = ~ 1,tau.formula = ~ 1,family = JSU(),data = S,trace = FALSE)
  BIC0=BIC(fit_jsu0);BIC1=BIC(fit_jsu1);fit_jsu=fit_jsu0;if(BIC0>BIC1)fit_jsu=fit_jsu1
  ##fit_jsu_final0 <- refit_jsu_explicit_fp(data = S,response = "y",age_var = "age",sex_var = "sex",study_var = "study1",
  ##	  fit_fp = fit_jsu,family_obj = JSU(),trace = FALSE);fit_jsu_final=fit_jsu_final0$fit; dat=fit_jsu_final0$data
  fit_jsu_final <- refit_jsu_explicit_fp(data = S,response = "y", age_var = "age",sex_var = "sex",study_var = "study1",
  fit_fp = fit_jsu,family_obj = JSU(),trace = FALSE,data_name = "S_augmented_fit")
  fit_nosex <- gamlss(formula = y ~ mu_fp1 + mu_fp2 + study1, sigma.formula = ~ sg_fp1 + sg_fp2 + study1,nu.formula = ~ 1,tau.formula = ~ 1,family = JSU(),data = S_augmented_fit,trace = FALSE)
  compare0<-capture.output(LR.test(fit_nosex,fit_jsu_final))
  compare0 <- as.numeric(sub(".*p-value=\\s*", "",compare0[grepl("p-value", compare0)]))
  peak_out <- peak_jsu(fit_jsu_final,study = "hcpa",p = 0.5,nsim = 500,hessian.fun = "PB",conf.level = 0.95,data_boot = S,B_boot=200)
  # -------- shift observed points to HCPA reference --------
  S_ref <- S;S_ref$study1 <- factor("hcpa", levels = levels(S$study1))
  q_obs <- predict_jsu_quantiles_explicit_fp(fit_jsu_final, S, probs = 0.5)[,1]
  q_ref <- predict_jsu_quantiles_explicit_fp(fit_jsu_final, S_ref, probs = 0.5)[,1]
  S$y1 <- S$y - (q_obs - q_ref)
  par(mgp = c(-1.0, -1.0, 0))
  plot(S$age, S$y1,ylim = c(0.15, 0.75),pch = 20, col = "white", cex = 0.08, main = "", cex.axis = 1.2, cex.main = 1.5, cex.lab = 1.5,
    xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(2, at = c(0.25, 0.7), labels = c("", ""), las = 1, tck = -0.02, cex.axis = 2.5)
  points(S$age[S$study == "hcp"],     S$y1[S$study == "hcp"],     col = adjustcolor(colors[1], alpha.f = 0.1), pch = 19, cex = 0.08)
  points(S$age[S$study == "hcpa"],    S$y1[S$study == "hcpa"],    col = adjustcolor(colors[2], alpha.f = 0.1), pch = 19, cex = 0.08)
  points(S$age[S$study == "pnc"],     S$y1[S$study == "pnc"],     col = adjustcolor(colors[4], alpha.f = 0.1), pch = 19, cex = 0.08)
  points(S$age[S$study == "ping"],    S$y1[S$study == "ping"],    col = adjustcolor(colors[5], alpha.f = 0.1), pch = 19, cex = 0.08)
  points(S$age[S$study == "ukb"],     S$y1[S$study == "ukb"],     col = adjustcolor(colors[6], alpha.f = 0.1), pch = 19, cex = 0.08)
  points(S$age[S$study == "adnigo2"], S$y1[S$study == "adnigo2"], col = adjustcolor(colors[7], alpha.f = 0.1), pch = 19, cex = 0.08)
  # -------- prediction grid --------
  newF <- data.frame(age = age_grid,sex = factor("Female", levels = levels(S$sex)),study1 = factor("hcpa", levels = levels(S$study1)))
  newM <- data.frame(age = age_grid,sex = factor("Male", levels = levels(S$sex)),study1 = factor("hcpa", levels = levels(S$study1)))
  predF <- predict_jsu_quantiles_explicit_fp(fit_jsu_final, newF, probs = c(0.025, 0.5, 0.975))
  predM <- predict_jsu_quantiles_explicit_fp(fit_jsu_final, newM, probs = c(0.025, 0.5, 0.975))
  # -------- median trajectory + 95% centile bands --------
  lines(age_grid, predF[, "q0.5"],   col = "red",  lwd = 2)
  lines(age_grid, predF[, "q0.025"], col = "red",  lty = 2)
  lines(age_grid, predF[, "q0.975"], col = "red",  lty = 2)
  lines(age_grid, predM[, "q0.5"],   col = "blue", lwd = 2)
  lines(age_grid, predM[, "q0.025"], col = "blue", lty = 2)
  lines(age_grid, predM[, "q0.975"], col = "blue", lty = 2)
  # -------- peak ages --------
  peakF <- which.max(predF[, "q0.5"]);peakM <- which.max(predM[, "q0.5"])
  ageF <- age_grid[peakF];yF   <- predF[peakF, "q0.5"];ageM <- age_grid[peakM];yM   <- predM[peakM, "q0.5"]
  points(ageF, yF, pch = 16, col = "red");ageF1=peak_out$female['peak'];ageM1=peak_out$male['peak']
  #text(ageF1, yF, labels = paste0(round(ageF1, 1),' [',as.numeric(round(peak_out$female['lwr'],1)),',',
  #as.numeric(round(peak_out$female['upr'],1))), col = "red", pos = 3, cex = 3.0)
  text(ageF1, yF, labels = paste0(round(ageF1, 1)), col = "red", pos = 3, cex = 3.0)
  points(ageM, yM, pch = 16, col = "blue")
  #text(ageM1, yM - 0.2, labels = paste0(round(ageM1, 1),' [',as.numeric(round(peak_out$male['lwr'],1)),',',
  #as.numeric(round(peak_out$male['upr'],1))), col = "blue", pos = 3, cex = 3.0)
  text(ageM1, yM - 0.2, labels = paste0(round(ageM1, 1)), col = "blue", pos = 3, cex = 3.0)
  peak.stat=c(peak_out$female,peak_out$male,peak_out$diff_male_minus_female,peak_out$p_value_diff,peak_out$nsim_used)
  names(peak.stat)=c('Peak_F','Peak_F_95%L','Peak_F_95%U','Peak_M','Peak_M_95%L','Peak_M_95%U','FM_Dif','FM_Dif_95%L','FM_Dif_95%U','FM_Dif_p','SimN')
  invisible(list(fit = fit_jsu_final,predF = predF,predM = predM,peak.stat=peak.stat,SexTest=compare0))
}
ModelOut='GAMLSS/Model_saved/'
system(paste0('mkdir -p ',ModelOut))
Name0=c('Peak_F','Peak_F_95%L','Peak_F_95%U','Peak_M','Peak_M_95%L','Peak_M_95%U','FM_Dif','FM_Dif_95%L','FM_Dif_95%U','FM_Dif_p','SimN')
MT=matrix(NA, sum(!is.na(idx_mat))/2+sum(!is.na(diag(idx_mat)))/2,length(Name0)+1)
colnames(MT)=c(Name0,'SexDif.pval');indd=1;rownames(MT)=paste0(1:dim(MT)[1])
print(c(r00,length(nets)-c00))
ii <- idx_mat[r00, c00]
if (!is.na(ii)) {
  filetemp=paste0(ModelOut,'/PeakStat_',rownames(idx_mat)[r00],'_',colnames(idx_mat)[c00],'.txt')
  if(file.exists(filetemp)){MT[indd,]=read.table(filetemp,head=T)[[1]];indd=indd+1;}
  if(!file.exists(filetemp)){
  Res=plot_one_pair(ii)
  saveRDS(Res$fit_jsu_final,file=paste0(ModelOut,'/Model_',rownames(idx_mat)[r00],'_',colnames(idx_mat)[c00],'.rds'))
  temp=cbind(Res$predF,Res$predM);colnames(temp)[1:3]=paste0('F_',colnames(temp)[1:3]);colnames(temp)[4:6]=paste0('M_',colnames(temp)[4:6])
  temp=cbind(age_grid,temp)
  write.csv(temp,file=paste0(ModelOut,'/Fitted_',rownames(idx_mat)[r00],'_',colnames(idx_mat)[c00],'.csv'),quote=F,row.names=F)
  MT[indd,]=c(Res$peak.stat,Res$SexTest);rownames(MT)[indd]=paste0(rownames(idx_mat)[r00],'_',colnames(idx_mat)[c00])
  write.csv(MT[indd,],file=paste0(ModelOut,'/PeakStat_',rownames(idx_mat)[r00],'_',colnames(idx_mat)[c00],'.txt'),quote=F,row.names=F)
  indd=indd+1}
}
