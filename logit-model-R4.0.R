# logit - log(dose) model£¨R 4.0 Compatible£©
rm(list=ls())
library(MASS)
library(ggplot2)
library(effects)

# read data, ensure your data in the R working directory 
Data <- read.table("Fairfax.txt", header = TRUE, quote = "\"")

# The data should be input as following format
# Trt	D	N	R
# Fairfax	0	30	0
# Fairfax	2	48	12
# Fairfax	3	50	15
# Fairfax	5	50	31
# Fairfax	7	48	31
# Fairfax	10	59	52
# Pixley	0	359	7
# Pixley	10	70	22
# Pixley	20	70	38
# Pixley	30	50	38
# Pixley	50	50	48
# Schaefer	0	600	0
# Schaefer	2	60	15
# Schaefer	3	120	41
# Schaefer	5	60	39
# Schaefer	10	120	110
# Schaefer	50	120	120
# where "Trt"=treatment; "D"=dose or concentreation; "N"=Total number of lavae;	"R"=response number; "Fairfax","Pixley" and "Schaefer" are 3 treatments

# Set pi (significant levels)
PIs <- c(0.10, 0.50, 0.90, 0.99)

# Data standardizing
Trt.1 <- factor(Data$Trt)
ind.cor <- which(Data$D == 0)   # Identify natural response rows
ind.trt <- match(Trt.1[ind.cor], levels(Trt.1))
trt.cor <- levels(Trt.1)[ind.trt]
const <- Data$R[ind.cor] / Data$N[ind.cor]  #Calculating natural response
data.new <- Data[-ind.cor, ]

data.1 <- na.omit(data.new[data.new$D != 0, ])  # Delete rows in which dead number is 0

# Regularize treatment variable
Trt.1 <- as.character(data.1$Trt)
Trt.2 <- factor(Trt.1, ordered = FALSE, levels = trt.cor)

# Define variables
ld <- log10(data.1$D)
N <- data.1$N
R <- data.1$R

# Initialize  variables
ld.pi <- rep(NA, length(levels(Trt.2)))
SE <- rep(NA, length(levels(Trt.2)))
alpha <- rep(NA, length(levels(Trt.2)))
beta <- rep(NA, length(levels(Trt.2)))
va <- rep(NA, length(levels(Trt.2)))
vb <- rep(NA, length(levels(Trt.2)))
vab <- rep(NA, length(levels(Trt.2)))
bind <- data.frame()

#  Identify valid treatments
trt.valid <- na.omit(levels(Trt.2))
trt.valid <- trt.valid[trt.valid != ""]

# ---- Modeling loops ----
for (i in seq_along(trt.valid)) {
    cat("\n===== Processing Treatment:", trt.valid[i], "=====\n")
    ind <- which(Trt.2 == trt.valid[i])
  if (length(ind) == 0) {
    cat("Warning: No data for treatment", trt.valid[i], "\n")
    next
  }
  
  ld.i <- ld[ind]
  N.i <- N[ind]
  R.i <- R[ind]
  const.i <- ifelse(i <= length(const), const[i], 0)
  
  if (any(N.i == 0)) {
    cat("Warning: N == 0 detected in", trt.valid[i], "\n")
    next
  }

  #Correct with the natural response rate using Abbott's equation
  r.i <- (R.i / N.i - const.i) / (1 - const.i)
  
  if (all(is.na(r.i)) || all(!is.finite(r.i))) {
    cat("Warning: invalid probability in", trt.valid[i], "\n")
    next
  }
  
  data.2 <- data.frame(log.dose = ld.i, N = N.i, R = R.i, probability = r.i)
  data.3 <- subset(data.2, probability >= 0)
  
  if (nrow(data.3) == 0) {
    cat("Warning: no valid probability >= 0 in", trt.valid[i], "\n")
    next
  }
  
  # Preliminary model fitting
  mod.1 <- glm(probability ~ log.dose, data = data.3, family = binomial(link = "logit"), weights = N)
  
  data.4 <- subset(data.2, probability < 0)
  if (nrow(data.4) > 0) {
    pred.1 <- predict(mod.1, newdata = data.4)
    data.5 <- data.frame(log.dose = data.4$log.dose, N = data.4$N,
                         R = pnorm(pred.1) * data.4$N, probability = pnorm(pred.1))
  } else {
    data.5 <- data.frame()
  }
  
  data.6 <- rbind(data.5, data.3)
  
  #Final model fitting
  mod <- glm(probability ~ log.dose, data = data.6, family = binomial(link = "logit"), weights = N)
  
  # Defining parameters
  alpha[i] <- summary(mod)$coefficients[1, 1]
  beta[i] <- summary(mod)$coefficients[2, 1]
  va[i] <- vcov(mod)[1, 1]
  vb[i] <- vcov(mod)[2, 2]
  vab[i] <- vcov(mod)[1, 2]
  
  chi <- sum(residuals(mod, type = "pearson")^2)
  Pchi <- 1 - pchisq(chi, mod$df.residual)
  h <- chi / mod$df.residual
  t <- ifelse(h < 1, -qnorm(0.025), -qt(0.025, mod$df.residual))
  h.1 <- ifelse(h < 1, 1, h)
  g <- t^2 * vb[i] * h.1 / beta[i]^2
  
  lethal.dose.f <- NULL
  lethal.dose.d <- NULL
  
  for (pi in PIs) {
    ld.pi[i] <- dose.p(mod, p = pi)
    SE[i] <- attr(dose.p(mod, p = pi), "SE")
    
    # Calculate confidence interval of lethal doses: Fieller's theorem
    ld.cil.f <- ld.pi[i] + g/(1-g)*(ld.pi[i]+vab[i]/vb[i]) - t/(abs(beta[i])*(1-g)) * sqrt((va[i]+2*ld.pi[i]*vab[i]+ld.pi[i]^2*vb[i]-g*(va[i]-vab[i]^2/vb[i]))*h.1)
    ld.ciu.f <- ld.pi[i] + g/(1-g)*(ld.pi[i]+vab[i]/vb[i]) + t/(abs(beta[i])*(1-g)) * sqrt((va[i]+2*ld.pi[i]*vab[i]+ld.pi[i]^2*vb[i]-g*(va[i]-vab[i]^2/vb[i]))*h.1)
    
    # Calculate confidence interval of lethal doses: Delta method
    ld.cil.d <- ld.pi[i] + qnorm(0.025) * SE[i]
    ld.ciu.d <- ld.pi[i] - qnorm(0.025) * SE[i]
    
    # Defining lethal doses
    lethal.dose.f <- rbind(lethal.dose.f, c(pi, 10^ld.pi[i], 10^ld.cil.f, 10^ld.ciu.f))
    lethal.dose.d <- rbind(lethal.dose.d, c(pi, 10^ld.pi[i], 10^ld.cil.d, 10^ld.ciu.d))
  }
  
  # Layout lethal.doses
  cat("\nTreatment", i, ":", trt.valid[i], "\n")
  cat("LDpi and 95% CLs (logit model, Fieller's Theorem):\n")
  colnames(lethal.dose.f) <- c("pi", "LDpi", "95%CL:lower", "upper")
  print(round(lethal.dose.f, 3))
  
  cat("\nLDpi and 95% CLs (logit model, Delta Method):\n")
  colnames(lethal.dose.d) <- c("pi", "LDpi", "95%CL:lower", "upper")
  print(round(lethal.dose.d, 3))
  
  # Goodness-of-fits
  view <- data.frame(
    Dose = 10^ld.i, N = N.i, R = R.i,
    Expected = fitted(mod) * N.i,
    Residual = R.i - fitted(mod) * N.i,
    Probability = fitted(mod),
    Std.resid = residuals(mod, type = "pearson")
  )
  
  cat("\nCounts and Residuals:\n")
  print(round(view, 3))
  
  cat("\nChi-square test of goodness-of-fit:\n")
  cat("chi2 =", round(chi, 3), "; d.f. =", mod$df.residual, "; p =", round(Pchi, 3), "; h =", round(h, 3), "; g =", round(g, 3), "\n\n")
  
  # binding subset data
  data.7 <- cbind(Trt.2 = trt.valid[i], data.6)
  bind <- rbind(bind, data.7)
  
  # Draw dose-response curves for each treatment
  if (nrow(data.7) > 0 && all(data.7$probability >= 0) && all(data.7$probability <= 1)) {
    dev.new()
    p <- ggplot(data = data.7, aes(x = log.dose, y = probability)) +
      geom_point(size = 2) +
      geom_smooth(method = "glm", method.args = list(family = binomial(link = "logit")), aes(weight = N), colour = "red", se = TRUE) +
      labs(title = paste("Response - log(dose) curve:", trt.valid[i]), x = "log10(Dose)", y = "Probability") +
      theme_minimal()
    print(p)
  } else {
    cat("Warning: skipped plot for treatment", trt.valid[i], "\n")
  }
}
# ---- Calculate lethal.dose and potency ratios ----
for (pi in PIs) {
  potency.ratio <- NULL
  ratio.test <- NULL
  pairs <- NULL
  
  for (i in 1:length(trt.valid)) {
    for (j in 1:length(trt.valid)) {
      if (!is.na(trt.valid[i]) && !is.na(trt.valid[j]) && trt.valid[i] != trt.valid[j]) {
        
        # Subset Data
        mod.i <- glm(probability ~ log.dose, family = binomial(link = "logit"), data = bind, weights = N, subset = (Trt.2 == trt.valid[i]))
        mod.j <- glm(probability ~ log.dose, family = binomial(link = "logit"), data = bind, weights = N, subset = (Trt.2 == trt.valid[j]))
        
        # Calculate lethal.dose at pi
        ld_pi_i <- dose.p(mod.i, p = pi)
        SE_i <- attr(dose.p(mod.i, p = pi), "SE")
        
        ld_pi_j <- dose.p(mod.j, p = pi)
        SE_j <- attr(dose.p(mod.j, p = pi), "SE")
        
        #  Calculate lethal.dose ratios
        a <- ld_pi_i - ld_pi_j
        ratio <- 10^a
        sigma <- sqrt(SE_i^2 + SE_j^2)
        ratio.low <- 10^(a - 2 * sigma)
        ratio.up <- 10^(a + 2 * sigma)
        
        potency.ratio <- rbind(potency.ratio, c(ratio, ratio.low, ratio.up))
        
        # calculate potency ratios and significance test
        sigma2 <- sqrt(va[i]/alpha[i]^2 + vb[i]/beta[i]^2 + va[j]/alpha[j]^2 + vb[j]/beta[j]^2 - 2*vab[i]/(alpha[i]*beta[i]) - 2*vab[j]/(alpha[j]*beta[j]))
        z <- abs(a)/sigma2
        pztest <- 2 * pnorm(-z)
        pz <- p.adjust(pztest, method = "BH", n = length(trt.valid) * (length(trt.valid) - 1) / 2)
        
        ratio.test <- rbind(ratio.test, c(z, pz))
        
        m <- paste(trt.valid[i], "/", trt.valid[j])
        pairs <- append(pairs, m)
      }
    }
  }
  
  # Layout lethal.dose ratios  (Robertson et al., 2017):
  cat("LD", pi*100, "Ratio and 95% CLs (Robertson et al., 2017):\n")
  dimnames(potency.ratio) <- list(pairs, c("LD ratio", "95%CL: Lower", "Upper"))
  print(round(potency.ratio, 3))
  cat("\n")
  
  #  Layout equal lethal.dose test  (Wheeler et al., 2006):
  cat("LD", pi*100, "Ratio test (Wheeler et al., 2006):\n")
  dimnames(ratio.test) <- list(pairs, c("z", "p"))
  print(round(ratio.test, 5))
  cat("\n")
}
# ---- Parallelism Test ----
cat("\n===== Parallelism Tests Among All Lines =====\n")

# Regularize data
bind.clean <- bind[!is.na(bind$probability) & !is.na(bind$N) & bind$N > 0, ]
bind.clean$Trt.2 <- droplevels(factor(bind.clean$Trt.2))  # Delete invalid levels
bind.clean <- bind.clean[bind.clean$Trt.2 %in% names(which(table(bind.clean$Trt.2) > 0)), ]
bind.clean$Trt.2 <- droplevels(bind.clean$Trt.2)

# Fitting unconstrained model£¨with interaction£©
mod.unconst <- glm(probability ~ Trt.2 + Trt.2:log.dose, data = bind.clean, 
                   family = binomial(link = "logit"), weights = N)
df.unconst <- mod.unconst$df.residual

# Fitting constrained model£¨no interaction£©
mod.const <- glm(probability ~ Trt.2 + log.dose, data = bind.clean, 
                 family = binomial(link = "logit"), weights = N)
df.const <- mod.const$df.residual

# ---Three types of parallelism tests---
# chi-square test
chi.parel <- mod.const$deviance - mod.unconst$deviance
p.chi <- 1 - pchisq(chi.parel, df.const - df.unconst)

# log likelihood ratio tests
LR <- 2 * (logLik(mod.unconst)[1] - logLik(mod.const)[1])
df.LR1 <- attr(logLik(mod.unconst), "df")
df.LR2 <- attr(logLik(mod.const), "df")
p.LR <- 1 - pchisq(LR, df.LR1 - df.LR2)

# F test
F.stat <- (df.unconst / (df.const - df.unconst)) * 
  ((mod.const$deviance - mod.unconst$deviance) / mod.unconst$deviance)
pF <- 1 - pf(F.stat, df.const - df.unconst, df.unconst)

# ---Layout  parallelism test results ---
cat("Chi-square:", round(chi.parel, 3), "; df:", df.const - df.unconst, "; p:", round(p.chi, 4), "\n")
cat("Likelihood Ratio:", round(LR, 3), "; df:", df.LR1 - df.LR2, "; p:", round(p.LR, 4), "\n")
cat("F-test:", round(F.stat, 3), "; df:", df.const - df.unconst, ",", df.unconst, "; p:", round(pF, 4), "\n\n")

# ---- Visualize the relationship between the lines ----
cat("\n===== Correct Interaction Plot (Treatment ¡Á log(dose)) =====\n")

# Regularize data
bind.clean <- bind[!is.na(bind$probability) & !is.na(bind$N) & bind$N > 0, ]
bind.clean$Trt.2 <- droplevels(factor(bind.clean$Trt.2))  # Delete invalid levels
bind.clean <- bind.clean[bind.clean$Trt.2 %in% names(which(table(bind.clean$Trt.2) > 0)), ]
bind.clean$Trt.2 <- droplevels(bind.clean$Trt.2)

#  Fitting unconstrained model£¨with interaction£©
mod.unconst <- glm(probability ~ Trt.2 * log.dose, data = bind.clean,
                   family = binomial(link = "logit"), weights = N)

# Generate predictive data
newdata <- expand.grid(
  log.dose = seq(min(bind.clean$log.dose, na.rm=TRUE), max(bind.clean$log.dose, na.rm=TRUE), length.out=100),
  Trt.2 = levels(bind.clean$Trt.2)
)

newdata$fit <- predict(mod.unconst, newdata = newdata, type = "link")
newdata$se <- predict(mod.unconst, newdata = newdata, type = "link", se.fit = TRUE)$se.fit

# Generate 95%CIs
newdata$upper <- newdata$fit + 1.96 * newdata$se
newdata$lower <- newdata$fit - 1.96 * newdata$se

# ---Draw interaction effect diagram---
library(ggplot2)
dev.new()
ggplot(newdata, aes(x = log.dose, y = fit, color = Trt.2)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trt.2), alpha = 0.2, color=NA) +
  labs(
    title = "Interaction plot: Treatment x log(dose)",
    x = "log10(Dose)",
    y = "Logit(Response)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank()
  )
