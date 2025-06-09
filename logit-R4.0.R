# logit - log(dose) model（R 4.0 兼容版）
rm(list=ls())
library(MASS)
library(ggplot2)
library(effects)

# 读取数据
Data <- read.table("C:/Users/pxw/Desktop/JerseyCity.txt", header = TRUE, quote = "\"")

# 设置目标P值
PIs <- c(0.10, 0.50, 0.90, 0.99)

# 标准化处理分组
Trt.1 <- factor(Data$Trt)
ind.cor <- which(Data$D == 0)  # 识别自然死亡组
ind.trt <- match(Trt.1[ind.cor], levels(Trt.1))
trt.cor <- levels(Trt.1)[ind.trt]
const <- Data$R[ind.cor] / Data$N[ind.cor]  # 计算自然反应率
data.new <- Data[-ind.cor, ]

# 删除死亡数为0的记录
data.1 <- na.omit(data.new[data.new$D != 0, ])

# 分组变量
Trt.1 <- as.character(data.1$Trt)
Trt.2 <- factor(Trt.1, ordered = FALSE, levels = trt.cor)

# 处理数值
ld <- log10(data.1$D)
N <- data.1$N
R <- data.1$R

# 初始化结果变量
ld.pi <- rep(NA, length(levels(Trt.2)))
SE <- rep(NA, length(levels(Trt.2)))
alpha <- rep(NA, length(levels(Trt.2)))
beta <- rep(NA, length(levels(Trt.2)))
va <- rep(NA, length(levels(Trt.2)))
vb <- rep(NA, length(levels(Trt.2)))
vab <- rep(NA, length(levels(Trt.2)))
bind <- data.frame()

# 提前清理有效分组
trt.valid <- na.omit(levels(Trt.2))
trt.valid <- trt.valid[trt.valid != ""]

# ---- 开始建模循环 ----
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
  
  # 初步拟合
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
  
  # 最终拟合
  mod <- glm(probability ~ log.dose, data = data.6, family = binomial(link = "logit"), weights = N)
  
  # 保存参数
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
    
    # Fieller法计算置信区间
    ld.cil.f <- ld.pi[i] + g/(1-g)*(ld.pi[i]+vab[i]/vb[i]) - t/(abs(beta[i])*(1-g)) * sqrt((va[i]+2*ld.pi[i]*vab[i]+ld.pi[i]^2*vb[i]-g*(va[i]-vab[i]^2/vb[i]))*h.1)
    ld.ciu.f <- ld.pi[i] + g/(1-g)*(ld.pi[i]+vab[i]/vb[i]) + t/(abs(beta[i])*(1-g)) * sqrt((va[i]+2*ld.pi[i]*vab[i]+ld.pi[i]^2*vb[i]-g*(va[i]-vab[i]^2/vb[i]))*h.1)
    
    # Delta法计算置信区间
    ld.cil.d <- ld.pi[i] + qnorm(0.025) * SE[i]
    ld.ciu.d <- ld.pi[i] - qnorm(0.025) * SE[i]
    
    # 保存
    lethal.dose.f <- rbind(lethal.dose.f, c(pi, 10^ld.pi[i], 10^ld.cil.f, 10^ld.ciu.f))
    lethal.dose.d <- rbind(lethal.dose.d, c(pi, 10^ld.pi[i], 10^ld.cil.d, 10^ld.ciu.d))
  }
  
  # 打印结果
  cat("\nTreatment", i, ":", trt.valid[i], "\n")
  cat("LDπ and 95% CLs (logit model, Fieller's Theorem):\n")
  colnames(lethal.dose.f) <- c("π", "LDπ", "95%CL:lower", "95%CL:upper")
  print(round(lethal.dose.f, 3))
  
  cat("\nLDπ and 95% CLs (logit model, Delta Method):\n")
  colnames(lethal.dose.d) <- c("π", "LDπ", "95%CL:lower", "95%CL:upper")
  print(round(lethal.dose.d, 3))
  
  # 残差分析
  view <- data.frame(
    Dose = 10^ld.i, N = N.i, R = R.i,
    Expected = fitted(mod) * N.i,
    Residual = R.i - fitted(mod) * N.i,
    Probability = fitted(mod),
    Std.resid = residuals(mod, type = "pearson")
  )
  
  # 输出回归系数与标准误
  cat("\nCoefficients and Standard Errors:\n")
  print(summary(mod)$coefficients)
  
  # 输出方差-协方差矩阵
  cat("\nVariance-Covariance Matrix:\n")
  print(vcov(mod))

  cat("\nCounts and Residuals:\n")
  print(round(view, 3))
  
  cat("\nChi-square test of goodness-of-fit:\n")
  cat("χ2 =", round(chi, 3), "; d.f. =", mod$df.residual, "; p =", round(Pchi, 3), "; h =", round(h, 3), "; g =", round(g, 3), "\n\n")
  
  # 绑定数据
  data.7 <- cbind(Trt.2 = trt.valid[i], data.6)
  bind <- rbind(bind, data.7)
  
  # 绘制每个Treatment的dose-response曲线
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
# ---- 计算potency ratios和比值检验 ----
for (pi in PIs) {
  potency.ratio <- NULL
  ratio.test <- NULL
  pairs <- NULL
  
  for (i in 1:length(trt.valid)) {
    for (j in 1:length(trt.valid)) {
      if (!is.na(trt.valid[i]) && !is.na(trt.valid[j]) && trt.valid[i] != trt.valid[j]) {
        
        # 子集数据
        mod.i <- glm(probability ~ log.dose, family = binomial(link = "logit"), data = bind, weights = N, subset = (Trt.2 == trt.valid[i]))
        mod.j <- glm(probability ~ log.dose, family = binomial(link = "logit"), data = bind, weights = N, subset = (Trt.2 == trt.valid[j]))
        
        # 估计各组LDπ
        ld_pi_i <- dose.p(mod.i, p = pi)
        SE_i <- attr(dose.p(mod.i, p = pi), "SE")
        
        ld_pi_j <- dose.p(mod.j, p = pi)
        SE_j <- attr(dose.p(mod.j, p = pi), "SE")
        
        # 计算LD比值
        a <- ld_pi_i - ld_pi_j
        ratio <- 10^a
        sigma <- sqrt(SE_i^2 + SE_j^2)
        ratio.low <- 10^(a - 2 * sigma)
        ratio.up <- 10^(a + 2 * sigma)
        
        potency.ratio <- rbind(potency.ratio, c(ratio, ratio.low, ratio.up))
        
        # 计算比值检验
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
  
  # 输出LD比值结果
  cat("LD", pi*100, "% Ratio and 95% CLs (Robertson et al., 2017):\n")
  dimnames(potency.ratio) <- list(pairs, c("LD ratio", "95%CL Lower", "95%CL Upper"))
  print(round(potency.ratio, 3))
  cat("\n")
  
  # 输出比值检验结果
  cat("LD", pi*100, "% Ratio test (Wheeler et al., 2006):\n")
  dimnames(ratio.test) <- list(pairs, c("z", "p"))
  print(round(ratio.test, 5))
  cat("\n")
}
# ---- 平行性检验 (Parallelism Test) ----
cat("\n===== Parallelism Tests Among All Lines =====\n")

# ---【新加步骤：bind清理】---
# 去除无效数据
bind.clean <- bind[!is.na(bind$probability) & !is.na(bind$N) & bind$N > 0, ]
bind.clean$Trt.2 <- droplevels(factor(bind.clean$Trt.2))  # 删除无效水平

# 确保每个处理组都有数据
bind.clean <- bind.clean[bind.clean$Trt.2 %in% names(which(table(bind.clean$Trt.2) > 0)), ]
bind.clean$Trt.2 <- droplevels(bind.clean$Trt.2)

# ---【开始建模】---
# 拟合无约束模型（有交互）
mod.unconst <- glm(probability ~ Trt.2 + Trt.2:log.dose, data = bind.clean, 
                   family = binomial(link = "logit"), weights = N)
df.unconst <- mod.unconst$df.residual

# 拟合有约束模型（无交互）
mod.const <- glm(probability ~ Trt.2 + log.dose, data = bind.clean, 
                 family = binomial(link = "logit"), weights = N)
df.const <- mod.const$df.residual

# ---【三种平行性检验】---
# 卡方检验
chi.parel <- mod.const$deviance - mod.unconst$deviance
p.chi <- 1 - pchisq(chi.parel, df.const - df.unconst)

# 似然比检验
LR <- 2 * (logLik(mod.unconst)[1] - logLik(mod.const)[1])
df.LR1 <- attr(logLik(mod.unconst), "df")
df.LR2 <- attr(logLik(mod.const), "df")
p.LR <- 1 - pchisq(LR, df.LR1 - df.LR2)

# F检验
F.stat <- (df.unconst / (df.const - df.unconst)) * 
  ((mod.const$deviance - mod.unconst$deviance) / mod.unconst$deviance)
pF <- 1 - pf(F.stat, df.const - df.unconst, df.unconst)

# ---【输出结果】---
cat("Chi-square:", round(chi.parel, 3), "; df:", df.const - df.unconst, "; p:", round(p.chi, 4), "\n")
cat("Likelihood Ratio:", round(LR, 3), "; df:", df.LR1 - df.LR2, "; p:", round(p.LR, 4), "\n")
cat("F-test:", round(F.stat, 3), "; df:", df.const - df.unconst, ",", df.unconst, "; p:", round(pF, 4), "\n\n")

# ---- 正确绘制交互效应图（改写版） ----
cat("\n===== Correct Interaction Plot (Treatment × log(dose)) =====\n")

# 清理无效数据
bind.clean <- bind[!is.na(bind$probability) & !is.na(bind$N) & bind$N > 0, ]
bind.clean$Trt.2 <- droplevels(factor(bind.clean$Trt.2))  # 删除无效水平

# 确保每组都有数据
bind.clean <- bind.clean[bind.clean$Trt.2 %in% names(which(table(bind.clean$Trt.2) > 0)), ]
bind.clean$Trt.2 <- droplevels(bind.clean$Trt.2)

# 拟合无约束模型（有交互）
mod.unconst <- glm(probability ~ Trt.2 * log.dose, data = bind.clean,
                   family = binomial(link = "logit"), weights = N)

# 生成新的预测数据
newdata <- expand.grid(
  log.dose = seq(min(bind.clean$log.dose, na.rm=TRUE), max(bind.clean$log.dose, na.rm=TRUE), length.out=100),
  Trt.2 = levels(bind.clean$Trt.2)
)

# 进行预测（用link scale）
newdata$fit <- predict(mod.unconst, newdata = newdata, type = "link")
newdata$se <- predict(mod.unconst, newdata = newdata, type = "link", se.fit = TRUE)$se.fit

# 计算上下95%置信区间
newdata$upper <- newdata$fit + 1.96 * newdata$se
newdata$lower <- newdata$fit - 1.96 * newdata$se

# ---【绘制交互效应图】---
library(ggplot2)
dev.new()
ggplot(newdata, aes(x = log.dose, y = fit, color = Trt.2)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trt.2), alpha = 0.2, color=NA) +
  labs(
    title = "Interaction plot: Treatment × log(dose)",
    x = "log10(Dose)",
    y = "Logit(Response)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank()
  )
