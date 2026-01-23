# Load necessary libraries
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(performance)
library(ggplot2)

# Load the dataset
data <- read.csv("C:/Users/Beth/Documents/Proc Roy Soc Paper/ProcRoySoc/for_R/filtered_df_with_colony_and_ITS.csv")

summary(data)
data$BeeID <- as.factor(data$BeeID)
data$event_type <- as.factor(data$event_type)
data$ColonyID <- as.factor(data$ColonyID)
summary(data)
colnames(data)

# MSMR distribution
ggplot(data, aes(x = mass_specific_mL_g_h)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() 

ggplot(data, aes(x = weight_g, y = log(mass_specific_mL_g_h), color = event_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ ColonyID) +
  theme_minimal()

ggplot(data, aes(x = ITS_mm, y = log(mass_specific_mL_g_h), color = event_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ ColonyID) +
  theme_minimal()

ggplot(data, aes(x = event_type, y = log(mass_specific_mL_g_h))) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_jitter(width = 0.1, alpha = 0.3) +
  facet_wrap(~ ColonyID)

# Models
model_log <- lmer(log10(mass_specific_mL_g_h) ~ event_type + ColonyID + ( 1|BeeID), data = data)

#FINAL MODEL
model_log2 <- lmer(log10(mass_specific_mL_g_h) ~ event_type + ColonyID  + (1+ event_type | BeeID), data = data)

#Check if variance slope improves model fit= sig
anova(model_log, model_log2, refit = FALSE)
isSingular(model_log2)

simulateResiduals(model_log2, plot = TRUE)
plotResiduals(model_log2, data$event_type)
plotResiduals(model_log2, data$ColonyID)
testDispersion(model_log2)
testUniformity(model_log2)

summary(model_log2)
anova(model_log2, type = "II")

emmeans(model_log2, pairwise ~ event_type, type = "response")

emm <- emmeans(model_log2, ~ event_type)
summary(emm)
contrast(emm, method = "pairwise", adjust = "none")  # No multiple testing correction needed here

emm_df <- as.data.frame(summary(emm, type = "response"))
emm_df

p_msmr <- ggplot(emm_df, aes(x = event_type, y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(x = "Event Type", y = "Mass-Specific Metabolic Rate (mL/g/h)") +
  theme_minimal(base_size = 14)
p_msmr

#Calculate ICC for model_log2

m2 <- model_log2  # your fitted model

vc <- VarCorr(m2)
G  <- as.matrix(vc$BeeID)          # random-effect (co)variance for BeeID: intercept/slope(s)
var_res <- attr(vc, "sc")^2        # residual variance

Xlev <- model.matrix(~ event_type, data = data)  # fixed-effect design for event_type (same coding used in model)
Zlev <- Xlev                                     # same columns correspond to random intercept + event_type slopes

# unique rows = one per event_type level (under current contrasts)
Zuniq <- unique(Zlev)

# between-bee variance for each event_type level:
var_bee_by_level <- apply(Zuniq, 1, function(z) as.numeric(t(z) %*% G %*% z))

ICC_by_level <- var_bee_by_level / (var_bee_by_level + var_res)

out <- data.frame(
  event_type = levels(data$event_type),
  var_bee = var_bee_by_level[seq_along(levels(data$event_type))],
  var_res = var_res,
  ICC = ICC_by_level[seq_along(levels(data$event_type))]
)

out

#Extract 95% CI

set.seed(123)

boot_icc <- bootMer(
  model_log2,
  FUN = function(fit) {
    vc <- VarCorr(fit)
    G  <- as.matrix(vc$BeeID)
    var_res <- attr(vc, "sc")^2
    
    # model matrix for event_type
    X <- model.matrix(~ event_type, data = data)
    Zuniq <- unique(X)
    
    icc <- apply(Zuniq, 1, function(z) {
      vb <- as.numeric(t(z) %*% G %*% z)
      vb / (vb + var_res)
    })
    
    names(icc) <- levels(data$event_type)
    icc
  },
  nsim = 1000,
  type = "parametric"
)

icc_ci <- t(apply(
  boot_icc$t,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
))

icc_df <- data.frame(
  event_type = rownames(icc_ci),
  ICC = icc_ci[, "50%"],
  lower = icc_ci[, "2.5%"],
  upper = icc_ci[, "97.5%"]
)

icc_df

#Plot ICC and CI

ggplot(icc_df, aes(x = event_type, y = ICC)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  ylim(0, 1) +
  labs(
    x = "Event type",
    y = "Individual-level repeatability (ICC)",
    title = "Repeatability of metabolic rate by behavior"
  ) +
  theme_classic(base_size = 14)

#Test if ICC sig different between behaviours = n.s. 
icc_diff <- boot_icc$t[, "buzz"] - boot_icc$t[, "flight"]
quantile(icc_diff, probs = c(0.025, 0.5, 0.975))
