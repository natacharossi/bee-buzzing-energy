# Load necessary libraries
library(lme4)
library(lmerTest)
library(DHARMa)
library(emmeans)
library(performance)
library(ggplot2)

# Load the dataset
data <- read.csv("C:/Users/Beth/Documents/Proc Roy Soc Paper/ProcRoySoc/for_R/filtered_df_with_colony_and_ITS.csv")

summary(data)
data$BeeID <- as.factor(data$BeeID)
data$event_type <- as.factor(data$event_type)
data$ColonyID <- as.factor(data$ColonyID)
summary(data)

# MR distribution
ggplot(data, aes(x = metabolic_rate_mL_h)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() 

#Model including body mass, body mass * behaviour interaction and slope variance
lmm1 <- lmer(log10(metabolic_rate_mL_h) ~ event_type + ColonyID + 
               (log10(weight_g)) + event_type*(log(weight_g)) + (1 +event_type| BeeID), data = data)
simulateResiduals(lmm1, plot = TRUE)
summary(lmm1)
anova(lmm1, type = "II")

#FINAL MODEL: drop behaviour x body mass interaction. 
lmm2 <- lmer(log10(metabolic_rate_mL_h) ~ 
               event_type +
               ColonyID +
               log10(weight_g) +
               (1 + event_type | BeeID),
             data = data)

simulateResiduals(lmm2, plot = TRUE)
plotResiduals(lmm2, data$event_type)
plotResiduals(lmm2, data$ColonyID)
plotResiduals(lmm2, data$log(weight_g))
testDispersion(lmm2)
testUniformity(lmm2)

#Check if removing the interaction had a significant effect on model fit= n.s.
anova(lmm1, lmm2, refit = FALSE)
isSingular(lmm2)
AIC(lmm1, lmm2)

#Model summary
summary(lmm2)
anova(lmm2, type = "II")

emmeans(lmm2, pairwise ~ event_type, type = "response")

# Estimated marginal means on the log scale
em <- emmeans(lmm2, ~ event_type)
em

# Back-transform to original scale
em_trans <- summary(em, type = "response")  # uses exp() under the hood if response was log-transformed
em_trans
ref_grid(lmm2)
#weight_g = 0.203g

#Calculate ICC for lmm2
library(lme4)

m2 <- lmm2  # your fitted model

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
  lmm2,
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

#Check if random effect and
#random slope had a significant effect
lmm3 <- lmer(log(metabolic_rate_mL_h) ~ 
               event_type +
               ColonyID +
               log(weight_g) +
               (1 | BeeID),
             data = data)

#Including random slope sig. improves fit so keep lmm2 as final model
anova(lmm2, lmm3, refit = FALSE)
isSingular(lmm3)


#Repeat analysis for ITS

cor(data$ITS_mm, data$weight_g, use = "complete.obs")
#relatively low correlation 0.3 compared to other studies

#Include both mass and ITS
lmm_size <- lmer(log(metabolic_rate_mL_h) ~ 
                   event_type + ColonyID + 
                   scale(log(weight_g)) + 
                   scale(log(ITS_mm)) + 
                   (1 + event_type | BeeID),
                 data = data)

#Model summary
summary(lmm_size)
anova(lmm_size, type = "II")

#Model with just ITS (but also interaction with behaviour)
lmm_its_interaction <- lmer(
  log10(metabolic_rate_mL_h) ~ event_type + ColonyID + log(ITS_mm) +
    event_type:log(ITS_mm) +
    (1 + event_type | BeeID),
  data = data
)

#Model summary- interaction NS
summary(lmm_its_interaction)
anova(lmm_its_interaction, type = "II")

#Model ITS without interaction

lmm_its <- lmer(
  log10(metabolic_rate_mL_h) ~ event_type + ColonyID + log10(ITS_mm) +
    (1 + event_type | BeeID),
  data = data
)

#Model summary
summary(lmm_its)
anova(lmm_its, type = "II")
#Compare models with and without  = n.s. effect
anova(lmm_its_interaction, lmm_its)
AIC (lmm_its_interaction, lmm_its)

#Model with just mass (lmm2)
lmm_mass <- lmer(log(metabolic_rate_mL_h) ~ 
                  event_type + ColonyID + 
                   (log(weight_g)) + 
                  (1 + event_type | BeeID),
                data = data)

#Model summary
summary(lmm_mass)
anova(lmm_mass, type = "II")

AIC(lmm_mass, lmm_its, lmm_size)
#lmm_mass is best fitting model with lowest AIC 
#and the other models are both >8 AIC units worse. 

#Plot MR data with every event
ggplot(data,
       aes(x = log10(weight_g),
           y = log10(metabolic_rate_mL_h))) +
  
  # connect events from the same bee
  geom_line(aes(group = BeeID),
            alpha = 0.25,
            linewidth = 0.4,
            colour = "grey40") +
  
  # points: same colours as lines
  geom_point(aes(shape = event_type,
                 colour = event_type,
                 fill = event_type),
             alpha = 0.8,
             size = 2.5,
             stroke = 0.6) +
  
  # behaviour-specific regression lines
  geom_smooth(aes(colour = event_type, fill = event_type),
              method = "lm",
              se = TRUE,
              alpha = 0.20,
              linewidth = 1) +
  
  theme_classic(base_size = 14) +
  labs(x = "log10(Body mass, g)",
       y = "log10 Metabolic rate (mL/h)") +
  
  # filled circle + triangle
  scale_shape_manual(values = c(21, 24))

#Plot ITS data with every event
#Plot data with every event
ggplot(data,
       aes(x = log10(ITS_mm),
           y = log10(metabolic_rate_mL_h))) +
  
  # connect events from the same bee
  geom_line(aes(group = BeeID),
            alpha = 0.25,
            size = 0.4,
            colour = "grey40") +
  
  # points: same colours as lines
  geom_point(aes(shape = event_type,
                 colour = event_type,
                 fill = event_type),
             alpha = 0.8,
             size = 2.5,
             stroke = 0.6) +
  
  # behaviour-specific regression lines
  geom_smooth(aes(colour = event_type, fill = event_type),
              method = "lm",
              se = TRUE,
              alpha = 0.20,
              size = 1) +
  
  theme_classic(base_size = 14) +
  labs(x = "log10(ITS mm)",
       y = "log10 Metabolic rate (mL/h)") +
  
  # filled circle + triangle
  scale_shape_manual(values = c(21, 24))

#Run SMA analysis to calculate scaling exponent

data_sma <- data %>%
  mutate(
    logMR   = log10(metabolic_rate_mL_h),
    logMass = log10(weight_g),
    event_type = factor(event_type)  # ensure factor
  ) %>%
  # remove any rows with missing values in the key variables
  filter(
    !is.na(logMR),
    !is.na(logMass),
    !is.na(event_type)
  )

sma_all <- sma(logMR ~ logMass, data = data_sma)
summary(sma_all)

sma_buzz_evt <- sma(logMR ~ logMass,
                    data = data_sma,
                    subset = event_type == "buzz")

sma_flight_evt <- sma(logMR ~ logMass,
                      data = data_sma,
                      subset = event_type == "flight")

summary(sma_buzz_evt)
summary(sma_flight_evt)

sma_evt_behav <- sma(
  logMR ~ logMass * event_type,
  data = data_sma
)

summary(sma_evt_behav)


# Common slope, different elevation model 
sma_evt_behav_common <- sma(
  logMR ~ logMass + event_type,  # no interaction term
  data = data_sma
)

summary(sma_evt_behav_common)


coef_table <- summary(sma_evt_behav)$coef
coef_table


slopes <- coef_table[, "slope"]
lower  <- coef_table[, "lower limit"]
upper  <- coef_table[, "upper limit"]

slopes
lower
upper
