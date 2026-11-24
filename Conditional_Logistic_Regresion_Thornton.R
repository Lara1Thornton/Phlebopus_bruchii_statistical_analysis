#Phlebopus_bruchii_statistical_analysis

rm(list = ls())
# Set the working directory
setwd("working directory") 
#Basic data matrix
data1 <-read.csv("NAME.csv")

#call packages
library(tidyverse)
library(survival)
library(broom)
library(GGally)
#First step: create pairs

data1 <- data %>%
  mutate(
    # Create the pair ID by removing "A" from the end
    par_id = str_remove(ID.punto, "A$"),
    # Check if it is an absence (ends in A)
    es_ausencia = str_ends(ID.punto, "A")
  )
complete_pairs <- data1 %>%
  group_by(pair_id) %>%
  summarise(
    n_sites = n(),
    presence = any(Presence.Absence == 1),
    absence = any(Presence.Absence == 0),
    is_pair_complete = (n_sites == 2) && presence && absence
  ) %>%
  filter(is_pair_complete)

#control number of sites
cat("Total number of sites in the original database:", nrow(data1), "\n")
cat("Paired complete sites", nrow(complete_pairs), "\n")

# Step 2: verify matching variables

# Paired Wilcoxon test for variables used in the matching
variables_matching <- c("Cob_canopia", "Cob_arbustos", "Pastos", "Latifoliadas") #use the name of your variables used for matching the paired sites

balance_results <- map_dfr(variables_matching, function(var) {
  # Extract paired values
  vals <- complete_pairs %>%
    select(par_id, Presence.Absence, all_of(var)) %>%
    pivot_wider(
      id_cols = par_id,
      names_from = Presence.Absence,
      values_from = all_of(var),
      names_prefix = "state_"
    ) %>%
    rename(presence = state_1, absence = state_0)
 # Paired test
 test <- wilcox.test(vals$presence, vals$absence, paired = TRUE)
 # Descriptive statistics
 summary <- complete_pairs %>%
    group_by(Presence.Absence) %>%
    summarise(
      media = mean(!!sym(var), na.rm = TRUE),
      sd = sd(!!sym(var), na.rm = TRUE),
      .groups = 'drop'
    )
  tibble(
    variable = var,
    mean_presence = summary$media[summary$Presence.Absence == 1],
    mean_absence = summary$media[summary$Presence.Absence == 0],
    p_value = test$p.value,
    balanced = test$p.value > 0.05
  )
})

print(balance_results)

# Step 3: Check colinearity 
model_data <- complete_pairs %>%
  select(
    par_id_num,
    Presence.Absence,
    humidity, pH, CE, SOC, C.N, Nitrogen, Phosphorus, #choose environmental predictors
    Ca, Mg, K, Na, CIC
  ) %>%
  rename(presence = Presence.Absence)

## Control NAs
print(colSums(is.na(model_data)))

cor_test <- model_data %>%
  select(V1, V2) %>%                    #V1=environmetal variable 1 #V2= environmetal variable 2
  cor.test(use = "complete.obs")
cor_test <-cor.test(model_data$V1, model_data$V2, method = "spearman")

var_exp <- subset(model_data, select= c("all the predictor you want to check" ))

ggpairs(var_exp)

#If r > 0.7, they are redundant and explain the same effect.
#remove collinear variables

# Step 4: Conditional Logistic Regression
model <- clogit(
  presence ~   V1 + V2 + V3 + V4 + Vz + #V1,V2..VZ = non collinear predictors
    strata(par_id_num),
  data = model_data,
  method = "exact"
)

summary(model)

# Step 5: Check model assumptions

##check the linearity assumption
res_mart <- residuals(model, type = "martingale")
### Check for linearity for each predictor

plot(model_data$V1, res_mart + coef(model)["V1"] * model_data$V1,
     xlab = "V1", ylab = "Residual component + adjustment")
lines(lowess(model_data$V1, res_mart + coef(model)["V1"] * model_data$V1), col = "red")

## check outliers

#### Extract the deviation residuals
res_dev <- residuals(model, type = "deviance")

### Plot the deviance residuals against the observation index
plot(res_dev, 
     ylab = "deviation residuals", 
     xlab = "Subject Index",
     main = "Deviance Residuals Plot",
     col = "blue", pch = 16)
abline(h = c(-2, 2), col = "red", lty = 2) # Standard warning lines

## DFBETAS: influence of each observation on each coefficient
dfbetas <- residuals(model, type = "dfbeta")
dim(dfbetas)
# 3. warning threshold: |dfbeta| > 2/sqrt(n)
# Para n=56: 2/sqrt(56) ≈ 0.267 (change)
umbral <- 2 / sqrt(nrow(dfbetas))


#Table for manuscript
library(knitr)
library(kableExtra)

table_results <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    # format numbers
    beta = round(estimate, 3),
    SE = round(std.error, 3),
    OR = round(estimate, 2),
    IC_inf = round(conf.low, 2),
    IC_sup = round(conf.high, 2),
    z = round(statistic, 3),
    p = round(p.value, 3),
    # Create a formatted Confidence Interval column
    CI_95 = paste0(IC_inf, "–", IC_sup),
    # Marcar significativos
    sig = ifelse(p.value < 0.05, "**", "")
  ) %>%
  select(term, beta, SE, OR, IC_95, z, p, sig) %>%
  # Renombrar variables para presentación
  mutate(
    Variable = case_when(
      term == "V1" ~ "name V1",
      term == "pH" ~ "pH",
      term == "Carbon" ~ "SOC",
      term == "Nitrogen" ~ "N",
      term == "V2" ~ "name V2",
      term == "CIC" ~ "CIC",
          TRUE ~ term
    )
  ) %>%
  select(Variable, beta, SE, OR, IC_95, z, p, sig)

## Save as CSV for Excel/Word
write.csv(tabla_final, "table_result.csv", row.names = FALSE)

