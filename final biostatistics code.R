################################################
# BIOSTATISTICS PROJECT — SURVIVAL ANALYSIS
# Dataset: Mayo Clinic Primary Biliary Cirrhosis
################################################

set.seed(123)

# ── LIBRARIES ──────────────────────────────────
library(survival)
library(survminer)
library(ggplot2)

# ── 1. IMPORT DATA ─────────────────────────────
data <- read.csv("cirrhosis.csv")

# ── 2. DATA CLEANING ───────────────────────────
# تحويل Status
data$Status <- ifelse(data$Status == "D", "Dead", "Censored")

# تحويل الأيام لسنين
data$N_Days <- data$N_Days / 365.25
names(data)[names(data) == "N_Days"] <- "Years"

# تحويل المتغيرات لـ factors
data$Status      <- factor(data$Status,
                           levels = c("Censored", "Dead"))

data$Sex         <- factor(data$Sex,
                           labels = c("Female", "Male"))

data$Drug        <- factor(data$Drug,
                           labels = c("D-penicillamine", "Placebo"))

data$Ascites     <- factor(data$Ascites,
                           labels = c("No", "Yes"))

data$Hepatomegaly <- factor(data$Hepatomegaly,
                            levels = c("N", "Y"),
                            labels = c("No", "Yes"))

data$Spiders     <- factor(data$Spiders,
                           levels = c("N", "Y"),
                           labels = c("No", "Yes"))

data$Edema       <- factor(data$Edema,
                           labels = c("No Edema","Slight Edema","Severe Edema"))

data$Stage       <- factor(data$Stage,
                           levels = c(1, 2, 3, 4),
                           labels = c("Stage 1","Stage 2","Stage 3","Stage 4"))

# ── 3. MISSING VALUES ──────────────────────────
cat("Total missing values in full dataset:", sum(is.na(data)), "\n")

final_data <- na.omit(data[, c("Years","Status","Age","Sex",
                               "Bilirubin","Albumin")])
cat("Rows after removing missing (key vars):", nrow(final_data), "\n")

# ── 4. DESCRIPTIVE ANALYSIS ────────────────────

# 4a. Numerical summary
num_vars <- data[, c("Years","Age","Bilirubin","Cholesterol","Albumin",
                     "Copper","Alk_Phos","SGOT","Tryglicerides",
                     "Platelets","Prothrombin")]
summary(num_vars)
sapply(num_vars, sd, na.rm = TRUE)

# 4b. Categorical frequency tables (counts + percentages)
cat_vars <- c("Status","Drug","Sex","Ascites",
              "Hepatomegaly","Spiders","Edema","Stage")

for (var in cat_vars) {
  cat("\n===", var, "===\n")
  tbl <- table(data[[var]], useNA = "ifany")
  pct <- round(prop.table(tbl) * 100, 1)
  print(cbind(Count = tbl, Percent = pct))
}

# 4c. Bar charts for categorical variables
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
for (var in cat_vars) {
  counts <- table(data[[var]], useNA = "ifany")
  barplot(counts,
          main = paste("Distribution of", var),
          col  = rainbow(length(counts)),
          ylab = "Number of Patients",
          las  = 2)
}
par(mfrow = c(1, 1))

# ── 5. EVENT CHART ─────────────────────────────
ggplot(data[1:50, ], aes(x = Years, y = as.factor(ID))) +
  geom_segment(aes(x = 0, xend = Years,
                   y = as.factor(ID), yend = as.factor(ID)),
               color = "lightgrey") +
  geom_point(aes(shape = Status, color = Status), size = 2) +
  scale_shape_manual(values = c("Censored" = 16, "Dead" = 17)) +
  scale_color_manual(values = c("Censored" = "blue", "Dead" = "red")) +
  labs(x = "Time (Years)", y = "Patient ID",
       title = "Patient Follow-up and Outcomes (First 50 Patients)") +
  theme_minimal()

# ── 6. KAPLAN-MEIER ────────────────────────────

# Overall KM
km_fit <- survfit(Surv(Years, Status == "Dead") ~ 1, data = data)

# Median survival + survival at key timepoints
print(km_fit)
summary(km_fit, times = c(1, 3, 5, 10))

ggsurvplot(km_fit,
           data = data,
           title = "Overall Kaplan-Meier Survival Estimate",
           xlab  = "Survival Time (Years)",
           surv.median.line = "hv",
           caption = "Cirrhosis data analysis")

# KM by Ascites
km_ascites <- survfit(Surv(Years, Status == "Dead") ~ Ascites,
                      conf.type = "log-log", data = data)
ggsurvplot(km_ascites, data = data, pval = TRUE,
           title    = "KM Survival by Ascites Status",
           xlab     = "Survival Time (Years)",
           surv.median.line = "hv",
           caption  = "Cirrhosis data analysis")
# KM by sex
km_sex <- survfit(Surv(Years, Status == "Dead") ~ Sex,
                  conf.type = "log-log", data = data)
ggsurvplot(km_sex, data = data, pval = TRUE,
           title   = "KM Survival by Sex",
           xlab    = "Survival Time (Years)",
           surv.median.line = "hv",
           caption = "Cirrhosis data analysis")

# KM by Hepatomegaly
km_hepato <- survfit(Surv(Years, Status == "Dead") ~ Hepatomegaly,
                     conf.type = "log-log", data = data)
ggsurvplot(km_hepato, data = data, pval = TRUE,
           title   = "KM Survival by Hepatomegaly",
           xlab    = "Survival Time (Years)",
           surv.median.line = "hv",
           legend.labs = c("No Hepatomegaly", "With Hepatomegaly"),
           palette = c("#00AFBB", "#FC4E07"),
           caption = "Cirrhosis data analysis")

# ── 7. COX REGRESSION ──────────────────────────

# Full model
result.cox_full <- coxph(
  Surv(Years, Status == "Dead") ~ Age + Sex + Drug + Ascites +
    Hepatomegaly + Spiders + Edema +
    Bilirubin + Cholesterol + Albumin +
    Copper + Alk_Phos + SGOT +
    Tryglicerides + Platelets + Prothrombin + Stage,
  data = data
)
summary(result.cox_full)

# Reduced model (significant predictors only)
result.cox_reduced <- coxph(
  Surv(Years, Status == "Dead") ~ Age + Edema + Bilirubin +
    Albumin + Copper + Prothrombin + Stage,
  data = data
)
summary(result.cox_reduced)

# Forest plot
ggforest(result.cox_reduced, data = data)

# ── 8. PH ASSUMPTION TEST ──────────────────────

# Test on reduced model
ph_test <- cox.zph(result.cox_reduced)
print(ph_test)

# Schoenfeld residual plots
png("schoenfeld_plots.png", width = 3500, height = 2500, res = 300)

ggcoxzph(ph_test,
         point.col   = "steelblue",
         point.size  = 1,
         point.alpha = 0.5,
         caption     = "Schoenfeld Residuals — Reduced Model",
         ggtheme     = theme_bw(base_size = 14),   # bigger cleaner font
         font.main   = 13)

dev.off()
shell.exec("schoenfeld_plots.png")
# ── 9. FIX PH VIOLATION — STRATIFICATION ───────

# Create groups for violating variables
data$Bilirubin_group <- cut(data$Bilirubin,
                            breaks = c(0, 1, 3, Inf),
                            labels = c("Normal","Elevated","High"),
                            right  = TRUE)

data$Prothrombin_group <- cut(data$Prothrombin,
                              breaks = c(0, 10.6, Inf),
                              labels = c("Normal","Prolonged"),
                              right  = TRUE)

# Stratified Cox model
cox_stratified <- coxph(
  Surv(Years, Status == "Dead") ~
    Age + Edema + Albumin + Copper +
    strata(Bilirubin_group) +
    strata(Prothrombin_group),
  data = data
)
summary(cox_stratified)

# Verify PH now satisfied
ph_test_stratified <- cox.zph(cox_stratified)
print(ph_test_stratified)
# Global p should be > 0.05 ✓

# Schoenfeld plots for stratified model
png("schoenfeld_stratified.png", width = 4000, height = 3000, res = 300)
ggcoxzph(ph_test_stratified,
         point.col   = "steelblue",
         point.size  = 1,
         point.alpha = 0.4,
         caption     = "Schoenfeld Residuals — Stratified Model (PH Fixed)",
         ggtheme     = theme_bw(base_size = 13) +
           theme(
             axis.title.y = element_text(size = 9),   # shrink y label
             axis.title.x = element_text(size = 10),
             plot.title   = element_text(size = 11),
             strip.text   = element_text(size = 10)
           ))
dev.off()
# Open it directly
shell.exec("schoenfeld_stratified.png")

# Stratified survival curves
km_strat <- survfit(cox_stratified)
ggsurvplot(km_strat,
           data    = data,
           title   = "Stratified Cox Model Survival Curves",
           xlab    = "Time (Years)",
           caption = "Cirrhosis data analysis")

################################################
# END OF SCRIPT
################################################