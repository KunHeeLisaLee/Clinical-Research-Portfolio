# ------------------------------------------------------------------------------
# This is a code sample from my previous work as a Clinical Data Analyst.
# It supports a published research study on the risk of ophthalmologic disease 
# (uveitis) in patients with dermatologic disease (psoriasis) in the Korean population.
#
# The original script included over 750 lines of code. Here, I’ve reorganized it 
# to focus on the key components that demonstrate my applied data skills:
#   - Study/control group preparation
#   - Key data cleaning and comorbidity tagging
#   - Propensity score matching
#   - Table 1: Baseline demographics
#   - Table 2: Uveitis incidence using Poisson regression
#   - Figure 2: Cumulative incidence curves over time
#
# This work was published in a high-impact dermatology journal.
# For holistic overview, please refer to the paper "Publication1_Psoriasis_Uveitis.pdf"
# uploaded in my `Clinical Research Portfolio` GitHub repository. 
# ------------------------------------------------------------------------------


library(data.table);library(magrittr);library(parallel);library(MatchIt)
library(tableone);library(survival);library(survminer);library(ggplot2)

# Load and filter data
t20 <- fread("input.csv")
setDT(t20)
t20[, usage_date := as.Date(as.character(RECU_FR_DD), format = "%Y%m%d")]

# Exclude subjects with prior diagnosis (e.g., psoriasis or uveitis)
excluded_ids <- t20[RECU_FR_DD < 20110101 & grepl("ICD-10 for uveitis|psoriasis", MAIN_SICK), unique(JID)]
t20 <- t20[!JID %in% excluded_ids]

# Define study and control groups (psoriasis vs urticaria)
study_group <- t20[grepl("ICD-10 for psoriasis", MAIN_SICK), .SD[1], by = JID]
control_group <- t20[grepl("ICD-10 for urticaria", MAIN_SICK) & !(JID %in% study_group$JID), .SD[1], by = JID]

# Combine and label groups
study_group[, Group := "Psoriasis"]
control_group[, Group := "Control"]
combined <- rbindlist(list(study_group, control_group), fill = TRUE)

# Add patient-level attributes
combined[, `:=`(Sex = SEX_TP_CD, Age = PAT_AGE, IndexDate = as.Date(as.character(RECU_FR_DD), "%Y%m%d"))]

# Add comorbidity flags (e.g., DM, HTN) using diagnosis history
comorb <- t20[grepl("ICD-10 for diabetes", MAIN_SICK), .(Prev_DM = 1), by = JID]
combined <- merge(combined, comorb, by = "JID", all.x = TRUE)
combined[is.na(Prev_DM), Prev_DM := 0]

# Propensity score matching (1:2) using demographics and comorbidities
combined[, Treatment := ifelse(Group == "Psoriasis", 1, 0)]
match_data <- matchit(Treatment ~ Age + Sex + Prev_DM, data = combined, ratio = 2)
matched_ids <- match.data(match_data)$JID
matched_data <- combined[JID %in% matched_ids]

# Outcome: First uveitis event after index date
uveitis_events <- t20[grepl("ICD-10 for uveitis", MAIN_SICK), .(dx_date = min(RECU_FR_DD)), by = JID]
matched_data <- merge(matched_data, uveitis_events, by = "JID", all.x = TRUE)
matched_data[, dx_date := as.Date(as.character(dx_date), format = "%Y%m%d")]
matched_data[, FU_year := as.numeric(as.Date("2020-12-31") - IndexDate) / 365.25]
matched_data[, time_to_event := as.numeric(dx_date - IndexDate) / 365.25]
matched_data[, uveitis_flag := as.integer(!is.na(dx_date) & dx_date > IndexDate)]

# ---Table 1: Baseline Demographics---
vars <- c("Sex", "Age", "Prev_DM")
tb1 <- CreateTableOne(vars = vars, strata = "Group", data = matched_data)
write.csv(print(tb1), "dataout/table1_baseline.csv")

# ---Table 2: Poisson Regression for Uveitis Incidence---
model <- glm(uveitis_flag ~ offset(log(FU_year)) + Group, data = matched_data, family = "poisson")
summary(model)

# Calculate incidence per 1000 PY and RR
inci <- matched_data[, .(
  N = .N,
  Events = sum(uveitis_flag),
  PY = sum(FU_year),
  Incidence_per_1000PY = round(sum(uveitis_flag)/sum(FU_year) * 1000, 2)
), by = Group]

rr <- exp(coef(model)[2])
ci <- exp(confint(model)[2, ])
inci[, `:=`(RR = round(rr, 2), CI = paste0("(", round(ci[1], 2), "~", round(ci[2], 2), ")"))]
fwrite(inci, "dataout/table2_inci_rr.csv")

# ---Figure 2: Cumulative Incidence Curve---
matched_data[, time_to_event := ifelse(is.na(time_to_event), FU_year, time_to_event)]
surv_obj <- survfit(Surv(time_to_event, uveitis_flag) ~ Group, data = matched_data)
ggsurvplot(surv_obj, data = matched_data, pval = TRUE, risk.table = TRUE,
           xlab = "Years", ylab = "Cumulative Incidence of Uveitis", 
           legend.title = "Group", fun = "event", break.x.by = 1)
ggsave("dataout/figure2_cumulative_incidence.png", width = 9, height = 6)

