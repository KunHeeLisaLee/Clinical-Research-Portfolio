# ------------------------------------------------------------------------------
# This is a code sample from my previous work as a Clinical Data Analyst.
# It supports a population-based cohort study on the risk of CCVD in Korean patients
# with psoriasis, based on systemic antipsoriatic therapy exposure. 
# The original script included over 1,000 lines of code. Here, I’ve reorganized it 
# to focus key parts of the analysis, including:

# - Data cleaning and treatment group assignment
# - Table 1: Baseline demographics 
# - Table 2: Impact of systemic therapy duration on CCVD risk
# - Table 3: Subgroup analysis of CCVD risk by treatment type
# - Figure 2: GAM model showing adjusted odds ratios by percent treatment period (PTP)
#
# This work was published in a high-impact dermatology journal.
# For holistic overview, please refer to the paper "Publication2_AntipsoriaticTherapy_CCVD.pdf"
# uploaded in my `Clinical Research Portfolio` GitHub repository. 
# ------------------------------------------------------------------------------


library(data.table);library(magrittr);library(tableone);library(mgcv)

# Load cohort
nested <- fread("input.csv")
nested[, CCVD := factor(CCVD, levels = c(0,1), labels = c("No", "Yes"))]

# ---Table 1: Baseline Demographics---
vars.tb1 <- c("Sex", "Age", "Prev_DM", "Prev_HTN", "Prev_Dyslipidemia", "Percent_both")
tb1 <- CreateTableOne(vars = vars.tb1, strata = "CCVD", data = nested, 
                      factorVars = setdiff(vars.tb1, "Percent_both"))
print(tb1, nonnormal = "Percent_both", smd = TRUE)

# ---Table 2: Impact of systemic therapy duration on CCVD risk---
nested[, Percent_both_10 := Percent_both / 10]
model.logit <- glm(CCVD ~ Percent_both_10 + Age + Sex + Prev_DM + Prev_HTN + Prev_Dyslipidemia,
                   data = nested, family = "binomial")
summary(model.logit)

# ---Table 3: Subgroup Analysis by Treatment Type---
nested[, `:=`(DrugGroup = fifelse(biologics != "No", "Biologics",
                                  fifelse(drug != "No", "Oral", "No treatment")))]
nested[, DrugGroup := factor(DrugGroup, levels = c("No treatment", "Oral", "Biologics"))]

res.subgroup <- glm(CCVD ~ DrugGroup + Age + Sex + Prev_DM + Prev_HTN + Prev_Dyslipidemia,
                    data = nested[DrugGroup != "No treatment"], family = "binomial")
summary(res.subgroup)

# ---Figure 2: GAM Model---
model.gam <- gam(CCVD ~ s(Percent_both) + Age + Sex + Prev_DM + Prev_HTN + Prev_Dyslipidemia, 
                 data = nested, family = "binomial")

# Plot: Log Odds and OR
par(mfrow = c(1, 2))
plot(model.gam, ylab = "Log(OR)", xlab = "Treatment Duration (%)", main = "GAM Log(OR)")
plot(model.gam, ylab = "Odds Ratio", xlab = "Treatment Duration (%)", trans = exp, main = "GAM OR")
