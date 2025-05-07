# ------------------------------------------------------------------------------
# This is a code sample from my previous work as a Clinical Data Analyst.
# It supports a published research study on the mass screening and associated 
# rules analysis for comorbidities of psoriasis in the Korean population. 
# Meta-analysis and association rules analysis are the major methods here.
#
# The original script included over 700 lines of code. Here, I’ve reorganized it 
# to focus on the key components that demonstrate my applied data skills:
#   - Table 1: Baseline demographics
#   - Table 2: Odds ratios of comorbid diseases (multivariable analysis)
#   - Table 3: Association rules of comorbidities 
#   - Figure 3: Network graph of comorbidities (support, lift, and confidence)
#
# This work was published in a high-impact dermatology journal.
# For holistic overview, please refer to the paper "Publication3_Psoriasis_Comorbid_MS.pdf"
# uploaded in my `Clinical Research Portfolio` GitHub repository. 
# ------------------------------------------------------------------------------

library(data.table);library(MatchIt);library(tableone)
library(meta);library(arules);library(arulesViz);library(magrittr)

setDTthreads(0)
set.seed(2)

# Load matched study population
info.mat <- fread("input.csv")

# ---Table 1: Baseline Demographics---
info.mat[, Age_group := ifelse(Age >= 80, 80, paste0(substr(Age, 1, 1), 
                              ifelse(as.integer(substr(Age, 2, 2)) < 5, 0, 5)))]
vars.tb1 <- c("Sex","Age", "Age_group", "Insurance")
tb1 <- CreateTableOne(vars = vars.tb1, strata = "Group", data = info.mat, 
                      factorVars = setdiff(vars.tb1, "Age"))
print(tb1)

# Load comorbidity codes per patient
info.code <- fread("input2.csv")
code <- setdiff(names(info.code), "JID")

# Logistic Regression for Each Code
info.est <- rbindlist(lapply(code, function(x) {
  data.logi <- merge(info.mat, info.code[, .(JID, dz = get(x))], by = "JID")
  model <- glm(dz ~ Group + Age + Sex + Insurance, data = data.logi, family = binomial)
  coef <- summary(model)$coefficients["GroupA", ]
  data.table(Code = x, Estimate = coef["Estimate"], Std.Error = coef["Std. Error"], Pvalue = coef["Pr(>|z|)"])
}))

# ---Table 2: Odds ratios of comorbid diseases---
info.est[, `:=`(OR = exp(Estimate),
                Lower = exp(Estimate - 1.96 * Std.Error),
                Upper = exp(Estimate + 1.96 * Std.Error))]
fwrite(info.est, "dataout/Table2_OddsRatios.csv")

# Meta-analysis
meta_input <- as.matrix(info.est[, .(Estimate, Std.Error)])
rownames(meta_input) <- info.est$Code
res.meta <- metagen(TE = meta_input[,1], seTE = meta_input[,2], sm = "OR", studlab = rownames(meta_input))
summary(res.meta)

# ---Table 3: Association rules of comorbidities---
code1 <- info.est[Pvalue < threshold & (OR < threshold | OR > threshold)]$Code
data1 <- as.matrix(info.code[, ..code1])
rownames(data1) <- info.code$JID
trans1 <- as(data1, "transactions")

rules <- apriori(trans1, parameter = list(supp = threshold, conf = threshold))
summary(rules)
rules.sorted <- sort(rules, by = "lift")
inspect(rules.sorted)
write(rules.sorted, file = "dataout/Table3_AssociationRules.txt")

# ---Figure 3: Network graph of comorbidities---
plot(rules.sorted, method = "graph", control = list(type = "items"))

