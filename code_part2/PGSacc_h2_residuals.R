# small analysis looking at the difference in residuals between groups for
# PGS accuracy ~ heritability

# Libraries and paths ####
library(tidyverse)
library(data.table)
library(psych)

# Code ####

# read table and filter to quantitative
tt <- as_tibble(fread("../generated_data/traits_table.txt"))
qtt <- tt %>% filter(PGS_trait_type == "quantitative",
                     GWAS_trait_type == "quantitative")
# make linear model
lm1 <- lm(pcor_United ~ ldpred2_h2, data=qtt)
ggplot(qtt, aes(x=ldpred2_h2, y=pcor_United)) +
  geom_point(aes(color=group_consolidated)) +
  geom_smooth(method='lm')

# check if residuals approximately normal and random
hist(lm1$residuals)
qtt$pcor_h2_resid <- lm1$residuals
# lowered variance at extreme values likely due to restriced axes? 
ggplot(qtt, aes(x=ldpred2_h2, y=pcor_h2_resid)) +
  geom_point(aes(color=group_consolidated)) +
  geom_smooth(method='lm')
mean(qtt$pcor_h2_resid)
# boxplots of residuals by group
ggplot(qtt, aes(x=group_consolidated, y=pcor_h2_resid)) +
  geom_boxplot(aes(fill=group_consolidated))


# two-sided t-test + Wilcoxon on lifestyle resids
resids <- qtt$pcor_h2_resid[qtt$group_consolidated=="lifestyle/psychological"]
t.test(x = resids)
hist(resids)
wilcox.test(x = resids)

# two-sided t-test + Wilcoxon on physical resids
resids <- qtt$pcor_h2_resid[qtt$group_consolidated=="physical measures"]
t.test(x = resids)
hist(resids)
wilcox.test(x = resids)

# two-sided t-test + Wilcoxon on biological resids
resids <- qtt$pcor_h2_resid[qtt$group_consolidated=="biological measures"]
t.test(x = resids)
hist(resids)
wilcox.test(x = resids)


##############

# difference of correlation analysis for portability

# Assume r12, r13, r23 are your correlation coefficients, and n is your sample size

r12 <- cor(qtt$portability_index, qtt$ldpred2_h2)
r13 <- cor(qtt$portability_index, qtt$pcor_United)
r23 <- cor(qtt$ldpred2_h2, qtt$pcor_United)

r_results <- paired.r(xy = r12, xz = r13, yz = r23, n = nrow(qtt))
r_results$p
