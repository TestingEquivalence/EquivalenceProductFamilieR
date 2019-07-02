 source("asymptotic_test.R")
source("bootstrap_test.R")

# This project contains asymptotic and bootstrap tests,
# which are based on the article:
# Vladimir Ostrovski "Testing equivalence to families of multinomial distributions
# with application to the independence model".
# Statistics and Probability Letters 139 (2018) 61–66
# https://doi.org/10.1016/j.spl.2018.03.014

# The project contains an equivalence test for approximate row column independence
# in two way contingency tables. 
# The test statistic is the minimum of the Euclidean distance 
# between the contingency table and a product measure.
# The minimum distance is computed using constrained optimization.
# The asymptotic test for the test statistic is the function "asymptotic_test".
# The bootstrap test is the function "bootstrap_test".
# The tests return the minimum tolerance parameter epsilon,
# for which the approximate independence can be shown.

# The real data examples, which are considered in the paper, are presented below.
# -------------------------------------------------------------------------------

# significance level
alpha=0.05

#table of tests results
results=matrix(data=NA, nrow=3,ncol=2)
rownames(results)=c("Gender and nitrendipine",
                    "Eye color and hair color",
                    "Children number and income")
colnames(results)=c("asymptotic","bootstrap")

# Nitren data set
# -------------------------------------------------------------------------------
# Contingency table relating gender and treatment outcome on nitrendipine mono-therapy
# in patients suffering from mild arterial hypertension

nitren_data  = matrix(data=c(9, 13, 13, 48,24, 18, 20, 72),
                      nrow=2, ncol=4, byrow=TRUE)

results[1,1]=asymptotic_test(nitren_data,alpha)
results[1,2]=bootstrap_test(nitren_data,alpha)

# Cross-classification of eye color and hair color
# -------------------------------------------------------------------------------

eye_hair_color_data=matrix(
  data=c(68, 119, 26, 7,20, 84, 17, 94,15, 54, 14, 10,5, 29, 14, 16),
  nrow=4,ncol=4, byrow=TRUE)

results[2,1]=asymptotic_test(eye_hair_color_data,alpha)
results[2,2]=bootstrap_test(eye_hair_color_data,alpha)

# Cross-classification of number of children by annual income
# -------------------------------------------------------------------------------

children_income_data=matrix(
  data=c(2161, 3577, 2184, 1636,
         2755, 5081, 2222, 1052,
         936, 1753, 640, 306,
         225, 419, 96, 38,
         39, 98, 31, 14),
  nrow=5, ncol=4, byrow = TRUE)


results[3,1]=asymptotic_test(children_income_data,alpha)
results[3,2]=bootstrap_test(children_income_data,alpha)
