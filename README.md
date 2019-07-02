This package provides new equivalence test for approximate independence in two-way contingency tables.
The package is based on the article:

Vladimir Ostrovski "Testing equivalence to families of multinomial distributions
with application to the independence model".
Statistics and Probability Letters 139 (2018) 61–66
https://doi.org/10.1016/j.spl.2018.03.014

Three examples are available in the script "examples.R“.
The program is rewritten in R and optimized for better performance and understanding. 
Hence the results for real data sets deviate slightly from these, published in the article. 

The goal is to show that two discrete random variables are approximately independent distributed. 
The approximate independence can be important for some applications or also greatly simplify calculations. 

The project contains an equivalence test for approximate row column independence
in two way contingency tables. 
The test statistic is the minimum of the Euclidean distance 
between the contingency table and a product measure.
The minimum distance is computed using constrained optimization.
The asymptotic test for the test statistic is the function "asymptotic_test".
The bootstrap test is the function "bootstrap_test".
The tests return the minimum tolerance parameter epsilon,
for which the approximate independence can be shown.

The asymptotic test is based on the asymptotic distribution of the test statistic. 
The asymptotic test needs some sufficiently large number of the observations
in any cell of the contingency table.
It should be used carefully because the test is approximate 
and may be anti-conservative at some points. 
In order to obtain a conservative test reducing of alpha  (usually halving) or
slight shrinkage of the tolerance parameter epsilon may be appropriate.
We prefer the slight shrinkage of the tolerance parameter 
because it is more effective and the significance level remains unchanged.

The bootstrap test is based on the re-sampling method called bootstrap.
The bootstrap test is more precise and reliable than the asymptotic test.
However, it should be used carefully because the test is approximate 
and may be anti-conservative. 
In order to obtain a conservative test reducing of alpha
(usually halving) or slight shrinkage of the tolerance parameter epsilon
may be appropriate. We prefer the slight shrinkage of the tolerance parameter 
because it is more effective and the significance level remains unchanged.
