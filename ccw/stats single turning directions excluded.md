# Summary
Out of a total of 28 runs, 28 turned both directions and 4 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 10 turned towards the shorter direction and 18 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ           σ
     │ String      Float64     Float64
─────┼─────────────────────────────────
   1 │ longer       -0.784312  52.9611
   2 │ shorter     -19.386     45.2451

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.37016

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.6457

Details:
    number of observations: [18, 10]
    F statistic:            1.3701594303608586
    degrees of freedom:     [17, 9]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          18.6017
    95% confidence interval: (-20.84, 58.05)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3382

Details:
    number of observations:   [18,10]
    t-statistic:              0.9796610813505572
    degrees of freedom:       21.363172107018677
    empirical standard error: 18.9878519442264


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          20.1703
    95% confidence interval: (-19.28, 59.62)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3000

Details:
    number of observations:   [18,10]
    t-statistic:              1.0622730490047354
    degrees of freedom:       21.363172107018677
    empirical standard error: 18.9878519442264


This suggests that beetles overshot or undershot by approximately the same amount (mean: -6.419360729982005°, standard deviation: 50.450716976770835°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          -6.41936
    95% confidence interval: (-25.98, 13.14)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.5065

Details:
    number of observations:   28
    t-statistic:              -0.6732920000073118
    degrees of freedom:       27
    empirical standard error: 9.53428932753143


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.