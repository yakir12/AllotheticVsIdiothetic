# Summary
Out of a total of 272 runs, 0 turned both directions and 10 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 150 turned towards the shorter direction and 122 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -17.9774  51.3568
   2 │ shorter      12.4813  42.0871

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.48901

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0209

Details:
    number of observations: [122, 150]
    F statistic:            1.4890091555897875
    degrees of freedom:     [121, 149]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -30.4587
    95% confidence interval: (-41.85, -19.07)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-06

Details:
    number of observations:   [122,150]
    t-statistic:              -5.268148372542848
    degrees of freedom:       232.86613093996755
    empirical standard error: 5.781679759390181


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          5.49613
    95% confidence interval: (-5.895, 16.89)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3428

Details:
    number of observations:   [122,150]
    t-statistic:              0.9506107944728073
    degrees of freedom:       232.86613093996755
    empirical standard error: 5.781679759390181


This suggests that beetles overshot or undershot by approximately the same amount (mean: 14.946484508040088°, standard deviation: 46.465445520633374°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          14.9465
    95% confidence interval: (9.4, 20.49)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-06

Details:
    number of observations:   272
    t-statistic:              5.305097899551337
    degrees of freedom:       271
    empirical standard error: 2.817381467984623


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.