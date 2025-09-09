# Summary
Out of a total of 272 runs, 0 turned both directions and 10 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 160 turned towards the shorter direction and 112 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -16.1135  51.4296
   2 │ shorter      10.1031  36.4216

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.99393

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-04

Details:
    number of observations: [112, 160]
    F statistic:            1.9939255130007543
    degrees of freedom:     [111, 159]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -26.2166
    95% confidence interval: (-37.36, -15.07)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   [112,160]
    t-statistic:              -4.64123453291693
    degrees of freedom:       186.56500011134597
    empirical standard error: 5.648624544009979


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          6.01031
    95% confidence interval: (-5.133, 17.15)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.2887

Details:
    number of observations:   [112,160]
    t-statistic:              1.0640310670489224
    degrees of freedom:       186.56500011134597
    empirical standard error: 5.648624544009979


This suggests that beetles overshot or undershot by approximately the same amount (mean: 12.577974001430485°, standard deviation: 43.2488006953867°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          12.578
    95% confidence interval: (7.415, 17.74)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   272
    t-statistic:              4.796462748591556
    degrees of freedom:       271
    empirical standard error: 2.622343727181851


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.