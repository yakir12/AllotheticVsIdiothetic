# Summary
The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction  μ         σ
     │ String     Float64   Float64
─────┼──────────────────────────────
   1 │ shorter     15.5727  51.6761
   2 │ longer     -10.7722  35.6455

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          2.10169

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-04

Details:
    number of observations: [113, 159]
    F statistic:            2.1016936368613517
    degrees of freedom:     [112, 158]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          0.459805
    95% confidence interval: (0.2662, 0.6534)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   [113,159]
    t-statistic:              4.684822364265552
    degrees of freedom:       185.51572656338445
    empirical standard error: 0.09814783349085907


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          0.0837849
    95% confidence interval: (-0.1098, 0.2774)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3944

Details:
    number of observations:   [113,159]
    t-statistic:              0.8536602684107915
    degrees of freedom:       185.51572656338445
    empirical standard error: 0.09814783349085907


This suggests that beetles overshot or undershot by approximately the same amount (mean: 12.766521236357251°, standard deviation: 43.01226310593082°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          0.222818
    95% confidence interval: (0.1332, 0.3124)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   272
    t-statistic:              4.89513564069442
    degrees of freedom:       271
    empirical standard error: 0.045518213580819686


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.