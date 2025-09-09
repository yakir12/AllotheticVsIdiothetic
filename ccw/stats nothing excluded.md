# Summary
Out of a total of 300 runs, 28 turned both directions and 14 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 171 turned towards the shorter direction and 129 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ          σ
     │ String      Float64    Float64
─────┼────────────────────────────────
   1 │ longer      -14.549    51.2095
   2 │ shorter       8.91767  37.0515

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.91025

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-04

Details:
    number of observations: [129, 171]
    F statistic:            1.910247920272187
    degrees of freedom:     [128, 170]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -23.4667
    95% confidence interval: (-33.96, -12.97)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-04

Details:
    number of observations:   [129,171]
    t-statistic:              -4.406780019068167
    degrees of freedom:       222.88782444543696
    empirical standard error: 5.325126193765271


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          5.63132
    95% confidence interval: (-4.863, 16.13)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.2914

Details:
    number of observations:   [129,171]
    t-statistic:              1.0574994539045113
    degrees of freedom:       222.88782444543696
    empirical standard error: 5.325126193765271


This suggests that beetles overshot or undershot by approximately the same amount (mean: 11.339137591919107°, standard deviation: 43.714627666118965°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          11.3391
    95% confidence interval: (6.372, 16.31)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-04

Details:
    number of observations:   300
    t-statistic:              4.492766717178301
    degrees of freedom:       299
    empirical standard error: 2.5238652050558046


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.