# Summary
Out of a total of 300 runs, 28 turned both directions and 14 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 151 turned towards the shorter direction and 111 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -16.5081  51.4922
   2 │ shorter      10.9632  36.9529

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.94172

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0002

Details:
    number of observations: [111, 151]
    F statistic:            1.9417158494313922
    degrees of freedom:     [110, 150]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -27.4713
    95% confidence interval: (-38.79, -16.15)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   [111,151]
    t-statistic:              -4.787216622757536
    degrees of freedom:       189.1710538205338
    empirical standard error: 5.738474098684832


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          5.54493
    95% confidence interval: (-5.775, 16.86)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3351

Details:
    number of observations:   [111,151]
    t-statistic:              0.9662718092888242
    degrees of freedom:       189.1710538205338
    empirical standard error: 5.738474098684832


This suggests that beetles overshot or undershot by approximately the same amount (mean: 13.312382522442789°, standard deviation: 43.70108867158032°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          13.3124
    95% confidence interval: (7.996, 18.63)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   262
    t-statistic:              4.930763560666724
    degrees of freedom:       261
    empirical standard error: 2.6998622746053402


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.