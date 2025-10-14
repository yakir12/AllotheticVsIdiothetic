# Summary
Out of a total of 262 runs, 0 turned both directions and 0 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 142 turned towards the shorter direction and 120 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -18.9991  51.1592
   2 │ shorter      13.9154  42.4565

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.45197

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0336

Details:
    number of observations: [120, 142]
    F statistic:            1.4519715627111736
    degrees of freedom:     [119, 141]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -32.9146
    95% confidence interval: (-44.49, -21.34)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-07

Details:
    number of observations:   [120,142]
    t-statistic:              -5.6033758298834915
    degrees of freedom:       231.6141716047414
    empirical standard error: 5.874060770230635


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          5.08371
    95% confidence interval: (-6.49, 16.66)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3877

Details:
    number of observations:   [120,142]
    t-statistic:              0.8654507160392173
    degrees of freedom:       231.6141716047414
    empirical standard error: 5.874060770230635


This suggests that beetles overshot or undershot by approximately the same amount (mean: 16.243846861297758°, standard deviation: 46.621313726911104°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          16.2438
    95% confidence interval: (10.57, 21.92)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-07

Details:
    number of observations:   262
    t-statistic:              5.639687304893081
    degrees of freedom:       261
    empirical standard error: 2.8802743810289524


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.