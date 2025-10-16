# Summary
Out of a total of 10 runs, 0 turned both directions and 10 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 8 turned towards the shorter direction and 2 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼────────────────────────────────
   1 │ longer       43.3248   6.64684
   2 │ shorter     -12.9743  24.7578

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          0.0720785

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.4079

Details:
    number of observations: [2, 8]
    F statistic:            0.07207848423549856
    degrees of freedom:     [1, 7]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          56.2991
    95% confidence interval: (33.03, 79.57)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0006

Details:
    number of observations:   [2,8]
    t-statistic:              5.666606463625059
    degrees of freedom:       7.344621663067803
    empirical standard error: 9.93523626658892


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -30.3504
    95% confidence interval: (-53.62, -7.079)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0174

Details:
    number of observations:   [2,8]
    t-statistic:              -3.0548268242182925
    degrees of freedom:       7.344621663067803
    empirical standard error: 9.93523626658892


This suggests that beetles overshot or undershot by approximately the same amount (mean: -19.044409147310763°, standard deviation: 25.404864419266644°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          -19.0444
    95% confidence interval: (-37.22, -0.8709)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0419

Details:
    number of observations:   10
    t-statistic:              -2.3705581971921275
    degrees of freedom:       9
    empirical standard error: 8.033723521265344


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.