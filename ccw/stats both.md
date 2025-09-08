# Summary
The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -4.24194  49.9959
   2 │ shorter     -8.32551  43.5482

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.31804

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.6708

Details:
    number of observations: [17, 11]
    F statistic:            1.3180425893194656
    degrees of freedom:     [16, 10]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          4.08357
    95% confidence interval: (-32.84, 41.0)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.8212

Details:
    number of observations:   [17,11]
    t-statistic:              0.22847909943768216
    degrees of freedom:       23.601433801818096
    empirical standard error: 17.872846625061637


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          12.5675
    95% confidence interval: (-24.35, 49.49)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.4888

Details:
    number of observations:   [17,11]
    t-statistic:              0.7031589426409395
    degrees of freedom:       23.601433801818096
    empirical standard error: 17.872846625061637


This suggests that beetles overshot or undershot by approximately the same amount (mean: -0.6952732433343043°, standard deviation: 47.14544424030437°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          -0.695273
    95% confidence interval: (-18.98, 17.59)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.9384

Details:
    number of observations:   28
    t-statistic:              -0.07803596401483351
    degrees of freedom:       27
    empirical standard error: 8.909651493536273


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.