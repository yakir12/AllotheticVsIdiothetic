# Summary
Out of a total of 300 runs, 28 turned both directions and 14 turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, 160 turned towards the shorter direction and 140 turned towards the longer direction.

The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

2×3 DataFrame
 Row │ direction1  μ         σ
     │ String      Float64   Float64
─────┼───────────────────────────────
   1 │ longer      -15.7669  51.695
   2 │ shorter      10.4896  42.8447

An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

Variance F-test
---------------
Population details:
    parameter of interest:   variance ratio
    value under h_0:         1.0
    point estimate:          1.4558

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0220

Details:
    number of observations: [140, 160]
    F statistic:            1.4557987741135667
    degrees of freedom:     [139, 159]


Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          -26.2565
    95% confidence interval: (-37.14, -15.37)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   [140,160]
    t-statistic:              -4.749532510771893
    degrees of freedom:       270.78663879186877
    empirical standard error: 5.528227677089429


This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          5.27729
    95% confidence interval: (-5.606, 16.16)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3406

Details:
    number of observations:   [140,160]
    t-statistic:              0.954607502221411
    degrees of freedom:       270.78663879186877
    empirical standard error: 5.528227677089429


This suggests that beetles overshot or undershot by approximately the same amount (mean: 12.95233895249136°, standard deviation: 47.17471314221572°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          12.9523
    95% confidence interval: (7.592, 18.31)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-05

Details:
    number of observations:   300
    t-statistic:              4.755536949411282
    degrees of freedom:       299
    empirical standard error: 2.723633333160162


The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.