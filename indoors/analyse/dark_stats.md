We investigated how the beetles' directedness deteriorated as a function of the following three factors: 1) light, 2) how we induced a dance in the beetles, and 3) the radial distance from the induction point (i.e. POI). As a measure for their directedness, we chose the mean resultant vector length (r-value). This length was calculated from the exit angles of all the beetles at a given radial distance from the POI. This meant that we got one r-value per group, which in our setup meant 4 r-values for each radial distance. Because r-values range between zero and one, we chose a beta regression to test the effect of our factors. While we could calculate the r-values for an arbitrary number of radial distances, each subsequent r-value would be highly correlated to the previous one. To avoid this auto-correlation issue we chose the most conservative number of radial distances, 3. Because this resulted in a sparse dataset (3 radial distances × 4 groups = 12 r-values), we bootstrapped (with replacement) our data 10,000 times, fitting a beta regression for each sample. Apart from collecting the p-values from each regression, we also predicted the response per group. Finally, we calculated some summary statistics for these p-values as well as the predicted responses. 

 Source             proportion  Q2.5    median  Q97.5   mode    mean   
───────────────────────────────────────────────────────────────────────
 dance_by: disrupt  98%         <0.001  <0.001  0.03    <0.001  0.0048
 dance_by: hold     97%         <0.001  <0.001  >0.05   <0.001  0.0076
 light: off         84%         <0.001  <0.001  >0.05   <0.001  0.04
 radial             90%         <0.001  <0.001  >0.05   <0.001  0.024

 Table X: Summary statistics for the p-values from the bootstrapped beta-regression. Each row in the table describes a different source of variation. Proportion shows the percent of p-values that were lower than 0.05. Q2.5 is the 2.5% lower quantile of the p-values, while Q97.5 is the 97.5% quantile. 

We can conclude that both induction methods (holding the ball and removing the beetle from the ball) are significantly different from the case were we did not induce the dance. Keeping the lights on was unsurprisingly also significant. Lastly, the radial distance from POI had also an effect (i.e. the slope was different from zero).
