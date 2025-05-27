library(betareg)
setwd("/home/yakir/new_projects/Elin/code/indoors/analyse")
data <- read.csv("data.csv")

m <- betareg(mean_resultant_vector ~ r*condition, data)

r <- seq(0, 30, length = 100)
condition <- c("remain", "dark", "dark induced")
newdata <- expand.grid(r, condition)
colnames(newdata) <- c("r", "condition")
f <- predict(m, newdata, type = "quantile", at = c(0.025, 0.5, 0.975))
f <- as.data.frame(f)
f$r = newdata$r
f$condition = newdata$condition

write.csv(f, "fitted mean_resultant_vector.csv")


library(emmeans)
emm_betareg <- emmeans(m, , type = 'quntile')
