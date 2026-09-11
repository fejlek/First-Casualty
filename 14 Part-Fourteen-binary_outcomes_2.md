# The First Casualty of Statistics: Part Fourteen
<big>**Binary Outcomes II**</big>

<br/>
Jiří Fejlek

2026-09-09
<br/>

<br/> We will complete our discussion of estimating marginal effects for
binary outcomes by demonstrating how to estimate them under observed
confounding. We will see that the computations are very similar to those
for continuous outcomes. <br/>


## Table of Contents

- [Estimating Marginal Risk Difference under
  Confounding](#estimating-marginal-risk-difference-under-confounding)
- [Estimating Marginal Risk Ratio under
  Confounding](#estimating-marginal-risk-ratio-under-confounding)
- [Estimating Marginal Odds Ratio under
  Confounding](#estimating-marginal-odds-ratio-under-confounding)
- [Conditional Odds Ratio under Confounding](#conditional-odds-ratio-under-confounding)
- [References](#references)

``` r
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(dagitty)
library(ggdag)
library(modelsummary)
library(GGally)
library(estimatr)
library(WeightIt)
library(MatchIt)
library(cobalt)
library(mgcv)
library(marginaleffects)
```

We will complete our discussion of estimating marginal effects for
binary outcomes by demonstrating how to estimate them under observed
confounding. We will see that the computations are very similar to those
for continuous outcomes.

## Estimating Marginal Risk Difference under Confounding

Let’s start with the marginal risk difference, which corresponds to the
standard ATE for continuous outcomes. We will use the same example as in
the previous part

``` r
set.seed(123)
n_pop <- 1000

X1 <- runif(n_pop, -5, 5)
X2 <- rnorm(n_pop, 2.5, 1)
X3 <- (runif(n_pop, 0, 1)>0.75) + 0

p0 <- plogis(-2 + 0.5*X1 + 0.75*X2 + 2*X3)
p1 <- plogis(-2 + 0.5*X1 + 0.75*X2 + 2*X3 + rep(1,n_pop))

y0 <- (runif(n_pop,0,1) > (1-p0)) + 0
y1 <- (runif(n_pop,0,1) > (1-p1)) + 0
```

but we will assume that treatment assignment depends on $`X`$. Let’s
estimate the marginal risk difference using a naive approach that
ignores confounding and regression adjustment (e.g., linear or logistic
regression). We use a linear regression model with treatment
interactions, since as we discussed earlier, the treatment effect is
usually heterogeneous on the risk difference scale.

``` r
set.seed(456)

n_sim <- 100
estimates_RD <- matrix(0,n_sim,3)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  estimates_RD[i,1] <- mean(y[treatment == 1]) - mean(y[treatment == 0])
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2)
  estimates_RD[i,2] <- avg_comparisons(ols_model, variables = "treatment")$estimate

  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  estimates_RD[i,3] <- avg_comparisons(logist_model, variables = "treatment")$estimate
 
}

results <- rbind(
  apply(estimates_RD,2,mean),
  apply(estimates_RD,2,sd))
```

``` r
results_RD <- rbind(
  c(mean(p1) - mean(p0), NA),
  c(mean(y1) - mean(y0), NA),
  t(results)
)

colnames(results_RD) <-  c('RD Estimate', 'sd')
rownames(results_RD) <-  c('True Marginal RD (mean)', 'True Marginal RD (population)', 'Naive RD', 'Linear Regression', 'Logistic Regression')
results_RD
```

    ##                               RD Estimate         sd
    ## True Marginal RD (mean)         0.1458635         NA
    ## True Marginal RD (population)   0.1710000         NA
    ## Naive RD                        0.2103759 0.02045429
    ## Linear Regression               0.1719565 0.01586666
    ## Logistic Regression             0.1724630 0.01520552

We see that the naive risk difference is biased due to confounding. Both
adjusted models estimated the true effect. However, we know more methods
based on weighting and matching. So let’s explore them. We will start
with propensity score weighting.

We again need to remember that the risk difference is just a usual ATE;
hence, the Hájek estimator and the augmented IPW will be the same as for
a continuous outcome.

``` r
set.seed(456)

n_sim <- 100
estimates_RD <- matrix(0,n_sim,7)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  data2_0 <- data2
  data2_1 <- data2
  data2_0$treatment <- 0
  data2_1$treatment <- 1
  
  
  prop_scores_model_ATE <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "glm", estimand = "ATE")
  
  estimates_RD[i,1] <- sum((treatment*y/prop_scores_model_ATE$ps))/sum(((treatment/prop_scores_model_ATE$ps))) - sum(((1-treatment)*y/(1-prop_scores_model_ATE$ps)))/sum((((1-treatment)/(1-prop_scores_model_ATE$ps))))
  
  ols_unadj_model <- lm(y ~ treatment, data = data2, weights = prop_scores_model_ATE$weights)
  estimates_RD[i,2] <- avg_comparisons(ols_unadj_model, variables = "treatment")$estimate
  
  logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_RD[i,3] <- avg_comparisons(logist_unadj_model, variables = "treatment")$estimate
  
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2)
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2)
  
  
  pred0 <- predict(poisson_model, newdata = data2_0, type = 'response')
  pred1 <- predict(poisson_model, newdata = data2_1, type = 'response')
  
  
  estimates_RD[i,4] <-  mean(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps) - mean(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))
  
  
  pred0 <- predict(logist_model, newdata = data2_0, type = 'response')
  pred1 <- predict(logist_model, newdata = data2_1, type = 'response')
  
  
  estimates_RD[i,5] <-  mean(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps) - mean(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))
  
  

  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_RD[i,6] <- avg_comparisons(ols_model, variables = "treatment")$estimate
  
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_RD[i,7] <- avg_comparisons(logist_model, variables = "treatment")$estimate
}  
  
results <- rbind(
  apply(estimates_RD,2,mean),
  apply(estimates_RD,2,sd))
```

``` r
results_RD2 <- rbind(
  t(results)
)

colnames(results_RD2) <-  c('RD Estimate', 'sd')
rownames(results_RD2) <-  c('IPW (Hájek)', 'IPW (Unadj. Linear Regression)', 'IPW (Unadj. Logistic Regression)', 'AIPW (Linear Regression)', 'AIPW (Logistic Regression)', 'IPWRA (Linear Regression)', 'IPWRA (Logistic Regression)')
results_RD2
```

    ##                                  RD Estimate         sd
    ## IPW (Hájek)                        0.1727857 0.01610720
    ## IPW (Unadj. Linear Regression)     0.1727857 0.01610720
    ## IPW (Unadj. Logistic Regression)   0.1727857 0.01610720
    ## AIPW (Linear Regression)           0.1727118 0.01748639
    ## AIPW (Logistic Regression)         0.1728885 0.01537158
    ## IPWRA (Linear Regression)          0.1729819 0.01602042
    ## IPWRA (Logistic Regression)        0.1728706 0.01536062

We observe that the Hájek estimator equals both the unadjusted
regression model and the model weighted by propensity scores, as is the
case for a continuous outcome. We can also use the doubly robust
estimators AIPW and IPWRA, which rely on both adjusted linear regression
and logistic regression, and all provide very similar estimates in this
case.

We can also move beyond simple inverse propensity scores weights.

``` r
set.seed(456)

n_sim <- 100
estimates_RD <- matrix(0,n_sim,9)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  
  prop_scores_model_cbps <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "cbps", estimand = "ATE", moments = 2, over = FALSE)
  
  ols_unadj_model <- lm(y ~ treatment, data = data2, weights = prop_scores_model_cbps$weights)
  estimates_RD[i,1] <- avg_comparisons(ols_unadj_model, variables = "treatment")$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2, weights = prop_scores_model_cbps$weights)
  estimates_RD[i,2] <- avg_comparisons(ols_model, variables = "treatment")$estimate
  
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_cbps$weights)
  estimates_RD[i,3] <- avg_comparisons(logist_model, variables = "treatment")$estimate
  
  
  prop_scores_model_ebal <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "ebal", estimand = "ATE", moments = 2)
  
  ols_unadj_model <- lm(y ~ treatment, data = data2, weights = prop_scores_model_ebal$weights)
  estimates_RD[i,4] <- avg_comparisons(ols_unadj_model, variables = "treatment")$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2, weights = prop_scores_model_ebal$weights)
  estimates_RD[i,5] <- avg_comparisons(ols_model, variables = "treatment")$estimate
  
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ebal$weights)
  estimates_RD[i,6] <- avg_comparisons(logist_model, variables = "treatment")$estimate
  
  
  prop_scores_model_eng <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "energy", estimand = "ATE")
  
  ols_unadj_model <- lm(y ~ treatment, data = data2, weights = prop_scores_model_eng$weights)
  estimates_RD[i,7] <- avg_comparisons(ols_unadj_model, variables = "treatment")$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2, weights = prop_scores_model_eng$weights)
  estimates_RD[i,8] <- avg_comparisons(ols_model, variables = "treatment")$estimate
  
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_eng$weights)
  estimates_RD[i,9] <- avg_comparisons(logist_model, variables = "treatment")$estimate
}  
  
results <- rbind(
  apply(estimates_RD,2,mean),
  apply(estimates_RD,2,sd))
```

``` r
results_RD3 <- rbind(
  t(results)
)

colnames(results_RD3) <-  c('RD Estimate', 'sd')
rownames(results_RD3) <-  c('CBPS (Unadjusted)', 'CBPS (Linear Regression)', 'CBPS (Logistic Regression)', 'Entropy Balancing (Unadjusted)', 'Entropy Balancing (Linear Regression)', 'Entropy Balancing (Logistic Regression)', 'Energy Balancing (Unadjusted)', 'Energy Balancing (Linear Regression)', 'Energy Balancing (Logistic Regression)')
results_RD3
```

    ##                                         RD Estimate         sd
    ## CBPS (Unadjusted)                         0.1729813 0.01592757
    ## CBPS (Linear Regression)                  0.1730122 0.01589331
    ## CBPS (Logistic Regression)                0.1728386 0.01549297
    ## Entropy Balancing (Unadjusted)            0.1725674 0.01577614
    ## Entropy Balancing (Linear Regression)     0.1725674 0.01577613
    ## Entropy Balancing (Logistic Regression)   0.1727843 0.01543630
    ## Energy Balancing (Unadjusted)             0.1717764 0.01886381
    ## Energy Balancing (Linear Regression)      0.1717291 0.01886609
    ## Energy Balancing (Logistic Regression)    0.1716928 0.01890319

We can also consider matching. Let us start with the propensity score
matching (PSM). We need to remember that PSM estimates ATT.

``` r
set.seed(456)

n_sim <- 100
estimates_RD <- matrix(0,n_sim,6)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  estimates_RD[i,1] <- mean(y1[treatment == 1]) - mean(y0[treatment == 1])
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = data2)
  estimates_RD[i,2] <- avg_comparisons(ols_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2)
  estimates_RD[i,3] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "nearest", m.order = "largest", distance = 'glm', replace = TRUE)
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,4] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = match_data, weights = weights)
  estimates_RD[i,5] <- avg_comparisons(ols_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,6] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
}
```

``` r
results <- rbind(
  apply(estimates_RD,2,mean),
  apply(estimates_RD,2,sd))
```

``` r
results_RD_ATT <- t(results)

colnames(results_RD_ATT) <-  c('RD Estimate', 'sd')
rownames(results_RD_ATT) <-  c('True Marginal RD for Treated (mean)', 'Linear Regression', 'Logistic Regression', 'PSM (Unadjusted)', 'PSM (Linear Regression)', 'PSM (Logistic Regression)')
results_RD_ATT
```

    ##                                     RD Estimate         sd
    ## True Marginal RD for Treated (mean)   0.1670557 0.02306800
    ## Linear Regression                     0.1717863 0.01658405
    ## Logistic Regression                   0.1693859 0.01559880
    ## PSM (Unadjusted)                      0.1702696 0.03263405
    ## PSM (Linear Regression)               0.1714441 0.02883986
    ## PSM (Logistic Regression)             0.1720494 0.02805571

Again, we can also consider other matching methods.

``` r
set.seed(456)

n_sim <- 100
estimates_RD <- matrix(0,n_sim,9)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "genetic", estimand = "ATT", replace = TRUE)
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,1] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = match_data, weights = weights)
  estimates_RD[i,2] <- avg_comparisons(ols_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,3] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "full", estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,4] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = match_data, weights = weights)
  estimates_RD[i,5] <- avg_comparisons(ols_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,6] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "cardinality",  ratio = 1, tols = .05, estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,7] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
  
  ols_model <- lm(y ~ treatment*(X1 + X2 + X3), data = match_data, weights = weights)
  estimates_RD[i,8] <- avg_comparisons(ols_model, variables = "treatment", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RD[i,9] <- avg_comparisons(logist_model, variables = "treatment", newdata = subset(treatment == 1))$estimate
}

results <- rbind(
  apply(estimates_RD,2,mean),
  apply(estimates_RD,2,sd))
```

``` r
results_RD_ATT2 <- t(results)

colnames(results_RD_ATT2) <-  c('RD Estimate', 'sd')
rownames(results_RD_ATT2) <-  c('Genetic Matching (Unadjusted)', 'Genetic Matching (Linear Regression)', 'Genetic Matching (Logistic Regression)', 'Full Matching (Unadjusted)', 'Full Matching (Linear Regression)', 'Full Matching (Logistic Regression)', 'Cardinality Matching (Unadjusted)', 'Cardinality Matching (Linear Regression)', 'Cardinality Matching (Logistic Regression)')
results_RD_ATT2
```

    ##                                            RD Estimate         sd
    ## Genetic Matching (Unadjusted)                0.1674638 0.02791502
    ## Genetic Matching (Linear Regression)         0.1673533 0.02791580
    ## Genetic Matching (Logistic Regression)       0.1680235 0.02788251
    ## Full Matching (Unadjusted)                   0.1725953 0.02460180
    ## Full Matching (Linear Regression)            0.1728639 0.02317221
    ## Full Matching (Logistic Regression)          0.1733122 0.02295646
    ## Cardinality Matching (Unadjusted)            0.1750795 0.01971016
    ## Cardinality Matching (Linear Regression)     0.1720869 0.02035655
    ## Cardinality Matching (Logistic Regression)   0.1696858 0.02043041

We see that all methods provide very similar results in this case.

## Estimating Marginal Risk Ratio under Confounding

Let’s move to the marginal risk ratio. We will consider Poisson models
(with interactions) and logistic models as the most natural models for
modeling ratios of risks.

``` r
set.seed(456)

n_sim <- 100
estimates_RR <- matrix(0,n_sim,3)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  estimates_RR[i,1] <- mean(y[treatment == 1])/mean(y[treatment == 0])
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2)
  estimates_RR[i,2] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio")$estimate

  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  estimates_RR[i,3] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio")$estimate
 
}

results <- rbind(
  apply(estimates_RR,2,mean),
  apply(estimates_RR,2,sd))
```

``` r
results_RR <- rbind(
  c(mean(p1)/mean(p0), NA),
  c(mean(y1)/mean(y0), NA),
  t(results)
)

colnames(results_RR) <-  c('RR Estimate', 'sd')
rownames(results_RR) <-  c('True Marginal RR (mean)', 'True Marginal RR (population)', 'Naive RR', 'Poisson Regression', 'Logistic Regression')
results_RR
```

    ##                               RR Estimate         sd
    ## True Marginal RR (mean)          1.260637         NA
    ## True Marginal RR (population)    1.313761         NA
    ## Naive RR                         1.399261 0.04438836
    ## Poisson Regression               1.317476 0.03587456
    ## Logistic Regression              1.317319 0.03114207

Again, using an unadjusted (naive) estimate results in a biased
estimator.

We can also use propensity score weights to estimate the marginal risk
ratio (Boughdiri et al. 2024).

``` math
\text{RR}_\text{Hájek} = \frac{\sum_{i=1}^n T_iY_i}{\sum_{i = 1}^n \hat{e}(X_i)T_i} / \frac{\sum_{i=1}^n (1-T_i)Y_i}{\sum_{i = 1}^n \hat{e}(X_i)(1-T_i)}
```

``` math
\text{RR}_\text{AIPW} = \left(\frac{\sum_{i=1}^n T_i(Y_i-\mu_1(X_i, \hat \beta_1))}{\sum_{i = 1}^n \hat{e}(X_i)T_i} + \mu_1(X_i, \hat \beta_1)\right) / \left(\frac{\sum_{i=1}^n (1-T_i)(Y_i - \mu_0(X_i, \hat \beta_0))}{\sum_{i = 1}^n \hat{e}(X_i)(1-T_i)} + \mu_0(X_i, \hat \beta_0)\right)
```

``` r
set.seed(456)

n_sim <- 100
estimates_RR <- matrix(0,n_sim,7)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  data2_0 <- data2
  data2_1 <- data2
  data2_0$treatment <- 0
  data2_1$treatment <- 1
  
  prop_scores_model_ATE <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "glm", estimand = "ATE")
  
  estimates_RR[i,1] <- (sum((treatment*y/prop_scores_model_ATE$ps))/sum(((treatment/prop_scores_model_ATE$ps))))/(sum(((1-treatment)*y/(1-prop_scores_model_ATE$ps)))/sum((((1-treatment)/(1-prop_scores_model_ATE$ps)))))
  
  poisson_unadj_model <- glm(y ~ treatment, data = data2, family = poisson(link = 'log'), weights = prop_scores_model_ATE$weights)
  estimates_RR[i,2] <- avg_comparisons(poisson_unadj_model, variables = "treatment", comparison = "ratio")$estimate
  
  logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_RR[i,3] <- avg_comparisons(logist_unadj_model, variables = "treatment", comparison = "ratio")$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2)
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  
  pred0 <- predict(poisson_model, newdata = data2_0, type = 'response')
  pred1 <- predict(poisson_model, newdata = data2_1, type = 'response')
  
  estimates_RR[i,4] <-  sum(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps)/sum(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))
  
  
  pred0 <- predict(logist_model, newdata = data2_0, type = 'response')
  pred1 <- predict(logist_model, newdata = data2_1, type = 'response')
  
  estimates_RR[i,5] <-  sum(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps)/sum(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))
  
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2, weights = prop_scores_model_ATE$weights)
  logist_model <-glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  
  estimates_RR[i,6] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio")$estimate
  estimates_RR[i,7] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio")$estimate
}  
  
results <- rbind(
  apply(estimates_RR,2,mean),
  apply(estimates_RR,2,sd))
```

``` r
results_RR2 <- rbind(
  t(results)
)

colnames(results_RR2) <-  c('RD Estimate', 'sd')
rownames(results_RR2) <-  c('IPW (Hájek)', 'IPW (Unadj. Poisson Regression)', 'IPW (Unadj. Logistic Regression)','AIPW (Poisson Regression)', 'AIPW (Logistic Regression)', 'IPWRA (Poisson Regression)', 'IPWRA (Logistic Regression)')
results_RR2
```

    ##                                  RD Estimate         sd
    ## IPW (Hájek)                         1.318116 0.03273140
    ## IPW (Unadj. Poisson Regression)     1.318116 0.03273140
    ## IPW (Unadj. Logistic Regression)    1.318116 0.03273140
    ## AIPW (Poisson Regression)           1.317976 0.03577273
    ## AIPW (Logistic Regression)          1.318220 0.03127328
    ## IPWRA (Poisson Regression)          1.318058 0.03568999
    ## IPWRA (Logistic Regression)         1.318180 0.03125198

We observe that, again, the Hájek estimator is equivalent to unadjusted
regression with inverse propensity score weights. Other than that, all
methods yield very similar results.

Other weighting methods work as well.

``` r
set.seed(456)

n_sim <- 100
estimates_RR <- matrix(0,n_sim,9)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  prop_scores_model_cbps <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "cbps", estimand = "ATE", moments = 2, over = FALSE)
  

  poisson_unadj_model <- glm(y ~ treatment, data = data2, family = poisson(link = 'log'), weights = prop_scores_model_cbps$weights)
  estimates_RR[i,1] <- avg_comparisons(poisson_unadj_model, variables = "treatment", comparison = "ratio")$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2, weights = prop_scores_model_cbps$weights)
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_cbps$weights)
  
  estimates_RR[i,2] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio")$estimate
  estimates_RR[i,3] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio")$estimate
  
  
  prop_scores_model_eng <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "energy", estimand = "ATE")
  
   poisson_unadj_model <- glm(y ~ treatment, data = data2, family = poisson(link = 'log'), weights = prop_scores_model_eng$weights)
  estimates_RR[i,4] <- avg_comparisons(poisson_unadj_model, variables = "treatment", comparison = "ratio")$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2, weights = prop_scores_model_eng$weights)
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_eng$weights)
  
  estimates_RR[i,5] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio")$estimate
  estimates_RR[i,6] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio")$estimate
  
  
  prop_scores_model_ebal <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "ebal", estimand = "ATE", moments = 2)
  
   poisson_unadj_model <- glm(y ~ treatment, data = data2, family = poisson(link = 'log'), weights = prop_scores_model_ebal$weights)
  estimates_RR[i,7] <- avg_comparisons(poisson_unadj_model, variables = "treatment", comparison = "ratio")$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2, weights = prop_scores_model_ebal$weights)
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ebal$weights)
  
  estimates_RR[i,8] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio")$estimate
  estimates_RR[i,9] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio")$estimate
}  
  
results <- rbind(
  apply(estimates_RR,2,mean),
  apply(estimates_RR,2,sd))
```

``` r
estimates_RR3 <- rbind(
  t(results)
)

colnames(estimates_RR3) <-  c('RD Estimate', 'sd')
rownames(estimates_RR3) <-  c('IPW CBPS (Unadjusted)', 'IPW CBPS (Poisson Regression)', 'IPW CBPS (Logistic Regression)','Entropy Balancing (Unadjusted)', 'Entropy Balancing (Poisson Regression)', 'Entropy Balancing (Logistic Regression)','Energy Balancing (Unadjusted)', 'Energy Balancing (Poisson Regression)', 'Energy Balancing (Logistic Regression)')
estimates_RR3
```

    ##                                         RD Estimate         sd
    ## IPW CBPS (Unadjusted)                      1.318477 0.03258232
    ## IPW CBPS (Poisson Regression)              1.318212 0.03407096
    ## IPW CBPS (Logistic Regression)             1.318141 0.03158517
    ## Entropy Balancing (Unadjusted)             1.315243 0.03907935
    ## Entropy Balancing (Poisson Regression)     1.314334 0.03894650
    ## Entropy Balancing (Logistic Regression)    1.315097 0.03915828
    ## Energy Balancing (Unadjusted)              1.317979 0.03213186
    ## Energy Balancing (Poisson Regression)      1.318765 0.03414451
    ## Energy Balancing (Logistic Regression)     1.318011 0.03144522

Lastly, we use matching to estimate the marginal risk ratio for the
treated.

``` r
set.seed(456)

n_sim <- 100
estimates_RR <- matrix(0,n_sim,12)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  
  estimates_RR[i,1] <- mean(y1[treatment == 1])/mean(y0[treatment == 1])
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = data2)
  estimates_RR[i,2] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment*(X1 + X2 + X3), family = binomial(link = 'logit'), data = data2)
  estimates_RR[i,3] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate

  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "nearest", m.order = "largest", distance = 'glm', replace = TRUE)
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,4] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = match_data, weights = weights)
  estimates_RR[i,5] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,6] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "full", estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,7] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = match_data, weights = weights)
  estimates_RR[i,8] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,9] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "cardinality",  ratio = 1, tols = .05, estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,10] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
  poisson_model <- glm(y ~ treatment*(X1 + X2 + X3), family = poisson(link = 'log'), data = match_data, weights = weights)
  estimates_RR[i,11] <- avg_comparisons(poisson_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_RR[i,12] <- avg_comparisons(logist_model, variables = "treatment", comparison = "ratio", newdata = subset(treatment == 1))$estimate
  
}
```

``` r
results <- rbind(
  apply(estimates_RR,2,mean),
  apply(estimates_RR,2,sd))
```

``` r
results_RR_ATT <- t(results)

colnames(results_RR_ATT) <-  c('RD Estimate', 'sd')
rownames(results_RR_ATT) <-  c('True Marginal RR in Treated (mean)', 'Poisson Regression', 'Logistic Regression', 'PSM (Unadjusted)', 'PSM (Poisson Regression)', 'PSM (Logistic Regression)', 'Full Matching (Unadjusted)', 'Full Matching (Poisson Regression)', 'Full Matching (Logistic Regression)', 'Cardinality Matching (Unadjusted))', 'Cardinality Matching (Linear Regression)', 'Cardinality Matching (Logistic Regression)')
results_RR_ATT
```

    ##                                            RD Estimate         sd
    ## True Marginal RR in Treated (mean)            1.293639 0.04733173
    ## Poisson Regression                            1.307349 0.04203552
    ## Logistic Regression                           1.299826 0.03420159
    ## PSM (Unadjusted)                              1.304113 0.07587983
    ## PSM (Poisson Regression)                      1.299850 0.06840758
    ## PSM (Logistic Regression)                     1.306886 0.06397779
    ## Full Matching (Unadjusted)                    1.307819 0.05504463
    ## Full Matching (Poisson Regression)            1.302327 0.04988240
    ## Full Matching (Logistic Regression)           1.306785 0.04575493
    ## Cardinality Matching (Unadjusted))            1.306661 0.04039661
    ## Cardinality Matching (Linear Regression)      1.311008 0.04637228
    ## Cardinality Matching (Logistic Regression)    1.295173 0.03718041

All estimates are very similar to each other.

## Estimating Marginal Odds Ratio under Confounding

Methods for estimating the marginal odds ratio under observed
confounding are very similar to those for estimating risk differences
and risk ratios.

``` r
set.seed(456)

n_sim <- 100
estimates_OD <- matrix(0,n_sim,2)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  estimates_OD[i,1] <- (mean(y[treatment == 1])/mean(1-y[treatment == 1]))/(mean(y[treatment == 0])/mean(1-y[treatment == 0]))
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  estimates_OD[i,2] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
 
}

results <- rbind(
  apply(estimates_OD,2,mean),
  apply(estimates_OD,2,sd))
```

``` r
results_OR <- rbind(
  c((mean(p1)/mean(1-p1))/(mean(p0)/mean(1-p0)), NA),
  c((mean(y1)/mean(1-y1))/(mean(y0)/mean(1-y0)), NA),
  t(results)
)

colnames(results_OR) <-  c('RD Estimate', 'sd')
rownames(results_OR) <-  c('True Marginal OR (mean)', 'True Marginal OR (population)', 'Naive OR', 'Logistic Regression')
results_OR
```

    ##                               RD Estimate        sd
    ## True Marginal OR (mean)          1.885031        NA
    ## True Marginal OR (population)    2.104794        NA
    ## Naive OR                         2.538854 0.2521414
    ## Logistic Regression              2.125749 0.1528363

Again, the unadjusted model provides biased estimates; the correctly
adjusted model is unbiased. We can get the same unbiased result using
weighting.

Let’s denote
``` math
 \hat Y_1 = \frac{\sum_{i=1}^n T_iY_i}{\sum_{i = 1}^n \hat{e}(X_i)T_i} 
```
and
``` math
\hat Y_0 = \frac{\sum_{i=1}^n (1-T_i)Y_i}{\sum_{i = 1}^n \hat{e}(X_i)(1-T_i)}
```

Then the Hájek estimator for odds ratios is
``` math
\text{OR}_\text{Hájek} = \frac{\hat Y_1}{1- \hat Y_1} / \frac{\hat Y_0}{1- \hat Y_0}.
```
We can also compute the augmented IPW estimator. We denote

``` math
\hat Y_1^\text{aug} = \left(\frac{\sum_{i=1}^n T_i(Y_i-\mu_1(X_i, \hat \beta_1))}{\sum_{i = 1}^n \hat{e}(X_i)T_i} + \mu_1(X_i, \hat \beta_1)\right)
```
``` math
\hat Y_0^\text{aug} = \left(\frac{\sum_{i=1}^n (1-T_i)(Y_i - \mu_0(X_i, \hat \beta_0))}{\sum_{i = 1}^n \hat{e}(X_i)(1-T_i)} + \mu_0(X_i, \hat \beta_0)\right)
```
Then,
``` math
\text{OR}_\text{AIPW} =  \frac{\hat Y_1^\text{aug}}{1-\hat Y_1^\text{aug}} / \frac{\hat Y_0^\text{aug}}{1-\hat Y_0^\text{aug}}.
```

``` r
set.seed(456)

n_sim <- 100
estimates_OR <- matrix(0,n_sim,4)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  data2_0 <- data2
  data2_1 <- data2
  data2_0$treatment <- 0
  data2_1$treatment <- 1
  
  
  prop_scores_model_ATE <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "glm", estimand = "ATE")
  
  estimates_OR[i,1] <- ((sum((treatment*y/prop_scores_model_ATE$ps))/sum(((treatment/prop_scores_model_ATE$ps))))/(1-sum((treatment*y/prop_scores_model_ATE$ps))/sum(((treatment/prop_scores_model_ATE$ps))))/((sum(((1-treatment)*y/(1-prop_scores_model_ATE$ps)))/sum((((1-treatment)/(1-prop_scores_model_ATE$ps)))))/(1-(sum(((1-treatment)*y/(1-prop_scores_model_ATE$ps)))/sum((((1-treatment)/(1-prop_scores_model_ATE$ps))))))))
  
  
  logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  
  estimates_OR[i,2] <- avg_comparisons(logist_unadj_model, comparison = "lnoravg", transform = "exp")$estimate
  
  logist_model <-glm(y ~ treatment*(X1 + X2 + X3), family = binomial(link = 'logit'), data = data2)
  pred0 <- predict(logist_model, newdata = data2_0, type = 'response')
  pred1 <- predict(logist_model, newdata = data2_1, type = 'response')
  
  
   estimates_OR[i,3] <- (mean(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps)/(1 - mean(pred1  + treatment*(y-pred1)/prop_scores_model_ATE$ps)))/
     (mean(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))/(1-mean(pred0 + (1-treatment)*(y-pred0)/(1-prop_scores_model_ATE$ps))))
   
   
   logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  
  estimates_OR[i,4] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate

}  
  
results <- rbind(
  apply(estimates_OR,2,mean),
  apply(estimates_OR,2,sd))
```

``` r
results_OR2 <- rbind(
  t(results)
)

colnames(results_OR2) <-  c('RD Estimate', 'sd')
rownames(results_OR2) <-  c('IPW (Hájek)', 'IPW (Unadj. Logistic Regression)',  'AIPW (Logistic Regression)', 'IPWRA (Logistic Regression)')
results_OR2
```

    ##                                  RD Estimate        sd
    ## IPW (Hájek)                         2.129634 0.1637508
    ## IPW (Unadj. Logistic Regression)    2.129634 0.1637508
    ## AIPW (Logistic Regression)          2.130402 0.1561759
    ## IPWRA (Logistic Regression)         2.130083 0.1562991

We notice that again the Hájek estimator is equivalent to the unadjusted
regression. Lastly, let’s check weighting and matching.

``` r
set.seed(758)

n_sim <- 100
estimates_OR <- matrix(0,n_sim,6)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  prop_scores_model_cbps <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "cbps", estimand = "ATE", moments = 2, over = FALSE)
  
  logist_unadj_model <- glm(y ~ treatment, data = data2, family = binomial(link = 'logit'), weights = prop_scores_model_cbps$weights)
  
  estimates_OR[i,1] <- avg_comparisons(logist_unadj_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
  

  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_cbps$weights)
  estimates_OR[i,2] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
  
   prop_scores_model_ebal <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "ebal", estimand = "ATE", moments = 2)
  
   logist_unadj_model <- glm(y ~ treatment, data = data2, family = binomial(link = 'logit'), weights = prop_scores_model_ebal$weights)
   
  estimates_OR[i,3] <- avg_comparisons(logist_unadj_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
  
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ebal$weights)
  
  estimates_OR[i,4] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
  
  prop_scores_model_eng <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "energy", estimand = "ATE")
  
   logist_unadj_model <- glm(y ~ treatment, data = data2, family = binomial(link = 'logit'), weights = prop_scores_model_eng$weights)
  estimates_OR[i,5] <- avg_comparisons(logist_unadj_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
  
  logist_model <-glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_eng$weights)
  
  estimates_OR[i,6] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp")$estimate
}  
  
results <- rbind(
  apply(estimates_OR,2,mean),
  apply(estimates_OR,2,sd))
```

``` r
results_OR3 <- rbind(
  t(results)
)

colnames(results_OR3) <-  c('RD Estimate', 'sd')
rownames(results_OR3) <-  c('IPW CBPS (Unadjusted)',  'IPW CBPS (Logistic Regression)', 'Entropy Balancing (Unadjusted)',  'Entropy Balancing (Logistic Regression)', 'Energy Balancing (Unadjusted)',  'Energy Balancing (Logistic Regression)')
results_OR3
```

    ##                                         RD Estimate        sd
    ## IPW CBPS (Unadjusted)                      2.105410 0.1717776
    ## IPW CBPS (Logistic Regression)             2.112206 0.1744043
    ## Entropy Balancing (Unadjusted)             2.099593 0.1712276
    ## Entropy Balancing (Logistic Regression)    2.111484 0.1745871
    ## Energy Balancing (Unadjusted)              2.109346 0.2207576
    ## Energy Balancing (Logistic Regression)     2.107656 0.2207476

``` r
set.seed(456)

n_sim <- 100
estimates_OR <- matrix(0,n_sim,8)

for (i in 1:n_sim){

  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  estimates_OR[i,1] <- (mean(y1[treatment == 1])/(1-mean(y1[treatment == 1])))/(mean(y0[treatment == 1])/(1-mean(y0[treatment == 1])))
    
  logist_model <- glm(y ~ treatment*(X1 + X2 + X3), family = binomial(link = 'logit'), data = data2)
  estimates_OR[i,2] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate

  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "nearest", m.order = "largest", distance = 'glm', replace = TRUE)
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,3] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate
  

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,4] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "full", estimand = "ATT")
  match_data <- match_data(match_pp)
  
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,5] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate
  
  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,6] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "cardinality",  ratio = 1, tols = .05, estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,7] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_OR[i,8] <- avg_comparisons(logist_model, variables = "treatment", comparison = "lnoravg", transform = "exp", newdata = subset(treatment == 1))$estimate
}
```

``` r
results <- rbind(
  apply(estimates_OR,2,mean),
  apply(estimates_OR,2,sd))
```

``` r
results_OR_ATT <- t(results)

colnames(results_OR_ATT) <-  c('RD Estimate', 'sd')
rownames(results_OR_ATT) <-  c('True Marginal OR in Treated (mean)', 'Logistic Regression', 'PSM (Unadjusted)', 'PSM (Logistic Regression)', 'Full Matching (Unadjusted)', 'Full Matching (Logistic Regression)', 'Cardinality Matching (Unadjusted)', 'Cardinality Matching (Logistic Regression)')
results_OR_ATT
```

    ##                                            RD Estimate        sd
    ## True Marginal OR in Treated (mean)            2.132379 0.2337930
    ## Logistic Regression                           2.150869 0.1596171
    ## PSM (Unadjusted)                              2.165819 0.2936468
    ## PSM (Logistic Regression)                     2.177852 0.2613129
    ## Full Matching (Unadjusted)                    2.180092 0.2178192
    ## Full Matching (Logistic Regression)           2.176765 0.1911722
    ## Cardinality Matching (Unadjusted)             2.177747 0.1861820
    ## Cardinality Matching (Logistic Regression)    2.133199 0.1710064

## Conditional Odds Ratio under Confounding

To conclude, let’s also check the conditional odds ratio.

``` r
set.seed(456)

n_sim <- 100
estimates_COR <- matrix(0,n_sim,2)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  logist_model_unadj <- glm(y ~ treatment, family = binomial(link = 'logit'), data = data2)
  estimates_COR[i,1] <- logist_model_unadj$coefficients[2]
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  estimates_COR[i,2] <- logist_model$coefficients[2]
 
}

results <- rbind(
  apply(estimates_COR,2,mean),
  apply(estimates_COR,2,sd))
```

``` r
results_COR <- rbind(
  c(1, NA),
  t(results)
)

colnames(results_COR) <-  c('OR Estimate (logarithm)', 'sd')
rownames(results_COR) <-  c('True COR', 'Naive (unadjusted) COR', 'COR (logistic regression)')
results_COR
```

    ##                           OR Estimate (logarithm)         sd
    ## True COR                                 1.000000         NA
    ## Naive (unadjusted) COR                   0.926866 0.09889487
    ## COR (logistic regression)                1.207666 0.11031204

We notice that both unadjusted and adjusted estimates appear biased.
However, the adjusted estimate is biased merely due to finite sampling.

``` r
set.seed(123)
n_pop <- 2000

X1 <- runif(n_pop, -5, 5)
X2 <- rnorm(n_pop, 2.5, 1)
X3 <- (runif(n_pop, 0, 1)>0.75) + 0

p0 <- plogis(-2 + 0.5*X1 + 0.75*X2 + 2*X3)
p1 <- plogis(-2 + 0.5*X1 + 0.75*X2 + 2*X3 + rep(1,n_pop))

y0 <- (runif(n_pop,0,1) > (1-p0)) + 0
y1 <- (runif(n_pop,0,1) > (1-p1)) + 0
```

``` r
set.seed(456)

n_sim <- 100
estimates_COR <- matrix(0,n_sim,2)

for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  logist_model_unadj <- glm(y ~ treatment, family = binomial(link = 'logit'), data = data2)
  estimates_COR[i,1] <- logist_model_unadj$coefficients[2]
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2)
  estimates_COR[i,2] <- logist_model$coefficients[2]
 
}

results <- rbind(
  apply(estimates_COR,2,mean),
  apply(estimates_COR,2,sd))
```

``` r
results_COR <- rbind(
  c(1, NA),
  t(results)
)

colnames(results_COR) <-  c('OR Estimate (logarithm)', 'sd')
rownames(results_COR) <-  c('True COR', 'Naive (unadjusted) COR', 'COR (logistic regression)')
results_COR
```

    ##                           OR Estimate (logarithm)         sd
    ## True COR                                 1.000000         NA
    ## Naive (unadjusted) COR                   0.819566 0.07691336
    ## COR (logistic regression)                1.034492 0.09530593

We see that doubling the number of observations makes the logistic
regression estimate more precise, whereas the unadjusted estimate
remains biased. However, we can clearly see how much more difficult it
is to estimate the conditional effect than the marginal effect: we
needed 2000 samples to obtain a reasonably accurate estimate, and this
is a simple problem with the treatment and 4 covariates, with no
interactions. With 1000, the estimate provided by the correctly
specified model was more biased than the unadjusted estimate (the effect
was overestimated, which is typical in finite samples (Nemes et al.
2009)).

Let’s have a look at weighted models.

``` r
set.seed(456)

n_sim <- 100
estimates_COR <- matrix(0,n_sim,6)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  prop_scores_model_ATE <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "glm", estimand = "ATE")
  
  logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_COR[i,1] <- logist_unadj_model$coefficients[2]
  
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ATE$weights)
  estimates_COR[i,2] <- logist_model$coefficients[2]
  
  prop_scores_model_cbps <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "cbps", estimand = "ATE", moments = 2, over = FALSE)
  
  logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_cbps$weights)
  estimates_COR[i,3] <- logist_unadj_model$coefficients[2]
  
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_cbps$weights)
  estimates_COR[i,4] <- logist_model$coefficients[2]
  
  
  prop_scores_model_ebal <- weightit(treatment ~ X1 + X2 + X3, data = data2, method = "ebal", estimand = "ATE", moments = 2)
  
    logist_unadj_model <-glm(y ~ treatment, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ebal$weights)
  estimates_COR[i,5] <- logist_unadj_model$coefficients[2]
  
  
  logist_model <- glm(y ~ treatment + X1 + X2 + X3, family = binomial(link = 'logit'), data = data2, weights = prop_scores_model_ebal$weights)
  estimates_COR[i,6] <- logist_model$coefficients[2]
  
  
}  
  
results <- rbind(
  apply(estimates_COR,2,mean),
  apply(estimates_COR,2,sd))
```

``` r
results_COR <- rbind(
  t(results)
)

colnames(results_COR) <-  c('OR Estimate (logarithm)', 'sd')
rownames(results_COR) <-  c('IPW (unadjusted)', 'IPWRA', 'CBPS (unadjusted)', 'CBPS', 'Energy Balancing (unadjusted)', 'Energy Balancing')
results_COR
```

    ##                               OR Estimate (logarithm)         sd
    ## IPW (unadjusted)                            0.6497808 0.06269997
    ## IPWRA                                       1.0343969 0.09857865
    ## CBPS (unadjusted)                           0.6507219 0.06238806
    ## CBPS                                        1.0348028 0.09880935
    ## Energy Balancing (unadjusted)               0.6479661 0.06206189
    ## Energy Balancing                            1.0330200 0.09848872

We make a very important observation: weighting balances the population,
but alone is not enough to obtain an unbiased conditional odds ratio due
to noncollapsibility. The same is true for matching.

``` r
set.seed(456)

n_sim <- 100
estimates_COR <- matrix(0,n_sim,6)


for (i in 1:n_sim){
  
  treatment <- (runif(n_pop,0,1) < plogis(-0.3 - 0.5*X3 + 0.1*X1)) + 0
  
  y <- y0
  y[treatment == 1] <- y1[treatment == 1]
  data2 <- data.frame(y = y, treatment = treatment, X1 = X1, X2 = X2, X3 = X3, p0 = p0, p1 = p1, y0 = y0, y1 = y1)
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "nearest", m.order = "largest", distance = 'glm', replace = TRUE)
  match_data <- match_data(match_pp)
  
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,1] <- logist_model$coefficients[2]
  
  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,2] <- logist_model$coefficients[2]
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "full", estimand = "ATT")
  match_data <- match_data(match_pp)
  
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,3] <- logist_model$coefficients[2]
  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,4] <- logist_model$coefficients[2]
  
  
  match_pp <- matchit(treatment ~ X1 + X2 + X3, data = data2, method = "cardinality",  ratio = 1, tols = .05, estimand = "ATT")
  match_data <- match_data(match_pp)
  
  logist_model <- glm(y ~ treatment, family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,5] <- logist_model$coefficients[2]

  logist_model <- glm(y ~ treatment + (X1 + X2 + X3), family = binomial(link = 'logit'), data = match_data, weights = weights)
  estimates_COR[i,6] <- logist_model$coefficients[2]
}
```

``` r
results <- rbind(
  apply(estimates_COR,2,mean),
  apply(estimates_COR,2,sd))
```

``` r
results_OR_ATT <- t(results)

colnames(results_OR_ATT) <-  c('RD Estimate', 'sd')
rownames(results_OR_ATT) <-  c('PSM (unadj.)', 'PSM (Logistic Regression)', 'Full Matching (unadj.)', 'Full Matching (Logistic Regression)','Cardinality Matching (unadj.)',  'Cardinality Matching (Logistic Regression)')
results_OR_ATT
```

    ##                                            RD Estimate         sd
    ## PSM (unadj.)                                 0.6668621 0.08764805
    ## PSM (Logistic Regression)                    1.0422284 0.12886368
    ## Full Matching (unadj.)                       0.6571451 0.07397881
    ## Full Matching (Logistic Regression)          1.0359805 0.11150806
    ## Cardinality Matching (unadj.)                0.6604632 0.06414917
    ## Cardinality Matching (Logistic Regression)   0.9948984 0.09736231

This implies that matching and weighting are less useful when estimating
conditional effects for binary data, since we have to use a
well-specified model even for data without confounding, as we discussed
in the previous part.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-boughdiri2024quantifying" class="csl-entry">

Boughdiri, Ahmed, Julie Josse, and Erwan Scornet. 2024. “Quantifying
Treatment Effects: Estimating Risk Ratios in Causal Inference.” *arXiv
Preprint arXiv:2410.12333*.

</div>

<div id="ref-nemes2009bias" class="csl-entry">

Nemes, Szilard, Junmei Miao Jonasson, Anna Genell, and Gunnar Steineck.
2009. “Bias in Odds Ratios by Logistic Regression Modelling and Sample
Size.” *BMC Medical Research Methodology* 9 (1): 56.

</div>

</div>
