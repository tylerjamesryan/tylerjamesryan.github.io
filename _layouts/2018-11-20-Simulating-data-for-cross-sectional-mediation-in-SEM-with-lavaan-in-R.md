Simulating data for cross-sectional mediation in SEM with lavaan in R
================
2018-11-20

We'll need the `lavaan` package for the actual analysis and the `ggplot2` package for visualizing the effect.

``` r
library(lavaan)
```

    ## This is lavaan 0.6-2

    ## lavaan is BETA software! Please report any bugs.

``` r
library(ggplot2)
```

Generate Data
=============

``` r
N <- 10000

X <- rnorm(N, 0, 1)
Med <- .7 * X  + rnorm(N, 0, 1)
Y <- .3 * X + .3 * Med + rnorm(N, 0, 1)
dat <- cbind.data.frame(X, Med, Y)
```

We then create our observed variables with loadings of .8 and residual errors of ~ .36. The observed variables will thus have a total variance of ~ 1.0.

``` r
dat$x1 <- .8 * X +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$x2 <- .8 * X +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$x3 <- .8 * X +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$m1 <- .8 * Med + (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$m2 <- .8 * Med + (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$m3 <- .8 * Med + (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$y1 <- .8 * Y +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$y2 <- .8 * Y +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))
dat$y3 <- .8 * Y +   (sqrt(1 - .8^2) * rnorm(N, 0, 1))

model <- '

            X =~ x1 + x2 + x3
            Med =~ m1 + m2 + m3
            Y =~ y1 + y2 + y3

            Med ~ a*X
            Y ~ c*X
            Y ~ b*Med

            ab := a*b

            abc := ab + c
          
'


mod <- sem(model = model, data = dat, std.lv = T)

summary(mod, fit.measures = T)
```

    ## lavaan 0.6-2 ended normally after 25 iterations
    ## 
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         21
    ## 
    ##   Number of observations                         10000
    ## 
    ##   Estimator                                         ML
    ##   Model Fit Test Statistic                      23.556
    ##   Degrees of freedom                                24
    ##   P-value (Chi-square)                           0.487
    ## 
    ## Model test baseline model:
    ## 
    ##   Minimum Function Test Statistic            49793.578
    ##   Degrees of freedom                                36
    ##   P-value                                        0.000
    ## 
    ## User model versus baseline model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000
    ##   Tucker-Lewis Index (TLI)                       1.000
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)             -109845.538
    ##   Loglikelihood unrestricted model (H1)     -109833.760
    ## 
    ##   Number of free parameters                         21
    ##   Akaike (AIC)                              219733.076
    ##   Bayesian (BIC)                            219884.493
    ##   Sample-size adjusted Bayesian (BIC)       219817.758
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000
    ##   90 Percent Confidence Interval          0.000  0.008
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.004
    ## 
    ## Parameter Estimates:
    ## 
    ##   Information                                 Expected
    ##   Information saturated (h1) model          Structured
    ##   Standard Errors                             Standard
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   X =~                                                
    ##     x1                0.801    0.009   88.555    0.000
    ##     x2                0.798    0.009   88.905    0.000
    ##     x3                0.794    0.009   89.095    0.000
    ##   Med =~                                              
    ##     m1                0.787    0.009   91.586    0.000
    ##     m2                0.796    0.009   92.249    0.000
    ##     m3                0.777    0.009   91.185    0.000
    ##   Y =~                                                
    ##     y1                0.795    0.008   93.850    0.000
    ##     y2                0.799    0.008   94.663    0.000
    ##     y3                0.788    0.009   92.591    0.000
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Med ~                                               
    ##     X          (a)    0.722    0.015   46.640    0.000
    ##   Y ~                                                 
    ##     X          (c)    0.322    0.016   19.640    0.000
    ##     Med        (b)    0.296    0.013   22.500    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .x1                0.372    0.008   48.161    0.000
    ##    .x2                0.363    0.008   47.782    0.000
    ##    .x3                0.357    0.007   47.575    0.000
    ##    .m1                0.366    0.008   47.089    0.000
    ##    .m2                0.357    0.008   45.962    0.000
    ##    .m3                0.366    0.008   47.740    0.000
    ##    .y1                0.357    0.008   46.023    0.000
    ##    .y2                0.343    0.008   44.762    0.000
    ##    .y3                0.377    0.008   47.844    0.000
    ##     X                 1.000                           
    ##    .Med               1.000                           
    ##    .Y                 1.000                           
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     ab                0.214    0.010   21.512    0.000
    ##     abc               0.535    0.014   38.133    0.000
