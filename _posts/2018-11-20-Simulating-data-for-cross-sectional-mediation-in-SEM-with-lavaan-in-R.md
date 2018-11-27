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

    ## lavaan 0.6-2 ended normally after 27 iterations
    ## 
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         21
    ## 
    ##   Number of observations                         10000
    ## 
    ##   Estimator                                         ML
    ##   Model Fit Test Statistic                      17.937
    ##   Degrees of freedom                                24
    ##   P-value (Chi-square)                           0.806
    ## 
    ## Model test baseline model:
    ## 
    ##   Minimum Function Test Statistic            51040.278
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
    ##   Loglikelihood user model (H0)             -109891.122
    ##   Loglikelihood unrestricted model (H1)     -109882.154
    ## 
    ##   Number of free parameters                         21
    ##   Akaike (AIC)                              219824.245
    ##   Bayesian (BIC)                            219975.662
    ##   Sample-size adjusted Bayesian (BIC)       219908.927
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000
    ##   90 Percent Confidence Interval          0.000  0.005
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.003
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
    ##     x1                0.804    0.009   89.692    0.000
    ##     x2                0.818    0.009   90.703    0.000
    ##     x3                0.813    0.009   90.620    0.000
    ##   Med =~                                              
    ##     m1                0.814    0.009   94.286    0.000
    ##     m2                0.812    0.009   93.715    0.000
    ##     m3                0.811    0.009   93.718    0.000
    ##   Y =~                                                
    ##     y1                0.791    0.008   94.033    0.000
    ##     y2                0.794    0.008   94.400    0.000
    ##     y3                0.804    0.009   94.406    0.000
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Med ~                                               
    ##     X          (a)    0.683    0.015   45.733    0.000
    ##   Y ~                                                 
    ##     X          (c)    0.299    0.016   18.899    0.000
    ##     Med        (b)    0.335    0.013   25.454    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .x1                0.363    0.008   48.165    0.000
    ##    .x2                0.358    0.008   47.044    0.000
    ##    .x3                0.355    0.008   47.137    0.000
    ##    .m1                0.351    0.008   46.194    0.000
    ##    .m2                0.362    0.008   47.125    0.000
    ##    .m3                0.361    0.008   47.119    0.000
    ##    .y1                0.358    0.008   46.628    0.000
    ##    .y2                0.353    0.008   46.081    0.000
    ##    .y3                0.362    0.008   46.072    0.000
    ##     X                 1.000                           
    ##    .Med               1.000                           
    ##    .Y                 1.000                           
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     ab                0.228    0.010   23.711    0.000
    ##     abc               0.527    0.014   37.634    0.000
