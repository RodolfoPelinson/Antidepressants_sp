Predictors of total concentration of antidepressants in streamwater
(excluding observations)
================
Rodolfo Pelinson
2024-11-15

Packages:

``` r
library(DHARMa)
library(AICcmodavg)
library(corrplot)
library(MuMIn)
library(vegan)
library(car)
library(lme4)
library(glmmTMB)
library(dplyr)
```

A few necessary functions to make the code more straightforward:

``` r
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/plot_models_one_var.R")
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/plot_models_two_var.R")
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/confidence_interval.R")
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/ajeitando_planilhas.R")
```

#### Excluding observations from which there was more than 10% variation in urbanization from 2010 to 2021.

Excluding observations from which there was more than 10% variation in
urbanization from 2010 to 2021.

``` r
diff_watersheds_10_urb
```

    ## [1] "b309" "b595" "b611" "b617"

``` r
new_predictors <- predictors
predictors_sem_10_urb <- predictors[-match(diff_watersheds_10_urb, rownames(new_predictors)),]
```

Performing the model selection procedure. Note that the name of the
variables are in portuguese.

Habitantes mean inhabitants casas_sem_saneamento means houses without
sanitation renda means income.

First we make a matrix to check which pairs of variables are correlated
(more than 50%) and therefore should not be together in the same model.

``` r
pred_mod_sel <- select(predictors_sem_10_urb,  casas_sem_saneamento, renda, habitantes)

#Standardizing response variables
Preditoras_stand <- decostand(pred_mod_sel, method = "stand")

smat <- abs(cor(Preditoras_stand)) <= .5
smat[!lower.tri(smat)] <- NA

colnames(smat) <- c(paste("cond(",colnames(smat)[1], ")", sep = ""),
                    paste("cond(",colnames(smat)[2], ")", sep = ""),
                    paste("cond(",colnames(smat)[3], ")", sep = ""))

rownames(smat) <- colnames(smat)

cor(Preditoras_stand)
```

    ##                      casas_sem_saneamento      renda habitantes
    ## casas_sem_saneamento            1.0000000 -0.3537566  0.6095417
    ## renda                          -0.3537566  1.0000000 -0.2672892
    ## habitantes                      0.6095417 -0.2672892  1.0000000

``` r
smat
```

    ##                            cond(casas_sem_saneamento) cond(renda)
    ## cond(casas_sem_saneamento)                         NA          NA
    ## cond(renda)                                      TRUE          NA
    ## cond(habitantes)                                FALSE        TRUE
    ##                            cond(habitantes)
    ## cond(casas_sem_saneamento)               NA
    ## cond(renda)                              NA
    ## cond(habitantes)                         NA

Now performing the model selection.

``` r
all <- responses$all[-match(diff_watersheds_10_urb, rownames(responses))]

Preditoras_stand$all <- all

mod_global_all <- glmmTMB(all ~ 
                            #area +
                            casas_sem_saneamento +
                            habitantes +
                            renda + 
                            casas_sem_saneamento:renda+
                            #casas_sem_saneamento:renda:area+
                            habitantes:renda,
                            #habitantes:renda:area+
                            #casas_sem_saneamento:area+
                            #habitantes:area+
                            #renda:area,
                          data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")


dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##    K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 3  4 497.08       0.00   0.48   0.48 -244.07
    ## 14 6 497.75       0.67   0.34   0.82 -241.83
    ## 23 6 498.99       1.90   0.18   1.00 -242.44

``` r
dredge_modelos_all
```

    ## Global model call: glmmTMB(formula = all ~ casas_sem_saneamento + habitantes + renda + 
    ##     casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt) cnd(rnd)
    ## 3       5.202          +                    0.3597         
    ## 14      5.846          +          1.13700           0.69490
    ## 23      5.304          +                    0.5228  0.05248
    ## 1       5.273          +                                   
    ## 7       5.199          +                    0.3442 -0.07811
    ## 5       5.260          +                           -0.17900
    ## 2       5.270          +         -0.08209                  
    ## 6       5.251          +         -0.15060          -0.22970
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 3                                     4 -244.066 497.1  0.00  0.330
    ## 14                2.123               6 -241.826 497.8  0.67  0.237
    ## 23                            0.6248  6 -242.443 499.0  1.90  0.128
    ## 1                                     3 -246.406 499.4  2.29  0.105
    ## 7                                     5 -243.996 499.5  2.37  0.101
    ## 5                                     4 -246.022 501.0  3.91  0.047
    ## 2                                     4 -246.314 501.6  4.50  0.035
    ## 6                                     5 -245.739 502.9  5.86  0.018
    ## Models ranked by AICc(x)

Selected models include houses without sanitaiton in interaction with
income OR just inhabitants OR inhabitants in interaction with income.

#### Excluding observations from which there was more than 5% variation in urbanization from 2010 to 2021.

Excluding observations from which there was more than 10% variation in
urbanization from 2010 to 2021.

``` r
diff_watersheds_5_urb
```

    ##  [1] "b030" "b035" "b039" "b040" "b067" "b309" "b595" "b611" "b617" "b636"
    ## [11] "b637"

``` r
new_predictors <- predictors
predictors_sem_10_urb <- predictors[-match(diff_watersheds_5_urb, rownames(new_predictors)),]
```

Performing the model selection procedure. Note that the name of the
variables are in portuguese.

Habitantes mean inhabitants casas_sem_saneamento means houses without
sanitation renda means income.

First we make a matrix to check which pairs of variables are correlated
(more than 50%) and therefore should not be together in the same model.

``` r
pred_mod_sel <- select(predictors_sem_10_urb,  casas_sem_saneamento, renda, habitantes)

#Standardizing response variables
Preditoras_stand <- decostand(pred_mod_sel, method = "stand")

smat <- abs(cor(Preditoras_stand)) <= .5
smat[!lower.tri(smat)] <- NA

colnames(smat) <- c(paste("cond(",colnames(smat)[1], ")", sep = ""),
                    paste("cond(",colnames(smat)[2], ")", sep = ""),
                    paste("cond(",colnames(smat)[3], ")", sep = ""))

rownames(smat) <- colnames(smat)

cor(Preditoras_stand)
```

    ##                      casas_sem_saneamento      renda habitantes
    ## casas_sem_saneamento            1.0000000 -0.3406279  0.6207516
    ## renda                          -0.3406279  1.0000000 -0.2814300
    ## habitantes                      0.6207516 -0.2814300  1.0000000

``` r
smat
```

    ##                            cond(casas_sem_saneamento) cond(renda)
    ## cond(casas_sem_saneamento)                         NA          NA
    ## cond(renda)                                      TRUE          NA
    ## cond(habitantes)                                FALSE        TRUE
    ##                            cond(habitantes)
    ## cond(casas_sem_saneamento)               NA
    ## cond(renda)                              NA
    ## cond(habitantes)                         NA

Now performing the model selection.

``` r
all <- responses$all[-match(diff_watersheds_5_urb, rownames(responses))]

Preditoras_stand$all <- all

mod_global_all <- glmmTMB(all ~ 
                            #area +
                            casas_sem_saneamento +
                            habitantes +
                            renda + 
                            casas_sem_saneamento:renda+
                            #casas_sem_saneamento:renda:area+
                            habitantes:renda,
                            #habitantes:renda:area+
                            #casas_sem_saneamento:area+
                            #habitantes:area+
                            #renda:area,
                          data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")


dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##    K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 3  4 400.56       0.00   0.54   0.54 -195.71
    ## 1  3 402.18       1.63   0.24   0.79 -197.76
    ## 23 6 402.43       1.87   0.21   1.00 -193.94

``` r
dredge_modelos_all
```

    ## Global model call: glmmTMB(formula = all ~ casas_sem_saneamento + habitantes + renda + 
    ##     casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt)  cnd(rnd)
    ## 3       5.204          +                    0.3440          
    ## 1       5.270          +                                    
    ## 23      5.311          +                    0.5682  0.005717
    ## 7       5.195          +                    0.3141 -0.138500
    ## 14      5.820          +          1.29300           0.554800
    ## 5       5.246          +                           -0.238000
    ## 2       5.269          +         -0.05136                   
    ## 6       5.239          +         -0.13370          -0.282100
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 3                                     4 -195.707 400.6  0.00  0.343
    ## 1                                     3 -197.758 402.2  1.63  0.152
    ## 23                            0.6631  6 -193.940 402.4  1.87  0.135
    ## 7                                     5 -195.504 402.8  2.21  0.113
    ## 14                2.094               6 -194.246 403.0  2.48  0.099
    ## 5                                     4 -197.125 403.4  2.84  0.083
    ## 2                                     4 -197.725 404.6  4.04  0.046
    ## 6                                     5 -196.914 405.6  5.03  0.028
    ## Models ranked by AICc(x)

Selected models include just inhabitants OR inhabitants in interaction
with income OR the model without effect of any predictors. The inclusion
of the this last model probably happened because of the significant
reduction in sample size (11 watersheds were removed).

#### Excluding observations from which there was more than 10% variation in subnormal constructions from 2010 to 2021.

Excluding observations from which there was more than 10% variation in
urbanization from 2010 to 2021.

``` r
diff_watersheds_10_AGSB
```

    ## [1] "b035" "b040" "b066" "b202" "b203" "b205" "b309" "b595" "b611"

``` r
new_predictors <- predictors
predictors_sem_10_urb <- predictors[-match(diff_watersheds_10_AGSB, rownames(new_predictors)),]
```

Performing the model selection procedure. Note that the name of the
variables are in portuguese.

Habitantes mean inhabitants casas_sem_saneamento means houses without
sanitation renda means income.

First we make a matrix to check which pairs of variables are correlated
(more than 50%) and therefore should not be together in the same model.

``` r
pred_mod_sel <- select(predictors_sem_10_urb,  casas_sem_saneamento, renda, habitantes)

#Standardizing response variables
Preditoras_stand <- decostand(pred_mod_sel, method = "stand")

smat <- abs(cor(Preditoras_stand)) <= .5
smat[!lower.tri(smat)] <- NA

colnames(smat) <- c(paste("cond(",colnames(smat)[1], ")", sep = ""),
                    paste("cond(",colnames(smat)[2], ")", sep = ""),
                    paste("cond(",colnames(smat)[3], ")", sep = ""))

rownames(smat) <- colnames(smat)

cor(Preditoras_stand)
```

    ##                      casas_sem_saneamento      renda habitantes
    ## casas_sem_saneamento            1.0000000 -0.3346052  0.6089342
    ## renda                          -0.3346052  1.0000000 -0.2233528
    ## habitantes                      0.6089342 -0.2233528  1.0000000

``` r
smat
```

    ##                            cond(casas_sem_saneamento) cond(renda)
    ## cond(casas_sem_saneamento)                         NA          NA
    ## cond(renda)                                      TRUE          NA
    ## cond(habitantes)                                FALSE        TRUE
    ##                            cond(habitantes)
    ## cond(casas_sem_saneamento)               NA
    ## cond(renda)                              NA
    ## cond(habitantes)                         NA

Now performing the model selection.

``` r
all <- responses$all[-match(diff_watersheds_10_AGSB, rownames(responses))]

Preditoras_stand$all <- all

mod_global_all <- glmmTMB(all ~ 
                            #area +
                            casas_sem_saneamento +
                            habitantes +
                            renda + 
                            casas_sem_saneamento:renda+
                            #casas_sem_saneamento:renda:area+
                            habitantes:renda,
                            #habitantes:renda:area+
                            #casas_sem_saneamento:area+
                            #habitantes:area+
                            #renda:area,
                          data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")


dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##    K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 3  4 429.46       0.00   0.53   0.53 -210.19
    ## 14 6 430.83       1.37   0.27   0.79 -208.22
    ## 1  3 431.32       1.85   0.21   1.00 -212.34

``` r
dredge_modelos_all
```

    ## Global model call: glmmTMB(formula = all ~ casas_sem_saneamento + habitantes + renda + 
    ##     casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt) cnd(rnd)
    ## 3       5.208          +                    0.3516         
    ## 14      5.825          +         1.304000           0.60800
    ## 1       5.279          +                                   
    ## 7       5.201          +                    0.3314 -0.12480
    ## 23      5.276          +                    0.5321 -0.03115
    ## 5       5.261          +                           -0.20630
    ## 2       5.279          +        -0.005801                  
    ## 6       5.259          +        -0.071610          -0.23010
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 3                                     4 -210.190 429.5  0.00  0.347
    ## 14                2.188               6 -208.217 430.8  1.37  0.175
    ## 1                                     3 -212.342 431.3  1.85  0.137
    ## 7                                     5 -210.023 431.7  2.25  0.112
    ## 23                            0.6094  6 -208.719 431.8  2.38  0.106
    ## 5                                     4 -211.878 432.8  3.38  0.064
    ## 2                                     4 -212.341 433.8  4.30  0.040
    ## 6                                     5 -211.818 435.3  5.84  0.019
    ## Models ranked by AICc(x)

Selected models include just inhabitants OR houses without sanitation in
interaction with income OR the model without effect of any predictors.
The inclusion of the this last model probably happened because of the
significant reduction in sample size (9 watersheds were removed).

#### Excluding observations from which there was more than 5% variation in subnormal constructions from 2010 to 2021.

Excluding observations from which there was more than 10% variation in
urbanization from 2010 to 2021.

``` r
diff_watersheds_5_AGSB
```

    ##  [1] "b035" "b039" "b040" "b066" "b202" "b203" "b205" "b309" "b573" "b579"
    ## [11] "b595" "b611" "b618" "b637"

``` r
new_predictors <- predictors
predictors_sem_10_urb <- predictors[-match(diff_watersheds_5_AGSB, rownames(new_predictors)),]
```

Performing the model selection procedure. Note that the name of the
variables are in portuguese.

Habitantes mean inhabitants casas_sem_saneamento means houses without
sanitation renda means income.

First we make a matrix to check which pairs of variables are correlated
(more than 50%) and therefore should not be together in the same model.

``` r
pred_mod_sel <- select(predictors_sem_10_urb,  casas_sem_saneamento, renda, habitantes)

#Standardizing response variables
Preditoras_stand <- decostand(pred_mod_sel, method = "stand")

smat <- abs(cor(Preditoras_stand)) <= .5
smat[!lower.tri(smat)] <- NA

colnames(smat) <- c(paste("cond(",colnames(smat)[1], ")", sep = ""),
                    paste("cond(",colnames(smat)[2], ")", sep = ""),
                    paste("cond(",colnames(smat)[3], ")", sep = ""))

rownames(smat) <- colnames(smat)

cor(Preditoras_stand)
```

    ##                      casas_sem_saneamento      renda habitantes
    ## casas_sem_saneamento            1.0000000 -0.2883068  0.7125644
    ## renda                          -0.2883068  1.0000000 -0.1962107
    ## habitantes                      0.7125644 -0.1962107  1.0000000

``` r
smat
```

    ##                            cond(casas_sem_saneamento) cond(renda)
    ## cond(casas_sem_saneamento)                         NA          NA
    ## cond(renda)                                      TRUE          NA
    ## cond(habitantes)                                FALSE        TRUE
    ##                            cond(habitantes)
    ## cond(casas_sem_saneamento)               NA
    ## cond(renda)                              NA
    ## cond(habitantes)                         NA

Now performing the model selection.

``` r
all <- responses$all[-match(diff_watersheds_5_AGSB, rownames(responses))]

Preditoras_stand$all <- all

mod_global_all <- glmmTMB(all ~ 
                            #area +
                            casas_sem_saneamento +
                            habitantes +
                            renda + 
                            casas_sem_saneamento:renda+
                            #casas_sem_saneamento:renda:area+
                            habitantes:renda,
                            #habitantes:renda:area+
                            #casas_sem_saneamento:area+
                            #habitantes:area+
                            #renda:area,
                          data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")


dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##    K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 3  4 360.07       0.00   0.54   0.54 -175.41
    ## 23 6 361.79       1.72   0.23   0.78 -173.50
    ## 1  3 361.85       1.77   0.22   1.00 -177.56

``` r
dredge_modelos_all
```

    ## Global model call: glmmTMB(formula = all ~ casas_sem_saneamento + habitantes + renda + 
    ##     casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt) cnd(rnd)
    ## 3       5.191          +                    0.3512         
    ## 23      5.200          +                    0.5930  -0.1905
    ## 1       5.270          +                                   
    ## 7       5.178          +                    0.3253  -0.1747
    ## 5       5.243          +                            -0.2508
    ## 14      5.641          +          1.18500            0.2639
    ## 2       5.268          +         -0.06608                  
    ## 6       5.235          +         -0.14270           -0.2874
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 3                                     4 -175.410 360.1  0.00  0.360
    ## 23                            0.5527  6 -173.495 361.8  1.72  0.152
    ## 1                                     3 -177.559 361.8  1.77  0.148
    ## 7                                     5 -175.115 362.2  2.09  0.126
    ## 5                                     4 -176.950 363.1  3.08  0.077
    ## 14                1.889               6 -174.295 363.4  3.32  0.068
    ## 2                                     4 -177.518 364.3  4.22  0.044
    ## 6                                     5 -176.765 365.5  5.40  0.024
    ## Models ranked by AICc(x)

Selected models include just inhabitants OR inhabitants in interaction
with income OR the model without effect of any predictors. The inclusion
of the this last model probably happened because of the significant
reduction in sample size (14 watersheds were removed).
