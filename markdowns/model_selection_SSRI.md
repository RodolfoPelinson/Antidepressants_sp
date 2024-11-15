Predictors of total concentration of SSRI antidepressants in streamwater
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
```

Performing the model selection procedure. Note that the name of the
variables are in portuguese.

Habitantes mean inhabitants casas_sem_saneamento means houses without
sanitation renda means income.

First we make a matrix to check which pairs of variables are correlated
(more than 50%) and therefore should not be together in the same model.

``` r
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/ajeitando_planilhas.R")

pred_mod_sel <- select(predictors,  casas_sem_saneamento, renda, habitantes)

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
    ## casas_sem_saneamento            1.0000000 -0.3425484  0.5866546
    ## renda                          -0.3425484  1.0000000 -0.2328670
    ## habitantes                      0.5866546 -0.2328670  1.0000000

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
Preditoras_stand$SSRI <- responses$SSRI

mod_global_SSRI <- glmmTMB(SSRI ~ 
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


dredge_modelos_SSRI <- MuMIn::dredge(mod_global_SSRI, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_SSRI, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##   K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 3 4 375.15          0      1      1 -183.14

``` r
dredge_modelos_SSRI
```

    ## Global model call: glmmTMB(formula = SSRI ~ casas_sem_saneamento + habitantes + 
    ##     renda + casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt) cnd(rnd)
    ## 3       4.194          +                    0.4503         
    ## 7       4.181          +                    0.4181  -0.1858
    ## 1       4.312          +                                   
    ## 14      4.999          +          1.30900            0.9555
    ## 5       4.270          +                            -0.3339
    ## 23      4.190          +                    0.4293  -0.1705
    ## 2       4.311          +         -0.04909                  
    ## 6       4.260          +         -0.14050           -0.3886
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 3                                     4 -183.140 375.1  0.00  0.429
    ## 7                                     5 -182.949 377.2  2.08  0.152
    ## 1                                     3 -185.492 377.5  2.35  0.133
    ## 14                2.796               6 -182.078 378.1  2.92  0.100
    ## 5                                     4 -184.867 378.6  3.45  0.076
    ## 23                           0.04223  6 -182.946 379.8  4.65  0.042
    ## 2                                     4 -185.470 379.8  4.66  0.042
    ## 6                                     5 -184.699 380.7  5.58  0.026
    ## Models ranked by AICc(x)

Selected models include only inhabitants as predictor.

Plotting the best models:

``` r
best_model <- models[[1]]

par(mar = c(3.5,5,1,0.1))
plot_models_one_var(best_model, var1 = pred_mod_sel$habitantes, xlab = "Number of inhabitants")
```

<img src="model_selection_SSRI_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />
