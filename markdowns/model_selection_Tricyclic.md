Predictors of total concentration of Tricyclic antidepressants in
streamwater
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
Preditoras_stand$Tricyclic <- responses$Tricyclic

mod_global_Tricyclic <- glmmTMB(Tricyclic ~ 
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


dredge_modelos_Tricyclic <- MuMIn::dredge(mod_global_Tricyclic, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_Tricyclic, subset = delta < 2)

aictab(models)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##    K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## 5  4 209.31       0.00   0.38   0.38 -100.22
    ## 14 6 210.19       0.88   0.25   0.63  -98.14
    ## 6  5 210.50       1.19   0.21   0.84  -99.58
    ## 1  3 211.09       1.78   0.16   1.00 -102.29

``` r
dredge_modelos_Tricyclic
```

    ## Global model call: glmmTMB(formula = Tricyclic ~ casas_sem_saneamento + habitantes + 
    ##     renda + casas_sem_saneamento:renda + habitantes:renda, data = Preditoras_stand, 
    ##     family = tweedie(link = "log"), na.action = "na.fail", ziformula = ~0, 
    ##     dispformula = ~1)
    ## ---
    ## Model selection table 
    ##    cnd((Int)) dsp((Int)) cnd(css_sem_snm) cnd(hbt) cnd(rnd)
    ## 5       2.291          +                            -0.6559
    ## 14      2.844          +           0.9537            0.2149
    ## 6       2.253          +          -0.3057           -0.7671
    ## 1       2.447          +                                   
    ## 7       2.280          +                    0.1155  -0.6429
    ## 3       2.429          +                    0.1809         
    ## 2       2.442          +          -0.1076                  
    ## 23      2.332          +                    0.2151  -0.5382
    ##    cnd(css_sem_snm:rnd) cnd(hbt:rnd) df   logLik  AICc delta weight
    ## 5                                     4 -100.219 209.3  0.00  0.298
    ## 14                2.274               6  -98.140 210.2  0.88  0.192
    ## 6                                     5  -99.582 210.5  1.19  0.165
    ## 1                                     3 -102.291 211.1  1.78  0.122
    ## 7                                     5 -100.090 211.5  2.21  0.099
    ## 3                                     4 -101.958 212.8  3.48  0.052
    ## 2                                     4 -102.203 213.3  3.97  0.041
    ## 23                            0.3219  6  -99.972 213.9  4.54  0.031
    ## Models ranked by AICc(x)

Selected models include many combinations of predictors, excluding
inhuabitants. It also includes the model with no predictors
