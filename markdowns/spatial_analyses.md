Spatial Analyses
================
Rodolfo Pelinson
2024-11-06

## Plots of the spatial structure of predictors and response variables

Functions necessary for the following plots:

``` r
plot_spatial_pattern <- function(y, main = ""){
  matrix <- matrix(data = c(0,.8,0,1,
                          .8,1,0,1), nrow = 2, ncol = 4, byrow = TRUE)

screens <- split.screen(matrix)

screen(screens[1])
par(mar = c(4,4,3,1))

pal <- col_numeric(
        palette = c(alpha("blue", 0.6), alpha("red", 0.6)),
        domain = y,
        na.color = "grey50",
        alpha = TRUE,
        reverse = FALSE)

y_col <- pal(y)

plot(coord$longitude, coord$latitude, xlab = "Longitud", ylab = "Latitud", pch = 21, bg = y_col, cex = 3, col = FALSE, main = main)

legend<-seq(min(y), max(y), length.out = 100)

screen(screens[2])

par(mar = c(4,2,3,4))
image(y = 1:length(legend),x = 1, z = t(as.matrix(1:length(legend))),
          col = pal(legend),
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
axis(4, at = c(1, length(legend)/2,length(legend)), labels = round(legend[c(1,length(legend)/2,length(legend))],1), las = 2,
         cex.axis = 1, tick = TRUE, line = 0)
close.screen(all.screens = TRUE)
}
```

``` r
colapsed_variogram <- function(variogram, dist_classes = 10){
  variogram <- variogram[order(variogram$dist),]
  new_seq <- seq(min(variogram$dist), max(variogram$dist), length.out = dist_classes)
  new_variogram <- rep(NA, length(new_seq)-1)
  for(i in 1:(length(new_seq) - 1)){
    new_variogram[i] <- mean(variogram$variog[variogram$dist >= new_seq[i] & variogram$dist < new_seq[i+1]])
  }
  new_variogram <- data.frame(variog = new_variogram, dist = new_seq[-length(new_seq)])
  class(new_variogram) <- class(variogram)
  return(new_variogram)
}
```

This is just the spatial pattern of the most important predictor
variables that we have.

``` r
plot_spatial_pattern(predictors$casas_sem_saneamento, main = "Houses without sanitation")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot_spatial_pattern(predictors$habitantes, main = "Inhanitants")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
plot_spatial_pattern(predictors$renda, main = "Per capta average income")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

As we can see, these variables do have some opposite spatial patterns.
Specifically, the total number of houses without sanitation and total
number of inhabitants exhibit opposite patterns to the per capta average
income.

``` r
plot_spatial_pattern(responses$all, main = "Antidepressant concentration")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plot_spatial_pattern(responses$SNRI, main = "SNRI Antidepressant concentration")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
plot_spatial_pattern(responses$Aminoketone, main = "Aminoketone Antidepressant concentration")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
plot_spatial_pattern(responses$SSRI, main = "SSRI Antidepressant concentration")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
plot_spatial_pattern(responses$Tricyclic, main = "Tricyclic Antidepressant concentration")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

These are the spatial patterns found for antidepressant concentration.
It seems that the ones with high concentration are concentrated in the
center of the distribution.

## Analysis for the total concentration of antidepressants

First lets find the best model to predict the total concentration of
antidepressants

``` r
#pred_mod_sel <- select(predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(predictors,  casas_sem_saneamento, renda, habitantes)

#Padronizando
Preditoras_stand <- decostand(pred_mod_sel, method = "stand")

#construindo uma matriz que indica quais combinações de variaveis tem correlação acima de 0.5 ou abaixo de -0.5
smat <- abs(cor(Preditoras_stand)) <= .5
smat[!lower.tri(smat)] <- NA

colnames(smat) <- c(paste("cond(",colnames(smat)[1], ")", sep = ""),
                    paste("cond(",colnames(smat)[2], ")", sep = ""),
                    paste("cond(",colnames(smat)[3], ")", sep = ""))
                    #paste("cond(",colnames(smat)[4], ")", sep = ""))

rownames(smat) <- colnames(smat)



###############################################
###############################################
####### TODOS OS ANTIDEPRESSIVOS ##############

Preditoras_stand$all <- responses$all

#Modelol global
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


#fazendo todas as combinações possíveis

dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)
```

    ## Fixed terms are "cond((Int))" and "disp((Int))"

``` r
models <- get.models(dredge_modelos_all, subset = delta < 2)

best_model <- models[[1]]

best_model
```

    ## Formula:          
    ## all ~ casas_sem_saneamento + renda + casas_sem_saneamento:renda +      1
    ## Data: Preditoras_stand
    ##       AIC       BIC    logLik  df.resid 
    ##  547.2504  558.8413 -267.6252        45 
    ## 
    ## Number of obs: 51
    ## 
    ## Dispersion parameter for tweedie family (): 27.4 
    ## 
    ## Tweedie power estimate:  1.4 
    ## 
    ## Fixed Effects:
    ## 
    ## Conditional model:
    ##                (Intercept)        casas_sem_saneamento  
    ##                     5.8056                      0.9669  
    ##                      renda  casas_sem_saneamento:renda  
    ##                     0.7120                      2.0802

With the best model found, we will check for spatial patterns in a model
without predictors, in our best model, and in our best model with the
addition of spatial autocorrelation.

``` r
set.seed(45)
library(DHARMa)
library(nlme)

Preditoras_stand$pos <- numFactor(coord$longitude, coord$latitude)
Preditoras_stand$group <- rep("A", nrow(Preditoras_stand))

mod_no_effect <- glmmTMB(all ~ 1,data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")
mod_auto_correlation <- glmmTMB(all ~ casas_sem_saneamento * renda + exp(pos + 0|group),data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")

resid_no_effect <- simulateResiduals(mod_no_effect, n = 10000, rotation = "estimated")
resid_all <- simulateResiduals(best_model, n = 10000, rotation = "estimated")
resid_auto_correlation <- simulateResiduals(mod_auto_correlation, n = 10000, rotation = "estimated")

plot(resid_no_effect)
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
plot(resid_all)
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
plot(resid_auto_correlation)
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
plot_spatial_pattern(resid_no_effect$scaledResiduals, main = "DHARMa residuals with no prectors")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
plot_spatial_pattern(resid_all$scaledResiduals, main = "DHARMa residuals for the best model")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->

``` r
plot_spatial_pattern(resid_auto_correlation$scaledResiduals, main = "DHARMa residuals modeling spatial autocorrelation")
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->

``` r
var_no_effect <- Variogram(resid_no_effect$scaledResiduals, distance = dist(coord[,2:3]))
var_no_effect_colapsed <- colapsed_variogram(var_no_effect)

var <- Variogram(resid_all$scaledResiduals, distance = dist(coord[,2:3]))
var_colapsed <- colapsed_variogram(var)

var_auto_cor <- Variogram(resid_auto_correlation$scaledResiduals, distance = dist(coord[,2:3]))
var_auto_cor_colapsed <- colapsed_variogram(var_auto_cor)

plot(var_no_effect_colapsed)
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->

``` r
plot(var_colapsed)
```

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->

``` r
plot(var_auto_cor_colapsed)
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = FALSE, :
    ## pseudoinverse used at 0.75633

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = FALSE, :
    ## neighborhood radius 0.25139

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = FALSE, :
    ## reciprocal condition number 0

![](spatial_analyses_files/figure-gfm/unnamed-chunk-7-9.png)<!-- -->

## Conclusion

As we can see, our model do eliminate most of the spatial patterns
observed in our response variable. When we look at the Semivariograms we
can see that our best model is sufficient to elimate a pattern of
increase in semivariance with distance. In fact, we are left with a
pattern of decrease in semivariance for samples that are very far apart.
Meaning that samples far apart are more similar that samples close to
each other. The probable explanation for that is that forested
watersheds mostly occur at the margins of our sample design.
Unfortunately this is could be accounted for the available spatial
aucorrelation structures.
