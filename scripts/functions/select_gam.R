select_gam <- function(y,x, family, drop.intercept = FALSE, weights = NULL){
  
  #if(family == "nb"){
  #  family_glm <- "nbinom1"
  #}
  #if(family == "tw"){
  #  family_glm <- "tweedie"
  #}
  
  if(drop.intercept == TRUE){
    #glm_0 <- glmmTMB::glmmTMB(y ~ 1, REML = TRUE, family = family_glm)
    #glm_1 <- glmmTMB::glmmTMB(y ~ x, REML = TRUE, family = family_glm)
    gam_0 <- mgcv::gam(y ~ 1, family = family, method = "REML", select = TRUE, data = data_farmacos, weights = weights)
    gam_1 <- mgcv::gam(y ~ s(x, k = 3, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = TRUE, weights = weights)
    gam_2 <- mgcv::gam(y ~ s(x, k = 6, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = TRUE, weights = weights)
    
    if(length(y) > 13){
      gam_3 <- mgcv::gam(y ~ s(x, k = 12, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = TRUE, weights = weights)
      if(length(y) > 25){
        gam_4 <- mgcv::gam(y ~ s(x, k = 24, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = TRUE, weights = weights)
      }
    }
  }
  
  if(drop.intercept == FALSE){
    #glm_0 <- glmmTMB::glmmTMB(y ~ 1, REML = TRUE, family = family_glm)
    #glm_1 <- glmmTMB::glmmTMB(y ~ x, REML = TRUE, family = family_glm)
    gam_0 <- mgcv::gam(y ~ 1, family = family, method = "REML", select = TRUE, data = data_farmacos, weights = weights)
    gam_1 <- mgcv::gam(y ~ s(x, k = 3, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = FALSE, weights = weights)
    gam_2 <- mgcv::gam(y ~ s(x, k = 6, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = FALSE, weights = weights)
    
    if(length(y) > 13){
      gam_3 <- mgcv::gam(y ~ s(x, k = 12, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = FALSE, weights = weights)
      if(length(y) > 25){
        gam_4 <- mgcv::gam(y ~ s(x, k = 24, fx = FALSE, bs = "tp"), family = family, method = "REML", select = TRUE, drop.intercept = FALSE, weights = weights)
      }
    }
  }

  if(length(y) > 13){
    models <- list(gam_0=gam_0, gam_1=gam_1, gam_2=gam_2, gam_3=gam_3)
    aictab <- data.frame(AICctab( gam_0, gam_1, gam_2, gam_3, sort = TRUE))
    
    if(length(y) > 25){
      models <- list(gam_0=gam_0, gam_1=gam_1, gam_2=gam_2, gam_3=gam_3, gam_4=gam_4) 
      aictab <- data.frame(AICctab( gam_0, gam_1, gam_2, gam_3, gam_4, sort = TRUE))    }
  }else{
    models <- list(gam_0=gam_0, gam_1=gam_1, gam_2=gam_2) 
    aictab <- data.frame(AICctab( gam_0, gam_1, gam_2, sort = TRUE))
  }
  
  best<- models[[match( rownames(aictab)[1], names(models) )]]
  
  second_best <-  models[[match( rownames(aictab)[2], names(models) )]]
  
  
  if(family == "binomial" | family == "betar"){
    new_data <- data.frame(x = c(0:100))
    predicted <- predict(best, se = TRUE, type = "link", newdata = new_data)
    
    fit <- predicted$fit
    upr <- fit + (1.96*predicted$se.fit)
    lwr <- fit - (1.96*predicted$se.fit)
    
    fit <- 1/(1+exp(-fit))
    upr <- 1/(1+exp(-upr))
    lwr <- 1/(1+exp(-lwr))
    
    estimates_best <- data.frame(fit, upr, lwr)
    
    
    new_data <- data.frame(x = c(0:100))
    predicted <- predict(second_best, se = TRUE, type = "link", newdata = new_data)
    
    fit <- predicted$fit
    upr <- fit + (1.96*predicted$se.fit)
    lwr <- fit - (1.96*predicted$se.fit)
    
    fit <- 1/(1+exp(-fit))
    upr <- 1/(1+exp(-upr))
    lwr <- 1/(1+exp(-lwr))
    
    estimates_second <- data.frame(fit, upr, lwr)
    
  }else{
    new_data <- data.frame(x = c(0:100))
    predicted <- predict(best, se = TRUE, type = "link", newdata = new_data)
    
    fit <- predicted$fit
    upr <- exp(fit + (1.96*predicted$se.fit))
    lwr <- exp(fit - (1.96*predicted$se.fit))
    fit <- exp(fit) 
    estimates_best <- data.frame(fit, upr, lwr)
    
    
    new_data <- data.frame(x = c(0:100))
    predicted <- predict(second_best, se = TRUE, type = "link", newdata = new_data)
    
    fit <- predicted$fit
    upr <- exp(fit + (1.96*predicted$se.fit))
    lwr <- exp(fit - (1.96*predicted$se.fit))
    fit <- exp(fit) 
    estimates_second <- data.frame(fit, upr, lwr)
  }

  
  
  
  result <- list(best_model_name = rownames(as.data.frame(aictab))[which(aictab$dAICc == 0)], best_model = best, models = models,
                 AICTAB = aictab, estimates_best = estimates_best, estimates_second = estimates_second)
  return(result)
  
}

