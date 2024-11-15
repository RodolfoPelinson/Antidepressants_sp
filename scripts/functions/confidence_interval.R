confidence_interval <- function(model, newdata,...){
  predicted_var <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  predicted_var_se <- predicted_var$se.fit
  
  z <- qnorm(.975, mean = 0, sd = 1)
  
  predicted_var_se_upper <- predicted_var$fit + predicted_var_se*z
  predicted_var_se_lower <- predicted_var$fit - predicted_var_se*z
  
  predicted_var <- exp(predicted_var$fit)
  
  predicted_var_upper <- exp(predicted_var_se_upper)
  predicted_var_lower <- exp(predicted_var_se_lower)
  
  confs <- list(fit = predicted_var, upper = predicted_var_upper, lower = predicted_var_lower)
  return(confs)
}
