loess_func <- function(y, x, span = 1, ...){
  
  loess <- loess(y ~ x, span=span, ...)
  newdata <- data.frame(x = seq(from = 0, to = 100, by = 1))
  loess_pred <- predict(loess, newdata = newdata)
  return(loess_pred)
}
