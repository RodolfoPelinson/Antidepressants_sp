
plot_models_one_var <- function(model, var1,
                                xlab = "Number of inhabitants",
                                ylab1 = "Antidepressant",
                                ylab2 = "concentration (ng/L)",
                                lwd = 2,
                                confint = FALSE){
library(scales)
response <- model$frame[,1]

predictors <- data.frame(var1 = var1)

library(vegan)
Preditoras_stand <- decostand(predictors, method = "stand")

rescale <- function(x, scale, center){
  (x * scale) + center
}


names_v <-colnames(model$frame)[-1]

scale_var1 <- attr(Preditoras_stand, "parameters")$scale
center_var1 <- attr(Preditoras_stand, "parameters")$center

var1<-Preditoras_stand$var1


new_data <- data.frame(var1 = seq(min(Preditoras_stand$var1), max(Preditoras_stand$var1), length.out = 100))

new_data_2 <- new_data
colnames(new_data_2) <- names_v

predicted_var1 <- confidence_interval(model = model, newdata = new_data_2)

new_data_var1 <- rescale(new_data$var1, scale_var1,center_var1)


var1_rescaled <- rescale(var1, scale_var1,center_var1)


  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(response ~ var1_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Number of inhabitants", line = 2., cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  if(isTRUE(confint)){
    polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var1$lower,predicted_var1$upper[100:1]), col = alpha("black", 0.2), border = FALSE)
    
  }
  
  lines(new_data_var1, predicted_var1$fit, lty = 1, lwd = lwd)
  
  #lines(new_data_var1, predicted_var1_50, lty = 1, lwd = 2)
  #lines(new_data_var1, predicted_var1_25, lty = 1, col ="red", lwd = 2)
  #lines(new_data_var1, predicted_var1_10, lty = 2, col ="red", lwd = 2)
  #lines(new_data_var1, predicted_var1_75, lty = 1, col ="blue", lwd = 2)
  #lines(new_data_var1, predicted_var1_90, lty = 2, col ="blue", lwd = 2)
  
  #cols <- rep(NA, length(area_rescaled))
  #for(i in 1:length(area_rescaled)){
  #  if(area_rescaled[i] < quantile50){
  #    cols[i] <- "red"
  #  }else{cols[i] <- "blue"}
  #}
  
  #area_size <- area_rescaled/100
  
  #points(responses$all ~ var1_rescaled, pch = 1, cex = log(area_size)+1, col = cols)
  
  points(response ~ var1_rescaled, pch = 16, cex = 2, col = alpha("black", 0.4), lwd = lwd *0.7)
  
  
  #legend(legend = c(paste(quantile10, "ha"),
  #                  paste(quantile25, "ha"),
  #                  paste(quantile50, "ha"),
  #                  paste(quantile75, "ha"),
  #                  paste(quantile90, "ha")), lty = c(2,1,1,1,2), lwd = 2, col = c("red", "red", "black", "blue","blue"),
  #       x = 40000 , y = 200, bty = "n")
  
  #dev.off()
  #par(new = TRUE, mar = c(0,0,0,0))
  #plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  #text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  
  
  #par(mar = margins)
  
}
