
plot_models_two_var <- function(model, var1, var2, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),trans_renda = 1,
                        xlab = "Number of houses with no sanitation",
                        ylab1 = "Antidepressant",
                        ylab2 = "concentration (ng/L)",
                        lwd = 2,
                        confint = FALSE,
                        legend_x = 2000,
                        legend_y = 700){

  library(scales)
  
  response <- model$frame[,1]
  
  Preditoras_stand <- data.frame(var1,var2)
  library(vegan)
  Preditoras_stand <- decostand(Preditoras_stand, method = "stand")
  
  rescale <- function(x, scale, center){
    (x * scale) + center
  }
  
  var1<-Preditoras_stand$var1
  var2<-Preditoras_stand$var2

  
  names_v <- names(attr(Preditoras_stand, "parameters")$scale)
  
  scale_var1 <- attr(Preditoras_stand, "parameters")$scale[names_v == "var1"]
  center_var1 <- attr(Preditoras_stand, "parameters")$center[names_v == "var1"]
  
  scale_var2 <- attr(Preditoras_stand, "parameters")$scale[names_v == "var2"]
  center_var2 <- attr(Preditoras_stand, "parameters")$center[names_v == "var2"]
  
  
  quantile10 <- quantile(var2, probs = probs[1])
  quantile25 <- quantile(var2, probs = probs[2])
  quantile50 <- quantile(var2, probs = probs[3])
  quantile75 <- quantile(var2, probs = probs[4])
  quantile90 <- quantile(var2, probs = probs[5])
  
  
  new_data_10 <- data.frame(var1 = seq(min(var1), max(var1), length.out = 100),
                            var2 = rep(quantile10, 100 ))
  new_data_25 <- data.frame(var1 = seq(min(var1), max(var1), length.out = 100),
                            var2 = rep(quantile25, 100 ))
  new_data_50 <- data.frame(var1 = seq(min(var1), max(var1), length.out = 100),
                            var2 = rep(quantile50, 100 ))
  new_data_75 <- data.frame(var1 = seq(min(var1), max(var1), length.out = 100),
                            var2 = rep(quantile75, 100 ))
  new_data_90 <- data.frame(var1 = seq(min(var1), max(var1), length.out = 100),
                            var2 = rep(quantile90, 100 ))
  
  
  var_names <- colnames(model$frame[,-1])
   
  new_data_10_2 <- new_data_10
  new_data_25_2 <- new_data_25
  new_data_50_2 <- new_data_50
  new_data_75_2 <- new_data_75
  new_data_90_2 <- new_data_90 
  
  colnames(new_data_10_2) <- var_names
  colnames(new_data_25_2) <- var_names
  colnames(new_data_50_2) <- var_names
  colnames(new_data_75_2) <- var_names
  colnames(new_data_90_2) <- var_names
  
  predicted_var2_10 <- confidence_interval(model, newdata = new_data_10_2)
  predicted_var2_25 <- confidence_interval(model, newdata = new_data_25_2)
  predicted_var2_50 <- confidence_interval(model, newdata = new_data_50_2)
  predicted_var2_75 <- confidence_interval(model, newdata = new_data_75_2)
  predicted_var2_90 <- confidence_interval(model, newdata = new_data_90_2)
  
  
  var1_rescaled <- rescale(var1, scale_var1,center_var1)
  var2_rescaled <- rescale(var2, scale_var2,center_var2)
  new_data_var1 <- rescale(new_data_10$var1, scale_var1,center_var1)
  
  quantile10 <- rescale(quantile10, scale_var2, center_var2) / trans_renda
  quantile25 <- rescale(quantile25, scale_var2, center_var2) / trans_renda
  quantile50 <- rescale(quantile50, scale_var2, center_var2) / trans_renda
  quantile75 <- rescale(quantile75, scale_var2, center_var2) / trans_renda
  quantile90 <- rescale(quantile90, scale_var2, center_var2) / trans_renda
  
  
  
  #png(paste(caminho, "Antidepressivos_vs_var1_by_var2_all.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
  
  
    par(mar = c(3.5,5,1,0.1), bty = "l")
    
    plot(response ~ var1_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    
    title(xlab = xlab, line = 2, cex.lab = 1.25)
    title(ylab = ylab1, line = 3, cex.lab = 1.25)
    title(ylab = ylab2, line = 2, cex.lab = 1.25)
    
    axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
    axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
    
    if(isTRUE(confint)){
      polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var2_50$lower,predicted_var2_50$upper[100:1]), col = alpha("black", 0.2), border = FALSE)
      polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var2_25$lower,predicted_var2_25$upper[100:1]), col = alpha("red", 0.2), border = FALSE)
      polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var2_10$lower,predicted_var2_10$upper[100:1]), col = alpha("red", 0.2), border = FALSE)
      polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var2_75$lower,predicted_var2_75$upper[100:1]), col = alpha("blue", 0.2), border = FALSE)
      polygon(x = c(new_data_var1,new_data_var1[100:1]), y = c(predicted_var2_90$lower,predicted_var2_90$upper[100:1]), col = alpha("blue", 0.2), border = FALSE)
    }
    

    
    lines(new_data_var1, predicted_var2_50$fit, lty = 1, lwd = lwd)
    lines(new_data_var1, predicted_var2_25$fit, lty = 1, col ="firebrick1", lwd = lwd)
    lines(new_data_var1, predicted_var2_10$fit, lty = 2, col ="firebrick1", lwd = lwd)
    lines(new_data_var1, predicted_var2_75$fit, lty = 1, col ="royalblue", lwd = lwd)
    lines(new_data_var1, predicted_var2_90$fit, lty = 2, col ="royalblue", lwd = lwd)
    
    cols <- rep(NA, length(var2_rescaled))
    bg <- rep(NA, length(var2_rescaled))
    
    for(i in 1:length(var2_rescaled)){
      if(var2_rescaled[i] < quantile50*trans_renda){
        cols[i] <- "firebrick1"
        bg[i] <- alpha("firebrick1", 0.4)
      }else{
        cols[i] <- "royalblue"
        bg[i] <- alpha("royalblue", 0.4)
      }
    }
    
    var2_size <- var2 + (min(var2)*-1)
    
    points(response ~ var1_rescaled, pch = 16, cex = var2_size+1, col = bg, lwd = lwd *0.7)
    
    legend(legend = c(paste("US$", sprintf("%.2f",quantile10)),
                      paste("US$", sprintf("%.2f",quantile25)),
                      paste("US$", sprintf("%.2f",quantile50)),
                      paste("US$", sprintf("%.2f",quantile75)),
                      paste("US$", sprintf("%.2f",quantile90))), lty = c(2,
                                                                         1,
                                                                         1,
                                                                         1,
                                                                         2), lwd = lwd, col = c("firebrick1",
                                                                                                        "firebrick1",
                                                                                                        "black",
                                                                                                        "royalblue",
                                                                                                        "royalblue"),
           x = legend_x, y = legend_y, bty = "n")
    
    #par(new = TRUE, mar = c(0,0,0,0))
    #plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
    #text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
    #par(mar = margins)
    
}


