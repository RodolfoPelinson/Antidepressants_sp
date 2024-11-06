#pred_mod_sel <- select(predictors, area, casas_sem_saneamento, renda, habitantes)
#pred_mod_sel <- select(predictors,  casas_sem_saneamento, renda, habitantes)
library(vegan)

#pred_mod_sel_low_income <- pred_mod_sel[pred_mod_sel$renda <= median(pred_mod_sel$renda),]

all <- responses$all
predictors_cor <- predictors
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "Impermeabilization", "Urban_cover", "Houses_without_sanitation", "Average_Income","Inhabitants" , "Inhabitants_without_sanitation")
colnames(predictors_cor) <- pred_names

#Padronizando
Preditoras_stand <- decostand(predictors_cor, method = "stand")
nrow(Preditoras_stand)


##################################################################################################
library(randomForest)

RFmodel <- randomForest(all ~ Average_Income + Inhabitants + Houses_without_sanitation,
                        data = Preditoras_stand, proximity = TRUE, ntree = 100000,
                        mtry = 3, importance = TRUE)
RFmodel

varImpPlot(RFmodel)

partialPlot(RFmodel, Preditoras_stand, "Houses_without_sanitation")
partialPlot(RFmodel, Preditoras_stand, "Average_Income")
partialPlot(RFmodel, Preditoras_stand, "Inhabitants")

###################################################################################################


names_v <- names(attr(Preditoras_stand, "parameters")$scale)

scale_casas <- attr(Preditoras_stand, "parameters")$scale[names_v == "casas_sem_saneamento"]
center_casas <- attr(Preditoras_stand, "parameters")$center[names_v == "casas_sem_saneamento"]

scale_renda <- attr(Preditoras_stand, "parameters")$scale[names_v == "renda"]
center_renda <- attr(Preditoras_stand, "parameters")$center[names_v == "renda"]

center_habitantes <- attr(Preditoras_stand, "parameters")$center[names_v == "habitantes"]
scale_habitantes <- attr(Preditoras_stand, "parameters")$scale[names_v == "habitantes"]



quantile20_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.2)
quantile50_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.5)
quantile80_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.8)

quantile20_renda <- quantile(Preditoras_stand$renda, probs = 0.2)
quantile50_renda <- quantile(Preditoras_stand$renda, probs = 0.5)
quantile80_renda <- quantile(Preditoras_stand$renda, probs = 0.8)

quantile20_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.2)
quantile50_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.5)
quantile80_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.8)


new_data_20_saneamento_20_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          casas_sem_saneamento = rep(quantile20_saneamento, 100 ),
                          renda = rep(quantile20_renda, 100 ))

new_data_20_saneamento_50_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile20_saneamento, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_20_saneamento_80_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile20_saneamento, 100 ),
                                              renda = rep(quantile80_renda, 100 ))



new_data_50_saneamento_20_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile50_saneamento, 100 ),
                                              renda = rep(quantile20_renda, 100 ))

new_data_50_saneamento_50_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile50_saneamento, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_50_saneamento_80_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile50_saneamento, 100 ),
                                              renda = rep(quantile80_renda, 100 ))



new_data_80_saneamento_20_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile80_saneamento, 100 ),
                                              renda = rep(quantile20_renda, 100 ))

new_data_80_saneamento_50_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile80_saneamento, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_80_saneamento_80_renda <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                                              casas_sem_saneamento = rep(quantile80_saneamento, 100 ),
                                              renda = rep(quantile80_renda, 100 ))










new_data_20_habitantes_20_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile20_habitantes, 100 ),
                                              renda = rep(quantile20_renda, 100 ))

new_data_20_habitantes_50_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile20_habitantes, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_20_habitantes_80_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile20_habitantes, 100 ),
                                              renda = rep(quantile80_renda, 100 ))



new_data_50_habitantes_20_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile50_habitantes, 100 ),
                                              renda = rep(quantile20_renda, 100 ))

new_data_50_habitantes_50_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile50_habitantes, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_50_habitantes_80_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile50_habitantes, 100 ),
                                              renda = rep(quantile80_renda, 100 ))



new_data_80_habitantes_20_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile80_habitantes, 100 ),
                                              renda = rep(quantile20_renda, 100 ))

new_data_80_habitantes_50_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile80_habitantes, 100 ),
                                              renda = rep(quantile50_renda, 100 ))

new_data_80_habitantes_80_renda <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                                              habitantes = rep(quantile80_habitantes, 100 ),
                                              renda = rep(quantile80_renda, 100 ))





predicted_saneamento_20_renda_20 <- predict(RFmodel, newdata = new_data_20_saneamento_20_renda, type = "response")
predicted_saneamento_20_renda_50 <- predict(RFmodel, newdata = new_data_20_saneamento_50_renda, type = "response")
predicted_saneamento_20_renda_80 <- predict(RFmodel, newdata = new_data_20_saneamento_80_renda, type = "response")
predicted_saneamento_50_renda_20 <- predict(RFmodel, newdata = new_data_50_saneamento_50_renda, type = "response")
predicted_saneamento_50_renda_50 <- predict(RFmodel, newdata = new_data_50_saneamento_50_renda, type = "response")
predicted_saneamento_50_renda_80 <- predict(RFmodel, newdata = new_data_50_saneamento_50_renda, type = "response")
predicted_saneamento_80_renda_20 <- predict(RFmodel, newdata = new_data_80_saneamento_20_renda, type = "response")
predicted_saneamento_80_renda_50 <- predict(RFmodel, newdata = new_data_80_saneamento_50_renda, type = "response")
predicted_saneamento_80_renda_80 <- predict(RFmodel, newdata = new_data_80_saneamento_20_renda, type = "response")



predicted_habitantes_20_renda_20 <- predict(RFmodel, newdata = new_data_20_habitantes_20_renda, type = "response")
predicted_habitantes_20_renda_50 <- predict(RFmodel, newdata = new_data_20_habitantes_50_renda, type = "response")
predicted_habitantes_20_renda_80 <- predict(RFmodel, newdata = new_data_20_habitantes_80_renda, type = "response")
predicted_habitantes_50_renda_20 <- predict(RFmodel, newdata = new_data_50_habitantes_50_renda, type = "response")
predicted_habitantes_50_renda_50 <- predict(RFmodel, newdata = new_data_50_habitantes_50_renda, type = "response")
predicted_habitantes_50_renda_80 <- predict(RFmodel, newdata = new_data_50_habitantes_50_renda, type = "response")
predicted_habitantes_80_renda_20 <- predict(RFmodel, newdata = new_data_80_habitantes_20_renda, type = "response")
predicted_habitantes_80_renda_50 <- predict(RFmodel, newdata = new_data_80_habitantes_50_renda, type = "response")
predicted_habitantes_80_renda_80 <- predict(RFmodel, newdata = new_data_80_habitantes_20_renda, type = "response")


casas_sem_saneamento_rescaled <- rescale(Preditoras_stand$casas_sem_saneamento, scale_casas,center_casas)
renda_rescaled <- rescale(Preditoras_stand$renda, scale_renda,center_renda)
habitantes_rescaled <- rescale(Preditoras_stand$habitantes, scale_habitantes,center_habitantes)

new_data_habitantes <- rescale(new_data_20_saneamento_20_renda$habitantes, scale_habitantes,center_habitantes)
new_data_casas_sem_saneamento <- rescale(new_data_20_habitantes_20_renda$casas_sem_saneamento, scale_casas,center_casas)




quantile20_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.2)
quantile50_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.5)
quantile80_saneamento <- quantile(Preditoras_stand$casas_sem_saneamento, probs = 0.8)

quantile20_renda <- quantile(Preditoras_stand$renda, probs = 0.2)
quantile50_renda <- quantile(Preditoras_stand$renda, probs = 0.5)
quantile80_renda <- quantile(Preditoras_stand$renda, probs = 0.8)

quantile20_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.2)
quantile50_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.5)
quantile80_habitantes <- quantile(Preditoras_stand$habitantes, probs = 0.8)

quantile20_saneamento <- rescale(quantile20_saneamento, scale_casas, center_casas)
quantile50_saneamento <- rescale(quantile50_saneamento, scale_casas, center_casas)
quantile80_saneamento <- rescale(quantile80_saneamento, scale_casas, center_casas)

quantile20_renda <- rescale(quantile20_renda, scale_renda, center_renda)
quantile50_renda <- rescale(quantile50_renda, scale_renda, center_renda)
quantile80_renda <- rescale(quantile80_renda, scale_renda, center_renda)

quantile20_habitantes <- rescale(quantile20_habitantes, scale_habitantes, center_habitantes)
quantile50_habitantes <- rescale(quantile50_habitantes, scale_habitantes, center_habitantes)
quantile80_habitantes <- rescale(quantile80_habitantes, scale_habitantes, center_habitantes)




#png(paste(caminho, "Antidepressivos_vs_casas_sem_saneamento_by_renda_all.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)


#plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  par(mfrow = c(3,1))

  plot(responses$all ~ habitantes_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Habitantes", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_habitantes, predicted_saneamento_20_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_20_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_20_renda_80, lty = 2, col ="blue", lwd = lwd)

  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ habitantes_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
#}
  
  
  
  
  #plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ habitantes_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Habitantes", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_habitantes, predicted_saneamento_50_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_50_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_50_renda_80, lty = 2, col ="blue", lwd = lwd)
  
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ habitantes_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
  #}

  
  
  
  
  
  #plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ habitantes_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Habitantes", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_habitantes, predicted_saneamento_80_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_80_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_habitantes, predicted_saneamento_80_renda_80, lty = 2, col ="blue", lwd = lwd)
  
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ habitantes_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
  #}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #png(paste(caminho, "Antidepressivos_vs_casas_sem_habitantes_by_renda_all.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
  
  
  #plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  par(mfrow = c(3,1))
  
  plot(responses$all ~ casas_sem_saneamento_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "casas_sem_saneamento", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_20_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_20_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_20_renda_80, lty = 2, col ="blue", lwd = lwd)
  
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ casas_sem_saneamento_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
  #}
  
  
  
  
  #plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ casas_sem_saneamento_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "casas_sem_saneamento", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_50_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_50_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_50_renda_80, lty = 2, col ="blue", lwd = lwd)
  
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ casas_sem_saneamento_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
  #}
  
  
  
  
  
  
  #plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ casas_sem_saneamento_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "casas_sem_saneamento", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_80_renda_20, lty = 1, col = "red", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_80_renda_50, lty = 1, col ="black", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_habitantes_80_renda_80, lty = 2, col ="blue", lwd = lwd)
  
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50_renda){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ casas_sem_saneamento_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile20_renda)),
                    paste("R$", sprintf("%.2f",quantile50_renda)),
                    paste("R$", sprintf("%.2f",quantile80_renda))), lty = c(1,1,1), lwd = lwd, col = c("red",  "black", "blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
  #}
  

