#Abordagem fazendo Seleção de modelos

#Selecionando as melhores variaveis locais

library(DHARMa)
library(AICcmodavg)
library(corrplot)
library(MuMIn)
library(vegan)
library(car)
library(lme4)
library(glmmTMB)
library(dplyr)

caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"


#pred_mod_sel <- select(predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(predictors, casas_sem_saneamento, renda, habitantes)

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

Preditoras_stand$SNRI <- responses$SNRI



#Modelol global
mod_global_SNRI <- glmmTMB(SNRI ~ 
                            # area +
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

dredge_modelos_SNRI <- MuMIn::dredge(mod_global_SNRI, subset = smat)

models <- get.models(dredge_modelos_SNRI, subset = delta < 2)

aictab(models)


write.table(data.frame(as.matrix(dredge_modelos_SNRI)), "Resultados/model_selection_SNRI.csv", sep = ",")







##################################################################################################
best_model_habitantes <- models$`23`

rescale <- function(x, scale, center){
  (x * scale) + center
}

names_v <- names(attr(Preditoras_stand, "parameters")$scale)

scale_habitantes <- attr(Preditoras_stand, "parameters")$scale[names_v == "habitantes"]
center_habitantes <- attr(Preditoras_stand, "parameters")$center[names_v == "habitantes"]

scale_area <- attr(Preditoras_stand, "parameters")$scale[names_v == "area"]
center_area <- attr(Preditoras_stand, "parameters")$center[names_v == "area"]

scale_renda <- attr(Preditoras_stand, "parameters")$scale[names_v == "renda"]
center_renda <- attr(Preditoras_stand, "parameters")$center[names_v == "renda"]


##################################################################################################




##################################################################################################

##################################################################################################

quantile10 <- quantile(Preditoras_stand$renda, probs = 0.1)
quantile25 <- quantile(Preditoras_stand$renda, probs = 0.25)
quantile50 <- quantile(Preditoras_stand$renda, probs = 0.5)
quantile75 <- quantile(Preditoras_stand$renda, probs = 0.75)
quantile90 <- quantile(Preditoras_stand$renda, probs = 0.9)



new_data_10 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          renda = rep(quantile10, 100 ))
new_data_25 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          renda = rep(quantile25, 100 ))
new_data_50 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          renda = rep(quantile50, 100 ))
new_data_75 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          renda = rep(quantile75, 100 ))
new_data_90 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          renda = rep(quantile90, 100 ))

predicted_renda_habitantes_10 <- predict(best_model_habitantes, newdata = new_data_10, type = "response")
predicted_renda_habitantes_25 <- predict(best_model_habitantes, newdata = new_data_25, type = "response")
predicted_renda_habitantes_50 <- predict(best_model_habitantes, newdata = new_data_50, type = "response")
predicted_renda_habitantes_75 <- predict(best_model_habitantes, newdata = new_data_75, type = "response")
predicted_renda_habitantes_90 <- predict(best_model_habitantes, newdata = new_data_90, type = "response")


habitantes_rescaled <- rescale(Preditoras_stand$habitantes, scale_habitantes,center_habitantes)
renda_rescaled <- rescale(Preditoras_stand$renda, scale_renda,center_renda)
new_data_habitantes <- rescale(new_data_10$habitantes, scale_habitantes,center_habitantes)

quantile10 <- rescale(quantile10, scale_renda, center_renda)
quantile25 <- rescale(quantile25, scale_renda, center_renda)
quantile50 <- rescale(quantile50, scale_renda, center_renda)
quantile75 <- rescale(quantile75, scale_renda, center_renda)
quantile90 <- rescale(quantile90, scale_renda, center_renda)



png(paste(caminho, "Antidepressivos_vs_habitantes_by_renda_SNRI.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mar = c(3.5,4,1,0.1), bty = "l")

plot(responses$SNRI ~ habitantes_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(xlab = "Number of inhabitants", line = 2, cex.lab = 1.25)
title(ylab = "SNRI Antidepressant", line = 2.5, cex.lab = 1.25)
title(ylab = "concentration (ng/L)", line = 1.5, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data_habitantes, predicted_renda_habitantes_50, lty = 1, lwd = 2)
lines(new_data_habitantes, predicted_renda_habitantes_25, lty = 1, col ="red", lwd = 2)
lines(new_data_habitantes, predicted_renda_habitantes_10, lty = 2, col ="red", lwd = 2)
lines(new_data_habitantes, predicted_renda_habitantes_75, lty = 1, col ="blue", lwd = 2)
lines(new_data_habitantes, predicted_renda_habitantes_90, lty = 2, col ="blue", lwd = 2)

cols <- rep(NA, length(renda_rescaled))
for(i in 1:length(renda_rescaled)){
  if(renda_rescaled[i] < quantile50){
    cols[i] <- "red"
  }else{cols[i] <- "blue"}
}

renda_size <- renda_rescaled/100

points(responses$SNRI ~ habitantes_rescaled, pch = 1, cex = log(renda_size), col = cols)

legend(legend = c(paste("R$", sprintf("%.2f",quantile10)),
                  paste("R$", sprintf("%.2f",quantile25)),
                  paste("R$", sprintf("%.2f",quantile50)),
                  paste("R$", sprintf("%.2f",quantile75)),
                  paste("R$", sprintf("%.2f",quantile90))), lty = c(2,1,1,1,2), lwd = 2, col = c("red", "red", "black", "blue","blue"),
       x = 40000, y = 300, bty = "n")

dev.off()







