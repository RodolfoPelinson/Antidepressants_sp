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

Preditoras_stand$Aminoketone <- responses$Aminoketone




#Modelol global
mod_global_Aminoketone <- glmmTMB(Aminoketone ~ 
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

dredge_modelos_Aminoketone <- MuMIn::dredge(mod_global_Aminoketone, subset = smat)

models <- get.models(dredge_modelos_Aminoketone, subset = delta < 2)

aictab(models)

write.table(data.frame(as.matrix(dredge_modelos_Aminoketone)), "Resultados/model_selection_Aminoketone.csv", sep = ",")



##################################################################################################
best_model_casas <- models$`14`

mod_best <- glmmTMB(Aminoketone ~   casas_sem_saneamento +
                      renda + 
                      casas_sem_saneamento:renda,
                    data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")

mod_best2 <- glmmTMB(Aminoketone ~   casas_sem_saneamento *renda,
                     data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")

rescale <- function(x, scale, center){
  (x * scale) + center
}

names_v <- names(attr(Preditoras_stand, "parameters")$scale)

scale_casas <- attr(Preditoras_stand, "parameters")$scale[names_v == "casas_sem_saneamento"]
center_casas <- attr(Preditoras_stand, "parameters")$center[names_v == "casas_sem_saneamento"]

scale_renda <- attr(Preditoras_stand, "parameters")$scale[names_v == "renda"]
center_renda <- attr(Preditoras_stand, "parameters")$center[names_v == "renda"]


##################################################################################################



quantile10 <- quantile(Preditoras_stand$renda, probs = 0.1)
quantile25 <- quantile(Preditoras_stand$renda, probs = 0.25)
quantile50 <- quantile(Preditoras_stand$renda, probs = 0.5)
quantile75 <- quantile(Preditoras_stand$renda, probs = 0.75)
quantile90 <- quantile(Preditoras_stand$renda, probs = 0.9)


new_data_10 <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                          renda = rep(quantile10, 100 ))
new_data_25 <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                          renda = rep(quantile25, 100 ))
new_data_50 <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                          renda = rep(quantile50, 100 ))
new_data_75 <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                          renda = rep(quantile75, 100 ))
new_data_90 <- data.frame(casas_sem_saneamento = seq(min(Preditoras_stand$casas_sem_saneamento), max(Preditoras_stand$casas_sem_saneamento), length.out = 100),
                          renda = rep(quantile90, 100 ))


predicted_renda_casas_10 <- predict(best_model_casas, newdata = new_data_10, type = "response")
predicted_renda_casas_25 <- predict(best_model_casas, newdata = new_data_25, type = "response")
predicted_renda_casas_50 <- predict(best_model_casas, newdata = new_data_50, type = "response")
predicted_renda_casas_75 <- predict(best_model_casas, newdata = new_data_75, type = "response")
predicted_renda_casas_90 <- predict(best_model_casas, newdata = new_data_90, type = "response")


casas_sem_saneamento_rescaled <- rescale(Preditoras_stand$casas_sem_saneamento, scale_casas,center_casas)
renda_rescaled <- rescale(Preditoras_stand$renda, scale_renda,center_renda)
new_data_casas_sem_saneamento <- rescale(new_data_10$casas_sem_saneamento, scale_casas,center_casas)

quantile10 <- rescale(quantile10, scale_renda, center_renda)
quantile25 <- rescale(quantile25, scale_renda, center_renda)
quantile50 <- rescale(quantile50, scale_renda, center_renda)
quantile75 <- rescale(quantile75, scale_renda, center_renda)
quantile90 <- rescale(quantile90, scale_renda, center_renda)



png(paste(caminho, "Antidepressivos_vs_casas_sem_saneamento_by_renda_Aminoketone.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mar = c(3.5,4,1,0.1), bty = "l")

plot(responses$Aminoketone ~ casas_sem_saneamento_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(xlab = "Number of houses with no sanitation", line = 2, cex.lab = 1.25)
title(ylab = "Aminoketone Antidepressant", line = 2.5, cex.lab = 1.25)
title(ylab = "concentration (ng/L)", line = 1.5, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data_casas_sem_saneamento, predicted_renda_casas_50, lty = 1, lwd = 2)
lines(new_data_casas_sem_saneamento, predicted_renda_casas_25, lty = 1, col ="red", lwd = 2)
lines(new_data_casas_sem_saneamento, predicted_renda_casas_10, lty = 2, col ="red", lwd = 2)
lines(new_data_casas_sem_saneamento, predicted_renda_casas_75, lty = 1, col ="blue", lwd = 2)
lines(new_data_casas_sem_saneamento, predicted_renda_casas_90, lty = 2, col ="blue", lwd = 2)

cols <- rep(NA, length(renda_rescaled))
for(i in 1:length(renda_rescaled)){
  if(renda_rescaled[i] < quantile50){
    cols[i] <- "red"
  }else{cols[i] <- "blue"}
}

renda_size <- renda_rescaled/100

points(responses$Aminoketone ~ casas_sem_saneamento_rescaled, pch = 1, cex = log(renda_size), col = cols)

legend(legend = c(paste("R$", sprintf("%.2f",quantile10)),
                  paste("R$", sprintf("%.2f",quantile25)),
                  paste("R$", sprintf("%.2f",quantile50)),
                  paste("R$", sprintf("%.2f",quantile75)),
                  paste("R$", sprintf("%.2f",quantile90))), lty = c(2,1,1,1,2), lwd = 2, col = c("red", "red", "black", "blue","blue"),
       x = 2500, y = 250, bty = "n")

dev.off()












