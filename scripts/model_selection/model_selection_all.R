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

colnames(predictors)
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

predictors_cor <- predictors[,-ncol(predictors)]
colnames(predictors_cor) <- pred_names
cors <- cor(predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(predictors_cor, use = "pairwise.complete.obs")

plot(predictors$habitantes ~ predictors$habitantes_sem_saneamento)
#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"

png(paste(caminho, "Correlation_among_variables.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower")
dev.off()

png(paste(caminho, "Correlation_among_variables_significat ones.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower", sig.level = 0.05, p.mat = cors2$p, insig = "blank")
dev.off()

png(paste(caminho, "Correlation_among_variables_greater_than_0.5.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors0, method = "circle", type = "lower")
dev.off()





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



Anova(mod_global_all)

#fazendo todas as combinações possíveis

dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)

models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all.csv", sep = ",")




##################################################################################################
best_model_casas <- models$`14`

mod_best <- glmmTMB(all ~   casas_sem_saneamento +
                            renda + 
                            casas_sem_saneamento:renda,
                          data = Preditoras_stand, family = tweedie(link = "log"), na.action = "na.fail")

mod_best2 <- glmmTMB(all ~   casas_sem_saneamento *renda,
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



#png(paste(caminho, "Antidepressivos_vs_casas_sem_saneamento_by_renda_all.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)


plot_casas_renda_all  <- function(letter = "", lwd = 1.6){
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ casas_sem_saneamento_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Number of houses with no sanitation", line = 2, cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_casas_sem_saneamento, predicted_renda_casas_50, lty = 1, lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_renda_casas_25, lty = 1, col ="red", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_renda_casas_10, lty = 2, col ="red", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_renda_casas_75, lty = 1, col ="blue", lwd = lwd)
  lines(new_data_casas_sem_saneamento, predicted_renda_casas_90, lty = 2, col ="blue", lwd = lwd)
  
  cols <- rep(NA, length(renda_rescaled))
  for(i in 1:length(renda_rescaled)){
    if(renda_rescaled[i] < quantile50){
      cols[i] <- "red"
    }else{cols[i] <- "blue"}
  }
  
  renda_size <- renda_rescaled/100
  
  points(responses$all ~ casas_sem_saneamento_rescaled, pch = 1, cex = log(renda_size), col = cols, lwd = lwd *0.7)
  
  legend(legend = c(paste("R$", sprintf("%.2f",quantile10)),
                    paste("R$", sprintf("%.2f",quantile25)),
                    paste("R$", sprintf("%.2f",quantile50)),
                    paste("R$", sprintf("%.2f",quantile75)),
                    paste("R$", sprintf("%.2f",quantile90))), lty = c(2,1,1,1,2), lwd = lwd, col = c("red", "red", "black", "blue","blue"),
         x = 2000, y = 700, bty = "n")
  
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
  
}

#dev.off()














##################################################################################################
best_model_habitantes <- models$`3`

rescale <- function(x, scale, center){
  (x * scale) + center
}

names_v <- names(attr(Preditoras_stand, "parameters")$scale)

scale_habitantes <- attr(Preditoras_stand, "parameters")$scale[names_v == "habitantes"]
center_habitantes <- attr(Preditoras_stand, "parameters")$center[names_v == "habitantes"]

scale_area <- attr(Preditoras_stand, "parameters")$scale[names_v == "area"]
center_area <- attr(Preditoras_stand, "parameters")$center[names_v == "area"]


##################################################################################################




##################################################################################################


#quantile10 <- quantile(Preditoras_stand$area, probs = 0.1)
#quantile25 <- quantile(Preditoras_stand$area, probs = 0.25)
#quantile50 <- quantile(Preditoras_stand$area, probs = 0.5)
#quantile75 <- quantile(Preditoras_stand$area, probs = 0.75)
#quantile90 <- quantile(Preditoras_stand$area, probs = 0.9)

new_data <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100))

#new_data_10 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          #area = rep(quantile10, 100 ))
#new_data_25 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          #area = rep(quantile25, 100 ))
#new_data_50 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          #area = rep(quantile50, 100 ))
#new_data_75 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          #area = rep(quantile75, 100 ))
#new_data_90 <- data.frame(habitantes = seq(min(Preditoras_stand$habitantes), max(Preditoras_stand$habitantes), length.out = 100),
                          #area = rep(quantile90, 100 ))

predicted_area_habitantes <- predict(best_model_habitantes, newdata = new_data, type = "response")

#predicted_area_habitantes_10 <- predict(best_model_habitantes, newdata = new_data_10, type = "response")
#predicted_area_habitantes_25 <- predict(best_model_habitantes, newdata = new_data_25, type = "response")
#predicted_area_habitantes_50 <- predict(best_model_habitantes, newdata = new_data_50, type = "response")
#predicted_area_habitantes_75 <- predict(best_model_habitantes, newdata = new_data_75, type = "response")
#predicted_area_habitantes_90 <- predict(best_model_habitantes, newdata = new_data_90, type = "response")


habitantes_rescaled <- rescale(Preditoras_stand$habitantes, scale_habitantes,center_habitantes)
area_rescaled <- rescale(Preditoras_stand$area, scale_area,center_area)
new_data_habitantes <- rescale(new_data$habitantes, scale_habitantes,center_habitantes)

#quantile10 <- rescale(quantile10, scale_area, center_area)
#quantile25 <- rescale(quantile25, scale_area, center_area)
#quantile50 <- rescale(quantile50, scale_area, center_area)
#quantile75 <- rescale(quantile75, scale_area, center_area)
#quantile90 <- rescale(quantile90, scale_area, center_area)



#png(paste(caminho, "Antidepressivos_vs_habitantes.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

plot_habitantes_all  <- function(letter = "", lwd = 1.6){
  
  #par(mar = c(3.5,5,1,0.1), bty = "l")
  
  plot(responses$all ~ habitantes_rescaled , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  title(xlab = "Number of inhabitants", line = 2., cex.lab = 1.25)
  title(ylab = "Antidepressant", line = 3, cex.lab = 1.25)
  title(ylab = "concentration (ng/L)", line = 2, cex.lab = 1.25)
  
  axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
  axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)
  
  #polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
  lines(new_data_habitantes, predicted_area_habitantes, lty = 1, lwd = lwd)
  
  #lines(new_data_habitantes, predicted_area_habitantes_50, lty = 1, lwd = 2)
  #lines(new_data_habitantes, predicted_area_habitantes_25, lty = 1, col ="red", lwd = 2)
  #lines(new_data_habitantes, predicted_area_habitantes_10, lty = 2, col ="red", lwd = 2)
  #lines(new_data_habitantes, predicted_area_habitantes_75, lty = 1, col ="blue", lwd = 2)
  #lines(new_data_habitantes, predicted_area_habitantes_90, lty = 2, col ="blue", lwd = 2)
  
  #cols <- rep(NA, length(area_rescaled))
  #for(i in 1:length(area_rescaled)){
  #  if(area_rescaled[i] < quantile50){
  #    cols[i] <- "red"
  #  }else{cols[i] <- "blue"}
  #}
  
  #area_size <- area_rescaled/100
  
  #points(responses$all ~ habitantes_rescaled, pch = 1, cex = log(area_size)+1, col = cols)
  
  points(responses$all ~ habitantes_rescaled, pch = 1, cex = 1.5, col = "black", lwd = lwd *0.7)
  
  
  #legend(legend = c(paste(quantile10, "ha"),
  #                  paste(quantile25, "ha"),
  #                  paste(quantile50, "ha"),
  #                  paste(quantile75, "ha"),
  #                  paste(quantile90, "ha")), lty = c(2,1,1,1,2), lwd = 2, col = c("red", "red", "black", "blue","blue"),
  #       x = 40000 , y = 200, bty = "n")
  
  #dev.off()
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = 0, y = 100, labels = letter, font = 2, cex = 2, adj = 0)
  par(mar = margins)
  
}



