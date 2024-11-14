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

source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/plot_models_one_var.R")
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/plot_models_two_var.R")
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/confidence_interval.R")

pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

predictors_cor <- predictors[,-ncol(predictors)]
colnames(predictors_cor) <- pred_names
cors <- cor(predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(predictors_cor, use = "pairwise.complete.obs")

plot(predictors$habitantes ~ predictors$habitantes_sem_saneamento)
plot(predictors$habitantes ~ predictors$casas_sem_saneamento)

#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/repos/Antidepressants_sp/Resultados/Resultados Antidepressivos/Figuras/PDF/"

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


#fazendo todas as combinações possíveis

dredge_modelos_all <- MuMIn::dredge(mod_global_all, subset = smat)

models <- get.models(dredge_modelos_all, subset = delta < 2)

aictab(models)

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all.csv", sep = ",")



#################################################################################################
best_model <- models[[1]]
second_model <- models[[2]]

pdf(paste(caminho, "Antidepressivos_vs_saneamento.pdf", sep = ""), width = 4, height = 3.5, pointsize = 7)
plot_models_two_var(model = best_model, var1 = pred_mod_sel$casas_sem_saneamento, var2 = pred_mod_sel$renda, trans_renda = 1.7485)
dev.off()

pdf(paste(caminho, "Antidepressivos_vs_habitantes.pdf", sep = ""), width = 4, height = 3.5, pointsize = 7)
plot_models_one_var(second_model, var1 = pred_mod_sel$habitantes)
dev.off()

pdf(paste(caminho, "Paper/Antidepressivos_vs_habitantes_saneamento.pdf", sep = ""), width = 8, height = 3.5, pointsize = 8)
par(mfrow = c(1,2), mar = c(3.5,5,1,0.1), bty = "l")
plot_models_two_var(model = best_model, var1 = pred_mod_sel$casas_sem_saneamento, var2 = pred_mod_sel$renda, trans_renda = 1.7485)
par(new = TRUE, mar = c(0,0,0,0))
plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
text(x = 0, y = 100, labels = "A)", font = 2, cex = 2, adj = 0)

par(mar = c(3.5,5,1,0.1))
plot_models_one_var(second_model, var1 = pred_mod_sel$habitantes)
par(new = TRUE, mar = c(0,0,0,0))
plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
text(x = 0, y = 100, labels = "B)", font = 2, cex = 2, adj = 0)
dev.off()

