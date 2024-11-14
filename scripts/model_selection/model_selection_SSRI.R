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

caminho <- "C:/Users/rodol/OneDrive/repos/Antidepressants_sp/Resultados/Resultados Antidepressivos/Figuras/PDF/"


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

Preditoras_stand$SSRI <- responses$SSRI




#Modelol global
mod_global_SSRI <- glmmTMB(SSRI ~ 
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



Anova(mod_global_SSRI)

#fazendo todas as combinações possíveis

dredge_modelos_SSRI <- MuMIn::dredge(mod_global_SSRI, subset = smat)

models <- get.models(dredge_modelos_SSRI, subset = delta < 2)

aictab(models)

write.table(data.frame(as.matrix(dredge_modelos_SSRI)), "Resultados/model_selection_SSRI.csv", sep = ",")



#################################################################################################
best_model <- models[[1]]

pdf(paste(caminho, "SSRI_Antidepressivos_vs_habitantes.pdf", sep = ""), width = 4, height = 3.5, pointsize = 7)
plot_models_one_var(best_model, var1 = pred_mod_sel$habitantes)
dev.off()

