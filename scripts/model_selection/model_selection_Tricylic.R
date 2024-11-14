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





###############################################
###############################################
####### TODOS OS ANTIDEPRESSIVOS ##############

Preditoras_stand$Tricyclic <- responses$Tricyclic





#Modelol global
mod_global_Tricyclic <- glmmTMB(Tricyclic ~ 
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

dredge_modelos_Tricyclic <- MuMIn::dredge(mod_global_Tricyclic, subset = smat)

models <- get.models(dredge_modelos_Tricyclic, subset = delta < 2)

aictab(models)

write.table(data.frame(as.matrix(dredge_modelos_Tricyclic)), "Resultados/model_selection_Tricyclic.csv", sep = ",")



