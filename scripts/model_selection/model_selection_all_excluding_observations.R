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

new_predictors <- predictors

predictors_sem_10_urb <- predictors[-match(diff_watersheds_10_urb, rownames(new_predictors)),]
predictors_sem_5_urb <- predictors[-match(diff_watersheds_5_urb, rownames(new_predictors)),]

predictors_sem_10_AGSB <- predictors[-match(diff_watersheds_10_AGSB, rownames(new_predictors)),]
predictors_sem_5_AGSB <- predictors[-match(diff_watersheds_5_AGSB, rownames(new_predictors)),]

nrow(new_predictors)
nrow(new_predictors) - nrow(predictors_sem_10_urb)
nrow(new_predictors) - nrow(predictors_sem_5_urb)
nrow(new_predictors) - nrow(predictors_sem_10_AGSB)
nrow(new_predictors) - nrow(predictors_sem_5_AGSB)

new_predictors <- predictors_sem_10_urb

colnames(new_predictors)
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

new_predictors_cor <- new_predictors[,-ncol(new_predictors)]
colnames(new_predictors_cor) <- pred_names
cors <- cor(new_predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(new_predictors_cor, use = "pairwise.complete.obs")

plot(new_predictors$habitantes ~ new_predictors$habitantes_sem_saneamento)
#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"

#png(paste(caminho, "Correlation_among_variables.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_significat ones.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower", sig.level = 0.05, p.mat = cors2$p, insig = "blank")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_greater_than_0.5.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors0, method = "circle", type = "lower")
#dev.off()





#pred_mod_sel <- select(new_predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(new_predictors,  casas_sem_saneamento, renda, habitantes)

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

all <- responses$all[-match(diff_watersheds_10_urb, rownames(responses))]

Preditoras_stand$all <- all

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

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all_10_urb.csv", sep = ",")




##############################################################################################
##############################################################################################
##############################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
##############################################################################################



new_predictors <- predictors_sem_5_urb

colnames(new_predictors)
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

new_predictors_cor <- new_predictors[,-ncol(new_predictors)]
colnames(new_predictors_cor) <- pred_names
cors <- cor(new_predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(new_predictors_cor, use = "pairwise.complete.obs")

plot(new_predictors$habitantes ~ new_predictors$habitantes_sem_saneamento)
#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"

#png(paste(caminho, "Correlation_among_variables.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_significat ones.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower", sig.level = 0.05, p.mat = cors2$p, insig = "blank")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_greater_than_0.5.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors0, method = "circle", type = "lower")
#dev.off()





#pred_mod_sel <- select(new_predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(new_predictors,  casas_sem_saneamento, renda, habitantes)

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

all <- responses$all[-match(diff_watersheds_5_urb, rownames(responses))]

Preditoras_stand$all <- all

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

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all_5_urb.csv", sep = ",")






##############################################################################################
##############################################################################################
##############################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
##############################################################################################



new_predictors <- predictors_sem_10_AGSB

colnames(new_predictors)
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

new_predictors_cor <- new_predictors[,-ncol(new_predictors)]
colnames(new_predictors_cor) <- pred_names
cors <- cor(new_predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(new_predictors_cor, use = "pairwise.complete.obs")

plot(new_predictors$habitantes ~ new_predictors$habitantes_sem_saneamento)
#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"

#png(paste(caminho, "Correlation_among_variables.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_significat ones.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower", sig.level = 0.05, p.mat = cors2$p, insig = "blank")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_greater_than_0.5.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors0, method = "circle", type = "lower")
#dev.off()





#pred_mod_sel <- select(new_predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(new_predictors,  casas_sem_saneamento, renda, habitantes)

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

all <- responses$all[-match(diff_watersheds_10_AGSB, rownames(responses))]

Preditoras_stand$all <- all

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

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all_10_AGSB.csv", sep = ",")





##############################################################################################
##############################################################################################
##############################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
##############################################################################################



new_predictors <- predictors_sem_5_AGSB

colnames(new_predictors)
pred_names <- c("Area", "Slope", "CI", "Discharge", "TOC", "% Impermeabilization", "% Urban cover", "HWNS", "Average Income", "Inhabitants")

new_predictors_cor <- new_predictors[,-ncol(new_predictors)]
colnames(new_predictors_cor) <- pred_names
cors <- cor(new_predictors_cor, use = "pairwise.complete.obs")
cors2 <- cor.mtest(new_predictors_cor, use = "pairwise.complete.obs")

plot(new_predictors$habitantes ~ new_predictors$habitantes_sem_saneamento)
#Pedindo para mostrar, em outra planilha, só as variaveis com correlação maior que 0.5 ou -0.5
cors0 <- cors
cors0[which(cors < 0.5 & cors > -0.5)] <- 0 


caminho <- "C:/Users/rodol/OneDrive/Trabalho/Projetos/Colaborações/Manuscrito fármacos/Análises_antidepressivos/Resultados/Resultados Antidepressivos/Figuras/PNG/"

#png(paste(caminho, "Correlation_among_variables.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_significat ones.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors, method = "circle", type = "lower", sig.level = 0.05, p.mat = cors2$p, insig = "blank")
#dev.off()

#png(paste(caminho, "Correlation_among_variables_greater_than_0.5.png", sep= ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)
corrplot(cors0, method = "circle", type = "lower")
#dev.off()





#pred_mod_sel <- select(new_predictors, area, casas_sem_saneamento, renda, habitantes)
pred_mod_sel <- select(new_predictors,  casas_sem_saneamento, renda, habitantes)

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

all <- responses$all[-match(diff_watersheds_5_AGSB, rownames(responses))]

Preditoras_stand$all <- all

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

write.table(data.frame(as.matrix(dredge_modelos_all)), "Resultados/model_selection_all_5_AGSB.csv", sep = ",")







##############################################################################################
##############################################################################################
##############################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
##############################################################################################



test <- function(n, n_remove){
  
  n_no_effect <- rep(0, n)
  
  for(i in 1:n){
    
    remove <- sample(1:nrow(predictors), n_remove)
    
    new_predictors <- predictors[-remove,]
    all <- responses$all[-remove]
    
    
    #pred_mod_sel <- select(new_predictors, area, casas_sem_saneamento, renda, habitantes)
    pred_mod_sel <- select(new_predictors,  casas_sem_saneamento, renda, habitantes)
    
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
    
    
    Preditoras_stand$all <- all
    
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
    
    
    if(any((names(models) == "1")) == TRUE){
      n_no_effect[i] <- 1 
    }
    message(paste("run",i))
    
    
  }
  percentage <- sum(n_no_effect)/length(n_no_effect)
  return(percentage)
}


percentage_1 <- test(n <- 100, n_remove <- 1)
percentage_4 <- test(n <- 100, n_remove <- 4)
percentage_9 <- test(n <- 100, n_remove <- 9)
percentage_11 <- test(n <- 100, n_remove <- 11)
percentage_14 <- test(n <- 100, n_remove <- 14)





