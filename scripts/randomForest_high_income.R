#pred_mod_sel <- select(predictors, area, casas_sem_saneamento, renda, habitantes)
#pred_mod_sel <- select(predictors,  casas_sem_saneamento, renda, habitantes)
library(vegan)

pred_mod_sel_low_income <- pred_mod_sel[pred_mod_sel$renda >= median(pred_mod_sel$renda),]
all <- responses$all[pred_mod_sel$renda >= median(pred_mod_sel$renda)]
nrow(pred_mod_sel_low_income)
#Padronizando
Preditoras_stand <- decostand(pred_mod_sel_low_income, method = "stand")


##################################################################################################
library(randomForest)

RFmodel <- randomForest(all ~ casas_sem_saneamento + renda + habitantes,
                        data = Preditoras_stand, proximity = TRUE, ntree = 100000,
                        mtry = 1, importance = TRUE)
RFmodel

varImpPlot(RFmodel)

par(mfrow = c(3,1))
partialPlot(RFmodel, Preditoras_stand, "casas_sem_saneamento")
partialPlot(RFmodel, Preditoras_stand, "renda")
partialPlot(RFmodel, Preditoras_stand, "habitantes")
###################################################################################################