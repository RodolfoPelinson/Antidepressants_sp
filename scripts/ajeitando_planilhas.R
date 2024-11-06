library(vegan)
library(shape)
library(scales)
library(dplyr)

#data_farmacos <- read.csv("Data/antidepressivos_r.csv")
#data_general <- read.csv("Data/planilha_master_r.csv")
data_agua <- read.csv("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/data/water_quality_R_ultima_versao.csv")
data_topografia <- read.csv("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/data/dados_topografia.csv")
data_saneamento <- read.csv("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/data/consolidacao_dados_socioeconomicos_20240624.csv")
data_coord <- read.csv("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/data/coordinates.csv")


data_agua <- data_agua[data_agua$Survey == "seca 2021" | data_agua$Survey == "Survey" | data_agua$Survey == "class",]

data_agua_interesse <- select(data_agua, Catchment, Sampling_date_order, Caffeine, Discharge_.L.s., TOC, TC)

data_agua_interesse <- data_agua_interesse[-c(1:2),]

data_agua_interesse[,3:ncol(data_agua_interesse)] <- apply(data_agua_interesse[,3:ncol(data_agua_interesse)], 2, as.numeric)

data_cafeina <- data_agua_interesse[data_agua_interesse$Sampling_date_order == "1+2+3",1:3]

data_agua_interesse <- data_agua_interesse[data_agua_interesse$Sampling_date_order != "1+2+3",-3]

data_agua_desc_toc_tc <- aggregate(data_agua_interesse[,3:ncol(data_agua_interesse)], by = list(data_agua_interesse$Catchment), FUN = "mean", na.rm = TRUE)

data_agua[3:nrow(data_agua),data_agua[1,] == "antidepressant"] <- data.frame(apply(data_agua[3:nrow(data_agua),data_agua[1,] == "antidepressant"], 2, as.numeric))
data_agua <- data_agua[is.na(data_agua$q_O_desmethylvenlafaxine)==FALSE,]

urb <- as.numeric(data_agua$X.urb)[-c(1:2)]
index <- data_agua$Catchment[-c(1:2)]


data_farmacos <- data_agua[,data_agua[1,] == "antidepressant"]
ncol(data_farmacos)
data_farmacos <- data_farmacos[,12:22]
data_farmacos <- data_farmacos[,which(colnames(data_farmacos)!="Norfluoxetina")]


data_farmacos_class <- data_farmacos[2,]

data_farmacos <- data_farmacos[-c(1:2),]

data_farmacos <- apply(data_farmacos,2,as.numeric)

data_farmacos_summed <- rowSums(data_farmacos)


data_farmacos_class <- unlist(c(data_farmacos_class))

classess <- unique(data_farmacos_class)

classess <- classess

data_farmacos_by_class <- list()

data_farmacos <- data.frame(data_farmacos)

for(i in 1:length(classess)){
  sums <- rowSums(data_farmacos[,data_farmacos_class == classess[i]]) 
  data_farmacos_by_class[[i]] <- sums
  names(data_farmacos_by_class)[i] <- classess[i]
}

data_farmacos_by_class <- data.frame(data_farmacos_by_class)



###################################################### SANEAMENTO

length(index)
length(data_topografia$riacho)

match(data_topografia$riacho, index)
match(index, data_cafeina$Catchment)
match(index, data_agua_desc_toc_tc$Group.1)

match(index, data_topografia$riacho)



data_topografia <- data_topografia[match(index, data_topografia$riacho),]
data_cafeina <- data_cafeina[match(index, data_cafeina$Catchment),]
data_agua_desc_toc_tc <- data_agua_desc_toc_tc[match(index, data_agua_desc_toc_tc$Group.1),]
data_saneamento <- data_saneamento[match(index, data_saneamento$ID_geral),]

nrow(data_topografia)
nrow(data_cafeina)
nrow(data_agua_desc_toc_tc)

data.frame(index, data_topografia$riacho, data_agua_desc_toc_tc$Group.1, data_cafeina$Catchment, data_saneamento$ID_geral)

######################################################

# V RESPOSTAS

# Area
# Declividade
# Forma da Bacia
# Descarga
# C Organico Total
# Cafeína
# Cobertura Urbana
# N casas (densidade)
# N casas sem água (densidade)
# N casas sem esgoto e fossa (densidade)
# Renda per capta

predictors <- data.frame(
  area = data_topografia$Area_ha,#area em hecatares
  declividade = data_topografia$Declividade_av,#Declividade
  IC = data_topografia$Ic, # indice de circularidade (forma)
  discharge = data_agua_desc_toc_tc$Discharge_.L.s., #Descarga
  toc = data_agua_desc_toc_tc$TOC, # Carbono orgânico total
  #tc = data_agua_desc_toc_tc$tc, # Carbono total
  #caffeine = data_cafeina$Caffeine, # Carbono total
  imperm = data_saneamento$perc_impermeabilizada, #Cobertura urbana
  urb = data_saneamento$perc_.area.urbana_2021, #Cobertura urbana
  #dens_dom = data_saneamento$Col8_AU_2021_ha, #Densidade de domicilios
  #dens_dom_sem_agua = data_saneamento$V002.V012_domicilios_sem_agua/data_topografia$Area_ha, #Densidade de domicilios sem acesso a agua
  casas_sem_saneamento = data_saneamento$Censo2010_Dom_Sem_Esgoto, #Densidade de domicilios sem esgoto ou fossa
  renda = data_saneamento$renda, #Renda média
  habitantes = data_saneamento$habitantes,
  habitantes_sem_saneamento = data_saneamento$habitantes_sem_saneamento
)

responses <- data.frame(all = data_farmacos_summed, data_farmacos_by_class)

rownames(predictors) <- index
rownames(responses) <- index

predictors
responses



diff_watersheds_10_urb <- data_saneamento$ID_geral[data_saneamento$incremento.2010_2021 >= 10 | data_saneamento$incremento.2010_2021 <= -10]
diff_watersheds_5_urb <- data_saneamento$ID_geral[data_saneamento$incremento.2010_2021 >= 5 | data_saneamento$incremento.2010_2021 <= -5]

diff_watersheds_10_AGSB <- data_saneamento$ID_geral[data_saneamento$incremento_2010_2019 >= 10 | data_saneamento$incremento_2010_2019 <= -10]
diff_watersheds_5_AGSB <- data_saneamento$ID_geral[data_saneamento$incremento_2010_2019 >= 5 | data_saneamento$incremento_2010_2019 <= -5]



##########################################
coord <- data_coord[,1:3]

nrow(coord)
nrow(predictors)


pos <- match(rownames(predictors),coord$bacia)

faltando <- which(is.na(pos))
rownames(predictors)[faltando]

coord$bacia[pos]

coord <- coord[pos,]

nrow(coord)
