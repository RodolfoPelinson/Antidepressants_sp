basico_sp1 <- read.csv("data/seleção para checagem de padroes RMSP/Basico_SP1_2.csv")
basico_sp2 <- read.csv("data/seleção para checagem de padroes RMSP/Basico_SP2_2.csv")

domicilio01_SP1 <- read.csv("data/seleção para checagem de padroes RMSP/Domicilio01_SP1_2.csv")
domicilio01_SP2 <- read.csv("data/seleção para checagem de padroes RMSP/Domicilio01_SP2_2.csv")

domicilio02_SP1 <- read.csv("data/seleção para checagem de padroes RMSP/Domicilio02_SP1_2.csv")
domicilio02_SP2 <- read.csv("data/seleção para checagem de padroes RMSP/Domicilio02_SP2_2.csv")

domicilio_renda_SP1 <- read.csv("data/seleção para checagem de padroes RMSP/DomicilioRenda_SP1_2.csv")
domicilio_renda_SP2 <- read.csv("data/seleção para checagem de padroes RMSP/DomicilioRenda_SP2_2.csv")


#EM domicilio01_SP1
#Cod_setor  será a articulação com as demais planilhas

#e se for em 
#Nome_da_Meso e selecionar "Metropolitana de São Paulo" JUNTO com 
#Nome_da_RM e selecionar "RM São Paulo"

nrow(basico_sp1)
#nrow(domicilio01_SP1)
nrow(domicilio02_SP1)

nrow(basico_sp2)
#nrow(domicilio01_SP2)
nrow(domicilio02_SP2)

nrow(domicilio_renda_SP2)

#terá os 38 municípios da Região Metropolitana de São Paulo MENOS a cidade de São Paulo, que terá que vir da outra planilha

ncol(domicilio01_SP1)
ncol(domicilio01_SP2)

cod_setor_sp <- domicilio01_SP1$Cod_setor
cod_setor_restante <- domicilio01_SP2$Cod_setor[]



setores_RM_SP <- which(basico_sp2$Nome_da_meso == "Metropolitana de Sao Paulo" & basico_sp2$Nome_da_RM == "RM Sao Paulo")

moradores_SP <- domicilio02_SP1$V002
moradores_RM_SP <- as.numeric(domicilio02_SP2$V002[setores_RM_SP])

#domicilio02_SP2$V002[setores_RM_SP][which(is.na(moradores_RM_SP))]

renda_SP <- domicilio_renda_SP1$V003
renda_RM_SP <- as.numeric(domicilio_renda_SP2$V003[setores_RM_SP])

moradores <- c(moradores_SP, moradores_RM_SP)
renda <- c(renda_SP, renda_RM_SP)


renda_per_capta <- na.omit(renda/moradores)

renda_per_capta <- renda_per_capta/1.7485

min(renda_per_capta)
max(renda_per_capta)

caminho <- "C:/Users/rodol/OneDrive/repos/Antidepressants_sp/Resultados/Resultados Antidepressivos/Figuras/PDF/"


pdf(paste(caminho, "HISTOGRAMA.pdf", sep = ""), width = 4, height = 3.5, pointsize = 7)
par(mar = c(4,4,0.1,0.1), cex.lab = 1.2)
hist(renda_per_capta, breaks = seq(from = 0, to = 26500/1.7485, by = 250/1.7485), xlim = c(0,5000), xlab = "Average Income (US$)", main = "", lwd = 0.5)
abline(v = min(predictors$renda)/1.7485, lty = 2, col = "red", lwd = 1)
abline(v = max(predictors$renda)/1.7485, lty = 2, col = "red", lwd = 1)
box(bty = "l")

dev.off()


total_de_setores <- length(renda_per_capta)
setores_cobertos <- length(renda_per_capta[renda_per_capta >= min(predictors$renda) & renda_per_capta <= max(predictors$renda)])
perc_setores_cobertos <- setores_cobertos/total_de_setores
perc_setores_cobertos




hist(renda_per_capta, breaks = seq(from = 0, to = 26500/1.7485, by = 250/1.7485), xlim = c(0,27000/1.7485), ylim = c(0,10), xlab = "Average Income (R$)", main = "")

