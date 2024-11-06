data_IQVIA <- read.csv("IQVIA_antidepressants_kgs.csv")

data_IQVIA_2021 <- data_IQVIA[data_IQVIA$year == 2021,]



data_farmacos <- read.csv("Data/antidepressivos_r.csv")
#data_farmacos <- data_farmacos[,which(colnames(data_farmacos)!="Norfluoxetina")]
data_general <- read.csv("Data/planilha_master_r.csv")

data_farmacos_reference <- data_farmacos[1,]
data_farmacos_class <- data_farmacos[2,]
data_farmacos <- data_farmacos[-c(1,2),]

data_farmacos$urb <- data_general$urb[match(data_farmacos$bacia,data_general$microbacia)]

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)

data_farmacos$urb[data_farmacos$bacia == "EBB_08" | 
                    data_farmacos$bacia == "EBB_09" |
                    data_farmacos$bacia == "EBB_11" |
                    data_farmacos$bacia == "EBB_12"] <- 0 

data_farmacos$urb[is.na(data_farmacos$urb)] <- data_farmacos$Decil[is.na(data_farmacos$urb)] + 5

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)


urb <- data_farmacos$urb
Decil <- data_farmacos$Decil
bacia <- data_farmacos$bacia

data_farmacos <- data_farmacos[,-c(1,2,ncol(data_farmacos))]
data_farmacos_reference <- data_farmacos_reference[,-c(1,2)]

data_farmacos <- data.frame(lapply(data_farmacos,as.numeric))

library(vegan)

data_farmacos_summing_met <- data_farmacos

data_farmacos_summing_met$Bupropion_Hydroxybupropion <- data_farmacos_summing_met$Bupropion + data_farmacos_summing_met$Hydroxybupropion
data_farmacos_summing_met <- data_farmacos_summing_met[,colnames(data_farmacos_summing_met) != "Hydroxybupropion" &  colnames(data_farmacos_summing_met) != "Bupropion"]

data_farmacos_summing_met$Citalopram_Desmethylcitalopram <- data_farmacos_summing_met$Citalopram + data_farmacos_summing_met$Desmethylcitalopram
data_farmacos_summing_met <- data_farmacos_summing_met[,colnames(data_farmacos_summing_met) != "Citalopram" &  colnames(data_farmacos_summing_met) != "Desmethylcitalopram"]

data_farmacos_summing_met$Fluoxetine_Norfluoxetina <- data_farmacos_summing_met$Fluoxetine + data_farmacos_summing_met$Norfluoxetina
data_farmacos_summing_met <- data_farmacos_summing_met[,colnames(data_farmacos_summing_met) != "Fluoxetine" &  colnames(data_farmacos_summing_met) != "Norfluoxetina"]

#data_farmacos_summing_met$Venlafaxine_O.desmethylvenlafaxine <- data_farmacos_summing_met$Venlafaxine + data_farmacos_summing_met$O.desmethylvenlafaxine
#data_farmacos_summing_met <- data_farmacos_summing_met[,colnames(data_farmacos_summing_met) != "Venlafaxine" &  colnames(data_farmacos_summing_met) != "O.desmethylvenlafaxine"]


data_farmacos_freq <- decostand(data_farmacos_summing_met, method = "pa")

colSums(data_farmacos)
colSums(data_farmacos_freq)/nrow(data_farmacos_freq)


#data_farmacos_reference <- data_farmacos_reference[,colSums(data_farmacos)>0]

#data_farmacos <- data_farmacos[,colSums(data_farmacos) > 0]

Sums_farmacos <- colSums(data_farmacos)

Sums_farmacos <- Sums_farmacos[order(Sums_farmacos, decreasing = TRUE)]



##################################################################

data_IQVIA_vs_encontrados <- read.csv("IQVIA_vs_moleculas_encontradas.csv")

data_IQVIA_vs_encontrados

plot(data_IQVIA_vs_encontrados$ng_L ~ data_IQVIA_vs_encontrados$kg)
plot(data_IQVIA_vs_encontrados$frequencia ~ data_IQVIA_vs_encontrados$kg)

cor.test(data_IQVIA_vs_encontrados$ng_L, data_IQVIA_vs_encontrados$kg)
cor.test(data_IQVIA_vs_encontrados$frequencia, data_IQVIA_vs_encontrados$kg)


library(glmmTMB)
library(DHARMa)
library(AICcmodavg)
library(car)

mod_gaus <- glmmTMB(ng_L~kg, family = "gaussian", data = data_IQVIA_vs_encontrados)
mod_tw <- glmmTMB(ng_L~kg, family = "tweedie", data = data_IQVIA_vs_encontrados)

plot(simulateResiduals(mod_gaus))
plot(simulateResiduals(mod_tw))
AICc(mod_gaus)
AICc(mod_tw)

mod_tw0 <- glmmTMB(ng_L~1, family = "tweedie", data = data_IQVIA_vs_encontrados)

aictab(list(mod_tw0, mod_tw))


mod_gaus <- glmmTMB(frequencia~kg, family = "gaussian", data = data_IQVIA_vs_encontrados)
mod_beta <- glmmTMB(frequencia~kg, family = ordbeta(link = "probit"), data = data_IQVIA_vs_encontrados, control = glmmTMBControl(optimizer = optim))
#mod_nb <- glmmTMB(n_riachos~kg, family = "nbinom2", data = data_IQVIA_vs_encontrados)
#mod_pois <- glmmTMB(n_riachos~kg, family = "poisson", data = data_IQVIA_vs_encontrados)


plot(simulateResiduals(mod_gaus))
plot(simulateResiduals(mod_beta))

#plot(simulateResiduals(mod_pois))
#plot(simulateResiduals(mod_nb))

AICc(mod_gaus)
AICc(mod_pois)
#AICc(mod_nb)

mod_gaus0 <- glmmTMB(frequencia~1, family = "gaussian", data = data_IQVIA_vs_encontrados)

aictab(list(mod_gaus0, mod_gaus))





new_data <- data.frame(kg = seq(min(data_IQVIA_vs_encontrados$kg), max(data_IQVIA_vs_encontrados$kg), length.out = 100))

predicted_ng_L <- predict(mod_tw, newdata = new_data, type = "response")

predicted_freq <- predict(mod_gaus, newdata = new_data, type = "response")




png(paste(caminho, "Vendas_VS_Encontrado_nomes.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mfrow = c(1,1),mar = c(3.5,3.55,1,0.1), bty = "l")


plot(data_IQVIA_vs_encontrados$ng_L~data_IQVIA_vs_encontrados$kg , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(ylab = "Total Concentration (ng/L)", line = 2., cex.lab = 1.25)
title(xlab = "Sales (t)", line = 2, cex.lab = 1.25)
#title(ylab = "", line = 2, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data$kg, predicted_ng_L, lty = 1, lwd = 2)

#points(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$ng_L, pch = 1, cex = 1.5, col = "black", lwd = 2 *0.7)

labels <- substr(data_IQVIA_vs_encontrados$Antidepressant, 1, 3)

text(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$ng_L, labels = labels)
dev.off()



png(paste(caminho, "Vendas_VS_Encontrado_pontos.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mfrow = c(1,1),mar = c(3.5,3.55,1,0.1), bty = "l")


plot(data_IQVIA_vs_encontrados$ng_L~data_IQVIA_vs_encontrados$kg , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(ylab = "Total Concentration (ng/L)", line = 2., cex.lab = 1.25)
title(xlab = "Sales (t)", line = 2, cex.lab = 1.25)
#title(ylab = "", line = 2, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data$kg, predicted_ng_L, lty = 1, lwd = 2)

points(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$ng_L, pch = 1, cex = 1.5, col = "black", lwd = 2 *0.7)

#labels <- substr(data_IQVIA_vs_encontrados$Antidepressant, 1, 3)

#text(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$ng_L, labels = labels)
dev.off()





png(paste(caminho, "Vendas_VS_Encontrado_pontos_nomes.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mfrow = c(1,1),mar = c(3.5,3.55,1,0.1), bty = "l")


plot(data_IQVIA_vs_encontrados$ng_L~data_IQVIA_vs_encontrados$kg , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(ylab = "Total Concentration (ng/L)", line = 2., cex.lab = 1.25)
title(xlab = "Sales (t)", line = 2, cex.lab = 1.25)
#title(ylab = "", line = 2, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data$kg, predicted_ng_L, lty = 1, lwd = 2)

points(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$ng_L, pch = 1, cex = 1.5, col = "black", lwd = 2 *0.7)

labels <- substr(data_IQVIA_vs_encontrados$Antidepressant, 1, 3)

data_IQVIA_vs_encontrados$kg

ycoord <- data_IQVIA_vs_encontrados$ng_L
ycoord[labels == "Clo"] <- ycoord[labels == "Clo"] + 150
ycoord[labels == "Ven"] <- ycoord[labels == "Ven"] - 300


xcoord <- data_IQVIA_vs_encontrados$kg
xcoord[labels == "Par"] <- xcoord[labels == "Par"] - 200
xcoord[labels == "Mir"] <- xcoord[labels == "Mir"] + 1400

text(xcoord, ycoord+150, labels = labels)
dev.off()




png(paste(caminho, "Vendas_VS_frequencia_pontos_nomes.png", sep = ""), width = 10, height = 9, units = "cm", pointsize = 7, res = 600)

par(mfrow = c(1,1),mar = c(3.5,3.55,1,0.1), bty = "l")

plot(data_IQVIA_vs_encontrados$frequencia~data_IQVIA_vs_encontrados$kg , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

title(ylab = "Frequency of antidepressants", line = 2., cex.lab = 1.25)
title(xlab = "Sales (t)", line = 2, cex.lab = 1.25)
#title(ylab = "", line = 2, cex.lab = 1.25)

axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

#polygon(x = c(0:100,100:0), y = c(predicted_lower,predicted_upper[101:1] ), col = "grey80", border = FALSE)
lines(new_data$kg, predicted_freq, lty = 1, lwd = 2)

points(data_IQVIA_vs_encontrados$kg, data_IQVIA_vs_encontrados$frequencia, pch = 1, cex = 1.5, col = "black", lwd = 2 *0.7)

labels <- substr(data_IQVIA_vs_encontrados$Antidepressant, 1, 3)


ycoord <- data_IQVIA_vs_encontrados$frequencia
ycoord[labels == "Clo"] <- ycoord[labels == "Clo"] + 0.025
ycoord[labels == "Ven"] <- ycoord[labels == "Ven"] - 0.05


xcoord <- data_IQVIA_vs_encontrados$kg
xcoord[labels == "Par"] <- xcoord[labels == "Par"] - 200
xcoord[labels == "Mir"] <- xcoord[labels == "Mir"] + 1400

text(xcoord, ycoord+0.025, labels = labels)
dev.off()


