library(vegan)
library(shape)
library(scales)
library(metacom)

data_farmacos <- read.csv("Data/antidepressivos_r.csv")
data_farmacos <- data_farmacos[,which(colnames(data_farmacos)!="Norfluoxetina")]
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

urb <- data_farmacos$urb
Decil <- data_farmacos$Decil
bacia <- data_farmacos$bacia

data_farmacos <- data_farmacos[,-c(1,2,ncol(data_farmacos))]
data_farmacos_reference <- data_farmacos_reference[,-c(1,2)]
data_farmacos_class <- data_farmacos_class[,-c(1,2)]

data_farmacos <- data.frame(lapply(data_farmacos,as.numeric))

data_farmacos_reference <- data_farmacos_reference[,colSums(data_farmacos)>0]
data_farmacos_class <- data_farmacos_class[,colSums(data_farmacos) > 0]
data_farmacos <- data_farmacos[,colSums(data_farmacos) > 0]

data_farmacos_class <- unlist(c(data_farmacos_class))


data_farmacos_pa <- decostand(data_farmacos, method = "pa")

data_farmacos_pa_no_zero <- data_farmacos_pa[rowSums(data_farmacos_pa) > 0,]
data_farmacos_pa_zero <- data_farmacos_pa[rowSums(data_farmacos_pa) == 0,]

urb_no_zero <- urb[rowSums(data_farmacos_pa) > 0]
urb_zero <- urb[rowSums(data_farmacos_pa) == 0]

bacia_no_zero <- bacia[rowSums(data_farmacos_pa) > 0]
bacia_zero <- bacia[rowSums(data_farmacos_pa) == 0]


######################################################################
rownames(data_farmacos_pa_zero) <- bacia_zero
rownames(data_farmacos_pa_no_zero) <- bacia_no_zero

names(urb_no_zero) <- bacia_no_zero
names(urb_zero) <- bacia_zero



set.seed(2); order <- sample(1:nrow(data_farmacos_pa_zero))
reordered_data_farmacos_pa_zero <- data_farmacos_pa_zero[order,]
reordered_urb_zero<- urb_zero[order]


reordered_no_zero <- OrderMatrix(data_farmacos_pa_no_zero)
order_urb_no_zero <- match(rownames(reordered_no_zero),names(urb_no_zero))
reordered_urb_no_zero <- urb_no_zero[order_urb_no_zero]

reordered_data_farmacos_pa_zero <- reordered_data_farmacos_pa_zero[,match(colnames(reordered_no_zero), colnames(reordered_data_farmacos_pa_zero))]

data_farmacos_pa_plot <- rbind(reordered_data_farmacos_pa_zero, reordered_no_zero)

urb_plot <- c(reordered_urb_no_zero[order(1:length(reordered_urb_no_zero), decreasing  = TRUE)], reordered_urb_zero[order(1:length(reordered_urb_zero), decreasing  = TRUE)]) #Aqui tem que inverter a ordem com relação à matriz, não lembro bem pq. Quando a função mesmo que reordena funciona, mas como ela nao está reordenando, precisa fazer manualmente


empty_com <- rep(0, nrow(data_farmacos_pa_plot))
riqueza <- rowSums(data_farmacos_pa_plot)
empty_com[riqueza == 0] <- 1

data_farmacos_pa_plot <- data.frame(data_farmacos_pa_plot, None_detected = empty_com)

farmacos_plus_empty <-  sub("_", " ",colnames(data_farmacos_pa_plot)) 

data.frame(names(urb_plot), rownames(data_farmacos_pa_plot))

pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/EMS_empty.pdf", width = 6, height = 6)
My_Imagine(data_farmacos_pa_plot, col = c(0,1, "grey60"), top_margin = 7, left_margin = 3, bottom_margin = 0.1, fill = FALSE,
           Env1 = urb_plot, Env.col_1 = c("white","black"), Env.label_1 = "Urban cover", order = FALSE, speciesnames = farmacos_plus_empty,
           cex.site = 0.6, cex.species = 0.7, xline = -0.5, yline = -0.5, gap.axis = -10, Empty = ncol(data_farmacos_pa_plot), Empty_col = "lightblue3", speciesfont = 1)
dev.off()

pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/EMS.pdf", width = 6, height = 6)
My_Imagine(data_farmacos_pa_no_zero, col = c(0,1, "black"), top_margin = 7, left_margin = 3, bottom_margin = 0.1, fill = FALSE, cex.site = 0.6, cex.species = 0.7, xline = -0.5, yline = -0.5, gap.axis = -10,
           Env1 = urb_no_zero, Env.col_1 = c("white","black"), Env.label_1 = "Urbanização", speciesnames = colnames(data_farmacos_pa_no_zero), speciesfont = 1)
dev.off()




######################################################################

rowSums(data_farmacos_pa_no_zero)

coerencia <- data.frame(Coherence(data_farmacos_pa_no_zero, method = "r0", sims = 10000, orderNulls = TRUE, order = TRUE, seed = 1, scores = 1, allowEmpty = FALSE))
coerencia$stat <- round(coerencia$stat, 3)

turnover <- Turnover(data_farmacos_pa_no_zero, method = "r0", sims = 10000, fill = FALSE)
turnover$stat <- round(turnover$stat, 3)

bc <- BoundaryClump(data_farmacos_pa_no_zero)
bc$stat <- round(bc$stat, 3)

scores <- OrderMatrix(data_farmacos_pa_no_zero, outputScores = T)$sitescores
cor.test(scores, urb_no_zero, method = "spearman")


scores2 <- OrderMatrix(data_farmacos_pa_no_zero, outputScores = T)$sitescores
scores2 <- scores2[order(scores2, decreasing = TRUE)]
scores_zero <- rep(scores2[1]*1.1, length(bacia_zero))
names(scores_zero) <- names(reordered_urb_zero)
scores2 <- c(scores_zero, scores2)
scores2 <- scores2[length(scores2):1]

data.frame(names(scores2), names(urb_plot))

cor.test(scores2, urb_plot, method = "spearman")

