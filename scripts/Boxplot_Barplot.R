data_farmacos <- read.csv("antidepressivos_r.csv")
data_farmacos <- data_farmacos[,which(colnames(data_farmacos)!="Norfluoxetina")]
data_general <- read.csv("planilha_master_r.csv")

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

colSums(data_farmacos)

data_farmacos_reference <- data_farmacos_reference[,colSums(data_farmacos)>0]

data_farmacos <- data_farmacos[,colSums(data_farmacos) > 0]

Sums_farmacos <- colSums(data_farmacos)

Sums_farmacos <- Sums_farmacos[order(Sums_farmacos, decreasing = TRUE)]




library(tidyr)

data_farmacos <- data_farmacos[,order(colSums(data_farmacos), decreasing = TRUE)]

data_farmacos <- data.frame(bacia, urb, Decil, data_farmacos)


data_long <- gather(data_farmacos, key = "API", value = "measurement", O.desmethylvenlafaxine:Trimipramine, factor_key=TRUE)

boxplot(measurement ~ API, data = data_long)


labels <- as.character(unique(data_long$API))
labels[1] <- "O-desmethylvenlafaxine"


pdf("Boxplot.pdf", width = 4, height = 3.5, pointsize = 10)

par(mar = c(8.5,2.5,0.75,0.5))

boxplot(measurement ~ API, data = data_long , xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
title(ylab = "Observed concentrations (ng/L)", line = line_labs, cex.lab = ylab.cex)

axis(1, at = c(1:length(labels)), tick = TRUE, line = 0, labels = FALSE,  tck = tck)
axis(1, at = c(1:length(labels)), tick = FALSE, cex.axis = 0.8, line = 0, gap.axis = -10, labels = labels, las = 2)
axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = 0.75, line = -0.6, gap.axis = -10)

box()

dev.off()


pdf("Barplot.pdf", width = 4, height = 3.5, pointsize = 10)

par(mar = c(8.5,2.5,0.75,0.5))

barplot(Sums_farmacos, xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16, ylim = c(0,4000), space = 0.25, width = 1)
title(ylab = "Total mass of API (ng/L)", line = line_labs, cex.lab = ylab.cex)

axis(1, at = c(0.75, 0.75+(1.25*1),
               0.75+(1.25*2),
               0.75+(1.25*3),
               0.75+(1.25*4),
               0.75+(1.25*5),
               0.75+(1.25*6),
               0.75+(1.25*7),
               0.75+(1.25*8),
               0.75+(1.25*9)), tick = TRUE, line = 0, labels = FALSE,  tck = tck)
axis(1, at = c(0.75, 0.75+(1.25*1),
                    0.75+(1.25*2),
                    0.75+(1.25*3),
                    0.75+(1.25*4),
                    0.75+(1.25*5),
                    0.75+(1.25*6),
                    0.75+(1.25*7),
                    0.75+(1.25*8),
                    0.75+(1.25*9)), tick = FALSE, cex.axis = 0.8, line = 0, gap.axis = -10, labels = labels, las = 2)
axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = 0.75, line = -0.6, gap.axis = -10)

box()

dev.off()