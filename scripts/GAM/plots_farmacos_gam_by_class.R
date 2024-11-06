data_farmacos <- read.csv("Data/antidepressivos_r.csv")
data_general <- read.csv("Data/planilha_master_r.csv")

data_farmacos_reference <- data_farmacos[1,]
data_farmacos_class <- data_farmacos[2,]
data_farmacos <- data_farmacos[-c(1,2),]

data_farmacos$urb <- data_general$urb[match(data_farmacos$bacia,data_general$microbacia)]

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)

data_farmacos$urb[data_farmacos$bacia == "ebb08" | 
                    data_farmacos$bacia == "ebb09" |
                    data_farmacos$bacia == "ebb11" |
                    data_farmacos$bacia == "ebb12"] <- 0 


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
data_farmacos_reference <- unlist(c(data_farmacos_reference))

classess <- unique(data_farmacos_class)

data_farmacos_by_class <- list()

for(i in 1:length(classess)){
  sums <- rowSums(data_farmacos[,data_farmacos_class == classess[i]]) 
  data_farmacos_by_class[[i]] <- sums
  names(data_farmacos_by_class)[i] <- classess[i]
}

data_farmacos_by_class <- data.frame(bacia, Decil, urb, data_farmacos_by_class)

data_farmacos_class

LQ <- tapply(data_farmacos_reference, INDEX = data_farmacos_class, unique)


######################################################

ylab.cex <- 0.9
xlab.cex <- 0.9
axis.cex <- 0.75
line_labs <- 1.25
line_axis <- -0.75
tck <- -0.03
cex.points <- 0.7
line_lwd <- 1.25
############### SNRI ###############

SNRI_model <- select_gam(y = data_farmacos_by_class$SNRI, x = data_farmacos_by_class$urb, drop.intercept = TRUE, family = "tw")

SNRI_plot <- function(){
  plot(data_farmacos_by_class$SNRI ~ data_farmacos_by_class$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "SNRI (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(SNRI_model$estimates_best$lwr, SNRI_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$SNRI, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), SNRI_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos_by_class$urb, y = data_farmacos_by_class$SNRI, pch = 16, cex = cex.points)
  box()
}


############### SSRI ###############

SSRI_model <- select_gam(y = data_farmacos_by_class$SSRI, x = data_farmacos_by_class$urb, drop.intercept = TRUE, family = "tw")

SSRI_plot <- function(){
  plot(data_farmacos_by_class$SSRI ~ data_farmacos_by_class$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "SSRI (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(SSRI_model$estimates_best$lwr, SSRI_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$SSRI, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), SSRI_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos_by_class$urb, y = data_farmacos_by_class$SSRI, pch = 16, cex = cex.points)
  box()
  
}




############### Aminoketone ###############

Aminoketone_model <- select_gam(y = data_farmacos_by_class$Aminoketone, x = data_farmacos_by_class$urb, drop.intercept = TRUE, family = "tw")

Aminoketone_plot <- function(){
  plot(data_farmacos_by_class$Aminoketone ~ data_farmacos_by_class$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Aminoketone (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Aminoketone_model$estimates_second$lwr, Aminoketone_model$estimates_second$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Aminoketone, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Aminoketone_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos_by_class$urb, y = data_farmacos_by_class$Aminoketone, pch = 16, cex = cex.points)
  box()
  
}



############### Tricyclic ###############

Tricyclic_model <- select_gam(y = data_farmacos_by_class$Tricyclic, x = data_farmacos_by_class$urb, drop.intercept = TRUE, family = "tw")

Tricyclic_plot <- function(){
  plot(data_farmacos_by_class$Tricyclic ~ data_farmacos_by_class$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Tricyclic (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Tricyclic_model$estimates_second$lwr, Tricyclic_model$estimates_second$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Tricyclic, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Tricyclic_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos_by_class$urb, y = data_farmacos_by_class$Tricyclic, pch = 16, cex = cex.points)
  box()
}



pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/Antidepressivos_GAM_by_class.pdf", width = 4, height = 3, pointsize = 9)
par(mfrow = c(2,2), mar = c(3.5,3.5,0.5,0.5))

SSRI_plot()
SNRI_plot()
Aminoketone_plot()
Tricyclic_plot()


dev.off()

