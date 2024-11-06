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

data_farmacos_pa <- decostand(data_farmacos, method = "pa")


#data_farmacos_reference <- data_farmacos_reference[,colSums(data_farmacos)>0]
#data_farmacos_class <- data_farmacos_class[,colSums(data_farmacos) > 0]
#data_farmacos <- data_farmacos[,colSums(data_farmacos) > 0]

data_farmacos_class <- unlist(c(data_farmacos_class))

classess <- unique(data_farmacos_class)

data_farmacos_by_class <- list()

for(i in 1:length(classess)){
  if(is.null(ncol(data_farmacos_pa[,data_farmacos_class == classess[i]])) == FALSE){
    sums <- rowSums(data_farmacos_pa[,data_farmacos_class == classess[i]]) 
    data_farmacos_by_class[[i]] <- sums
    names(data_farmacos_by_class)[i] <- classess[i]
  }else{
    sums <- rep(0, nrow(data_farmacos_pa))
    data_farmacos_by_class[[i]] <- sums
    names(data_farmacos_by_class)[i] <- classess[i]
  }
}

data_farmacos_by_class <- data.frame(data_farmacos_by_class)

n_pharma <- rowSums(data_farmacos_by_class)

data_farmacos <- data.frame(bacia, Decil, urb, data_farmacos_by_class, n_pharma)


ylab.cex <- 0.9
xlab.cex <- 0.9
axis.cex <- 0.75
line_labs <- 1.25
line_axis <- -0.75
tck <- -0.03
cex.points <- 0.7
line_lwd <- 1.25

############### N Farmacos ###############

n_pharma_model <- select_gam(y = data_farmacos$n_pharma, x = data_farmacos$urb, family = "nb", drop.intercept = TRUE)

n_pharma_plot <- function(){
  plot(data_farmacos$n_pharma ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Number of API", line = line_labs, cex.lab = ylab.cex)
  polygon(y = c(n_pharma_model$estimates_best$lwr, n_pharma_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$n_pharma, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), n_pharma_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$n_pharma, pch = 16, cex = cex.points)
  box()
}


############### SNRI ###############

SNRI_model <- select_gam(y = data_farmacos$SNRI, x = data_farmacos$urb, family = "nb", drop.intercept = TRUE)

SNRI_plot <- function(){
  plot(data_farmacos$SNRI ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,2.5))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Number of API (SNRI)", line = line_labs, cex.lab = ylab.cex)
  polygon(y = c(SNRI_model$estimates_best$lwr, SNRI_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$SNRI, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, at = c(0,1,2),line = 0, labels = FALSE,  tck = tck); 
  axis(2, tick = FALSE, at = c(0,1,2), cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), SNRI_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$SNRI, pch = 16, cex = cex.points)
  box()
}



############### SSRI ###############

SSRI_model <- select_gam(y = data_farmacos$SSRI, x = data_farmacos$urb, family = "nb", drop.intercept = TRUE)

SSRI_plot <- function(){
  plot(data_farmacos$SSRI ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Number of API (SSRI)", line = line_labs, cex.lab = ylab.cex)
  polygon(y = c(SSRI_model$estimates_best$lwr, SSRI_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$SSRI, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), SSRI_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$SSRI, pch = 16, cex = cex.points)
  box()
}



############### Tricyclic ###############

Tricyclic_model <- select_gam(y = data_farmacos$Tricyclic, x = data_farmacos$urb, family = "nb", drop.intercept = TRUE)

Tricyclic_plot <- function(){
  plot(data_farmacos$Tricyclic ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Number of API (Tricyclic)", line = line_labs, cex.lab = ylab.cex)
  #polygon(y = c(Tricyclic_model$estimates_second$lwr, Tricyclic_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  #polygon(y = c(Tricyclic_model$estimates_best$lwr, Tricyclic_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  abline(h = data_farmacos_reference$Tricyclic, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0,labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0,  at = c(0,1,2), labels = FALSE,  tck = tck); axis(2,  at = c(0,1,2),tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Tricyclic_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  #lines(c(0:100), Tricyclic_model$estimates_best$fit, lty = 1, lwd  = line_lwd, col = "black")
  
  points(x = data_farmacos$urb, y = data_farmacos$Tricyclic, pch = 16, cex = cex.points)
  box()
}



############### Aminoketone ###############

Aminoketone_model <- select_gam(y = data_farmacos$Aminoketone, x = data_farmacos$urb, family = "nb", drop.intercept = TRUE)

Aminoketone_plot <- function(){
  plot(data_farmacos$Aminoketone ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Number of API", line = 2, cex.lab = ylab.cex)
  title(ylab = "(Aminoketone)", line = 1.25, cex.lab = ylab.cex)
  #polygon(y = c(Aminoketone_model$estimates_second$lwr, Aminoketone_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  #polygon(y = c(Aminoketone_model$estimates_best$lwr, Aminoketone_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  abline(h = data_farmacos_reference$Aminoketone, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, at = c(0,1,2), line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, at = c(0,1,2), cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Aminoketone_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  #lines(c(0:100), Aminoketone_model$estimates_best$fit, lty = 1, lwd  = line_lwd, col = "black")
  
  points(x = data_farmacos$urb, y = data_farmacos$Aminoketone, pch = 16, cex = cex.points)
  box()
}




pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/Antidepressivos_GAM_N.pdf", width = 6, height = 3)
par(mfrow = c(2,3), mar = c(3.5,3.5,0.5,0.5))

n_pharma_plot()
SSRI_plot()
SNRI_plot()
Aminoketone_plot()
Tricyclic_plot()

dev.off()
