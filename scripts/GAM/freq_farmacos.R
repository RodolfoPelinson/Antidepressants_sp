data_farmacos <- read.csv("data/antidepressivos_r.csv")
data_farmacos <- data_farmacos[,which(colnames(data_farmacos)!="Norfluoxetina")]
data_general <- read.csv("data/planilha_master_r.csv")

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

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)


urb <- data_farmacos$urb
Decil <- data_farmacos$Decil
bacia <- data_farmacos$bacia

data_farmacos <- data_farmacos[,-c(1,2,ncol(data_farmacos))]
data_farmacos_reference <- data_farmacos_reference[,-c(1,2)]
LQ <- data_farmacos_reference

data_farmacos <- data.frame(lapply(data_farmacos,as.numeric))

data_farmacos <- data_farmacos[,colSums(data_farmacos) > 0]

library(vegan)

data_farmacos_pa <- decostand(data_farmacos, method = "pa")

freq <- function(x){
  freq <- sum(x)/length(x)
  return(freq)
}

frequencies <- list()

for(i in 1:ncol(data_farmacos_pa)){
  frequencies[[i]] <- tapply(data_farmacos_pa[,i], INDEX = Decil, freq)
  names(frequencies)[i] <- colnames(data_farmacos_pa)[i]
}

frequencies <- lapply(frequencies, as.numeric)

Decil_unique <- unique(Decil)


library(bbmle)
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/select_gam.R")

ylab.cex <- 0.9
xlab.cex <- 0.9
axis.cex <- 0.75
line_labs <- 1.25
line_axis <- -0.75
tck <- -0.03
cex.points <- 0.7
line_lwd <- 1.25

############### Sertraline ###############

Sertraline_model <- select_gam(frequencies$Sertraline, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Sertraline_plot <- function(letter = ""){
  plot(frequencies$Sertraline ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "sertraline", line = 1.25, cex.lab = ylab.cex)

  polygon(y = c(Sertraline_model$estimates_best$lwr, Sertraline_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Sertraline, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Sertraline_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Sertraline, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



############### O.desmethylvenlafaxine ###############

O.desmethylvenlafaxine_model <- select_gam(frequencies$O.desmethylvenlafaxine, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

O.desmethylvenlafaxine_plot <- function(letter = ""){
  plot(frequencies$O.desmethylvenlafaxine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "O-desmethylvenlafaxine", line = 1.25, cex.lab = ylab.cex)
  
  #polygon(y = c(O.desmethylvenlafaxine_model$estimates_best$lwr, O.desmethylvenlafaxine_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  #polygon(y = c(O.desmethylvenlafaxine_model$estimates_second$lwr, O.desmethylvenlafaxine_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  
  abline(h = LQ$O.desmethylvenlafaxine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), O.desmethylvenlafaxine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  #lines(c(0:100), O.desmethylvenlafaxine_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  
  points(x = Decil_unique, y = frequencies$O.desmethylvenlafaxine, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}


############### Hydroxybupropion ###############

Hydroxybupropion_model <- select_gam(frequencies$Hydroxybupropion, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Hydroxybupropion_plot <- function(letter = ""){
  plot(frequencies$Hydroxybupropion ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Hydroxybupropion", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Hydroxybupropion_model$estimates_best$lwr, Hydroxybupropion_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Hydroxybupropion, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Hydroxybupropion_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Hydroxybupropion, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}




############### Bupropion ###############

Bupropion_model <- select_gam(frequencies$Bupropion, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Bupropion_plot <- function(letter = ""){
  plot(frequencies$Bupropion ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Bupropion", line = 1.25, cex.lab = ylab.cex)
  
  #polygon(y = c(Bupropion_model$estimates_best$lwr, Bupropion_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Bupropion, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Bupropion_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Bupropion, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



############### Venlafaxine ###############

Venlafaxine_model <- select_gam(frequencies$Venlafaxine, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Venlafaxine_plot <- function(letter = ""){
  plot(frequencies$Venlafaxine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Venlafaxine", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Venlafaxine_model$estimates_best$lwr, Venlafaxine_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Venlafaxine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Venlafaxine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Venlafaxine, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}




############### Desmethylcitalopram ###############

Desmethylcitalopram_model <- select_gam(frequencies$Desmethylcitalopram, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Desmethylcitalopram_plot <- function(letter = ""){
  plot(frequencies$Desmethylcitalopram ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Desmethylcitalopram", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Desmethylcitalopram_model$estimates_best$lwr, Desmethylcitalopram_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Desmethylcitalopram, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Desmethylcitalopram_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Desmethylcitalopram, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}




############### Citalopram ###############

Citalopram_model <- select_gam(frequencies$Citalopram, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Citalopram_plot <- function(letter = ""){
  plot(frequencies$Citalopram ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Citalopram", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Citalopram_model$estimates_best$lwr, Citalopram_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Citalopram, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Citalopram_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Citalopram, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



############### Norfluoxetina ###############

#Norfluoxetina_model <- select_gam(frequencies$Norfluoxetina, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

#Norfluoxetina_plot <- function(){
#  plot(frequencies$Norfluoxetina ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
#  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
#  
#  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
#  title(ylab = "Norfluoxetine", line = 1.25, cex.lab = ylab.cex)
  
  #polygon(y = c(Norfluoxetina_model$estimates_best$lwr, Norfluoxetina_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
#  abline(h = LQ$Norfluoxetina, lty = 2, col = "brown")
#  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
#  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Norfluoxetina_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
#  points(x = Decil_unique, y = frequencies$Norfluoxetina, pch = 16, cex = cex.points)
#  box()
#}



############### Fluoxetine ###############

Fluoxetine_model <- select_gam(frequencies$Fluoxetine, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Fluoxetine_plot <- function(letter = ""){
  plot(frequencies$Fluoxetine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Fluoxetine", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Fluoxetine_model$estimates_best$lwr, Fluoxetine_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Fluoxetine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Fluoxetine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Fluoxetine, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



############### Amitriptyline ###############

Amitriptyline_model <- select_gam(frequencies$Amitriptyline, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Amitriptyline_plot <- function(letter = ""){
  plot(frequencies$Amitriptyline ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Amitriptyline", line = 1.25, cex.lab = ylab.cex)
  
  polygon(y = c(Amitriptyline_model$estimates_best$lwr, Amitriptyline_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Amitriptyline, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Amitriptyline_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Amitriptyline, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



############### Trimipramine ###############

Trimipramine_model <- select_gam(frequencies$Trimipramine, x =Decil_unique,  family = "betar", drop.intercept = FALSE)

Trimipramine_plot <- function(letter = ""){
  plot(frequencies$Trimipramine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)

  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Trimipramine", line = 1.25, cex.lab = ylab.cex)
  
  #polygon(y = c(Trimipramine_model$estimates_best$lwr, Trimipramine_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = LQ$Trimipramine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Trimipramine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Trimipramine, pch = 16, cex = cex.points)
  margins <- par()$mar
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = cex.letter, adj = 0)
  par(mar = margins)
  box()
}



pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/Antidepressivos_GAM_frequency.pdf", width = 4, height = 7.5)
par(mfrow = c(5,2), mar = c(2.5,3.5,1.2,0.25))

Venlafaxine_plot(letter = "A)")
O.desmethylvenlafaxine_plot(letter = "B)")
Citalopram_plot(letter = "C)")
Desmethylcitalopram_plot(letter = "D)")
Bupropion_plot(letter = "E)")
Hydroxybupropion_plot(letter = "F)")
Amitriptyline_plot(letter = "G)")
Sertraline_plot(letter = "H)")
Fluoxetine_plot(letter = "I)")
Trimipramine_plot(letter = "J)")
#Norfluoxetina_plot()

dev.off()



model_tabs <- list(
  O.desmethylvenlafaxine_model$AICTAB,
  Desmethylcitalopram_model$AICTAB,
  Hydroxybupropion_model$AICTAB,
  Amitriptyline_model$AICTAB,
  Sertraline_model$AICTAB,
  Venlafaxine_model$AICTAB,
  Fluoxetine_model$AICTAB,
  Citalopram_model$AICTAB,
  Bupropion_model$AICTAB,
  Trimipramine_model$AICTAB)

names(model_tabs) <- c("O.desmethylvenlafaxine",
                       "Desmethylcitalopram",
                       "Hydroxybupropion",
                       "Amitriptyline",
                       "Sertraline",
                       "Venlafaxine",
                       "Fluoxetine",
                       "Citalopram",
                       "Bupropion",
                       "Trimipramine")


order <- c("gam_0","gam_1","gam_2","gam_3","gam_4")

for(i in 1:length(model_tabs)){
  model_tabs[[i]] <- model_tabs[[i]][match(order,rownames(model_tabs[[i]])),]
}



write.csv(model_tabs,"Resultados/model_selection_GAM_FREQ.csv")


