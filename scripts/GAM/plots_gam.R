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

data_farmacos <- data.frame(bacia, urb, Decil, data_farmacos)


ylab.cex <- 0.9
xlab.cex <- 0.9
axis.cex <- 0.75
line_labs <- 1.25
line_axis <- -0.75
tck <- -0.03
cex.points <- 0.7
line_lwd <- 1.25
############### O-desmethylvenlafaxine ###############

O.desmethylvenlafaxine_model <- select_gam(y = data_farmacos$O.desmethylvenlafaxine, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

O.desmethylvenlafaxine_plot <- function(){
  plot(data_farmacos$O.desmethylvenlafaxine ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "O-desmethylvenlafaxine", line = 2, cex.lab = ylab.cex)
  title(ylab = "(ng/L)", line = 1.25, cex.lab = ylab.cex)
  polygon(y = c(O.desmethylvenlafaxine_model$estimates_best$lwr, O.desmethylvenlafaxine_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$O.desmethylvenlafaxine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), O.desmethylvenlafaxine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$O.desmethylvenlafaxine, pch = 16, cex = cex.points)
  box()
}


############### Hydroxybupropion ###############

Hydroxybupropion_model <- select_gam(y = data_farmacos$Hydroxybupropion, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Hydroxybupropion_plot <- function(){
  plot(data_farmacos$Hydroxybupropion ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "O-desmethylvenlafaxine", line = 2, cex.lab = ylab.cex)
  title(ylab = "(ng/L)", line = 1.25, cex.lab = ylab.cex)
  polygon(y = c(Hydroxybupropion_model$estimates_best$lwr, Hydroxybupropion_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Hydroxybupropion, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Hydroxybupropion_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$Hydroxybupropion, pch = 16, cex = cex.points)
  box()
}



############### Bupropion ###############

Bupropion_model <- select_gam(y = data_farmacos$Bupropion, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Bupropion_plot <- function(){
  plot(data_farmacos$Bupropion ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Bupropion (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  #polygon(y = c(Bupropion_model$estimates_best$lwr, Bupropion_model$estimates_best$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Bupropion, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Bupropion_model$estimates_best$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  points(x = data_farmacos$urb, y = data_farmacos$Bupropion, pch = 16, cex = cex.points)
}




############### Venlafaxine ###############

Venlafaxine_model <- select_gam(y = data_farmacos$Venlafaxine, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Venlafaxine_plot <- function(){
  plot(data_farmacos$Venlafaxine ~ data_farmacos$urb , type = "n", ylim = c(0,60), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Venlafaxine (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  #polygon(y = c(Venlafaxine_model$estimates_second$lwr, Venlafaxine_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Venlafaxine, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Venlafaxine_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  points(x = data_farmacos$urb, y = data_farmacos$Venlafaxine, pch = 16, cex = cex.points)
}



############### Desmethylcitalopram ###############

Desmethylcitalopram_model <- select_gam(y = data_farmacos$Desmethylcitalopram, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Desmethylcitalopram_plot <- function(){
  plot(data_farmacos$Desmethylcitalopram ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Desmethylcitalopram", line = 2, cex.lab = ylab.cex)
  title(ylab = "(ng/L)", line = 1.25, cex.lab = ylab.cex)
  #polygon(y = c(Desmethylcitalopram_model$estimates_second$lwr, Desmethylcitalopram_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Desmethylcitalopram, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Desmethylcitalopram_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  points(x = data_farmacos$urb, y = data_farmacos$Desmethylcitalopram, pch = 16, cex = cex.points)
  box()
}



############### Citalopram ###############

Citalopram_model <- select_gam(y = data_farmacos$Citalopram, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Citalopram_plot <- function(){
  plot(data_farmacos$Citalopram ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Citalopram (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Citalopram_model$estimates_best$lwr, Citalopram_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Citalopram, lty = 2, col = "brown")
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Citalopram_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$Citalopram, pch = 16, cex = cex.points)
}





############### Norfluoxetina ###############

Norfluoxetina_model <- select_gam(y = data_farmacos$Norfluoxetina, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Norfluoxetina_plot <- function(){
  plot(data_farmacos$Norfluoxetina ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Norfluoxetina (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  #polygon(y = c(Norfluoxetina_model$estimates_second$lwr, Norfluoxetina_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Norfluoxetina, lty = 2, col = "brown") #LQ Talvez esteja errado
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  #lines(c(0:100), Norfluoxetina_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  points(x = data_farmacos$urb, y = data_farmacos$Norfluoxetina, pch = 16, cex = cex.points)
}



############### Fluoxetine ###############

Fluoxetine_model <- select_gam(y = data_farmacos$Fluoxetine, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Fluoxetine_plot <- function(){
  plot(data_farmacos$Fluoxetine ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Fluoxetine (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Fluoxetine_model$estimates_best$lwr, Fluoxetine_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Fluoxetine, lty = 2, col = "brown") #LQ Talvez esteja errado
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Fluoxetine_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$Fluoxetine, pch = 16, cex = cex.points)
}




############### Amitriptyline ###############

Amitriptyline_model <- select_gam(y = data_farmacos$Amitriptyline, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Amitriptyline_plot <- function(){
  plot(data_farmacos$Amitriptyline ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Amitriptyline (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Amitriptyline_model$estimates_best$lwr, Amitriptyline_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Amitriptyline, lty = 2, col = "brown") #LQ Talvez esteja errado
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Amitriptyline_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$Amitriptyline, pch = 16, cex = cex.points)
}





############### Trimipramine ###############

Trimipramine_model <- select_gam(y = data_farmacos$Trimipramine, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Trimipramine_plot <- function(){
  plot(data_farmacos$Trimipramine ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Trimipramine (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  #polygon(y = c(Trimipramine_model$estimates_second$lwr, Trimipramine_model$estimates_second$upr[101:1]),
  #        x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Trimipramine, lty = 2, col = "brown") #LQ Talvez esteja errado
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = -0.5, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = -0.5, gap.axis = -10)
  #lines(c(0:100), Trimipramine_model$estimates_second$fit, lty = 2, lwd  = line_lwd, col = "grey40")
  points(x = data_farmacos$urb, y = data_farmacos$Trimipramine, pch = 16, cex = cex.points)
}




############### Sertraline ###############

Sertraline_model <- select_gam(y = data_farmacos$Sertraline, x = data_farmacos$urb, drop.intercept = TRUE, family = "tw")

Sertraline_plot <- function(){
  plot(data_farmacos$Sertraline ~ data_farmacos$urb , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Sertraline (ng/L)", line = line_labs, cex.lab = ylab.cex)
  
  polygon(y = c(Sertraline_model$estimates_best$lwr, Sertraline_model$estimates_best$upr[101:1]),
          x = c(0:100,100:0), col = "grey90", border = FALSE)
  abline(h = data_farmacos_reference$Sertraline, lty = 2, col = "brown") #LQ Talvez esteja errado
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = -0.5, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = -0.5, gap.axis = -10)
  lines(c(0:100), Sertraline_model$estimates_best$fit, lty = 1, lwd  = line_lwd)
  points(x = data_farmacos$urb, y = data_farmacos$Sertraline, pch = 16, cex = cex.points)
}


pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/Antidepressivos_GAM.pdf", width = 6, height = 6)
par(mfrow = c(4,3), mar = c(3.5,3.5,0.5,0.5))

O.desmethylvenlafaxine_plot()
Desmethylcitalopram_plot()
Hydroxybupropion_plot()
Amitriptyline_plot()
Sertraline_plot()
Venlafaxine_plot()
Fluoxetine_plot()
Citalopram_plot()
Bupropion_plot()
Trimipramine_plot()
Norfluoxetina_plot()

dev.off()

