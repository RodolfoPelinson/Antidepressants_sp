Frequency of molecules across the gradient in urban cover
================
Rodolfo Pelinson
2024-11-14

The frequency data:

``` r
data.frame(urban_cover = Decil_unique, as.data.frame(frequencies))
```

    ##    urban_cover O.desmethylvenlafaxine Desmethylcitalopram Hydroxybupropion
    ## 1            0                   0.00                0.00             0.00
    ## 2           10                   0.75                0.25             0.75
    ## 3           20                   0.75                0.25             0.50
    ## 4           30                   1.00                0.00             0.75
    ## 5           40                   0.75                0.50             0.75
    ## 6           50                   1.00                0.00             0.50
    ## 7           60                   1.00                0.50             1.00
    ## 8           70                   0.75                0.25             0.50
    ## 9           80                   0.80                0.40             0.40
    ## 10          90                   1.00                0.75             0.75
    ## 11         100                   1.00                1.00             0.75
    ##    Amitriptyline Sertraline Venlafaxine Fluoxetine Citalopram Bupropion
    ## 1           0.00       0.00        0.00       0.00       0.00      0.00
    ## 2           0.00       0.00        0.25       0.00       0.00      0.25
    ## 3           0.25       0.25        0.50       0.00       0.25      0.00
    ## 4           0.25       0.25        0.00       0.25       0.00      0.00
    ## 5           0.75       0.50        0.75       0.75       0.00      0.00
    ## 6           0.00       0.25        0.25       0.00       0.00      0.25
    ## 7           0.50       0.50        0.50       0.25       0.50      0.25
    ## 8           0.50       0.50        0.50       0.50       0.25      0.00
    ## 9           0.60       0.80        0.60       0.60       0.20      0.00
    ## 10          0.75       1.00        0.75       0.75       0.50      0.25
    ## 11          0.50       0.75        0.50       0.50       0.25      0.25
    ##    Trimipramine
    ## 1          0.00
    ## 2          0.00
    ## 3          0.00
    ## 4          0.00
    ## 5          0.25
    ## 6          0.00
    ## 7          0.00
    ## 8          0.00
    ## 9          0.00
    ## 10         0.00
    ## 11         0.00

Fitting the loess functions:

``` r
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/functions/loess_func.R")

span <- 1
weights <- c(100,rep(1,10))
degree <- 2


Sertraline_model <- loess_func(frequencies$Sertraline, x =Decil_unique,  span = span, weights = weights, degree = degree) 
O.desmethylvenlafaxine_model <- loess_func(frequencies$O.desmethylvenlafaxine, x =Decil_unique, span = span, weights = weights, degree = degree)
Hydroxybupropion_model <- loess_func(frequencies$Hydroxybupropion, x =Decil_unique,  span = span, weights = weights, degree = degree)
Bupropion_model <- loess_func(frequencies$Bupropion, x =Decil_unique,  span = span, weights = weights, degree = degree)
Venlafaxine_model <- loess_func(frequencies$Venlafaxine, x =Decil_unique,  span = span, weights = weights, degree = degree)
Desmethylcitalopram_model <- loess_func(frequencies$Desmethylcitalopram, x =Decil_unique,  span = span, weights = weights, degree = degree)
Citalopram_model <- loess_func(frequencies$Citalopram, x =Decil_unique,  span = span, weights = weights, degree = degree)
Fluoxetine_model <- loess_func(frequencies$Fluoxetine, x =Decil_unique,  span = span, weights = weights, degree = degree)
Amitriptyline_model <- loess_func(frequencies$Amitriptyline, x =Decil_unique,  span = span, weights = weights, degree = degree)
Trimipramine_model <- loess_func(frequencies$Trimipramine, x =Decil_unique,  span = span, weights = weights, degree = degree)
```

Code to plot these curves:

``` r
############### Sertraline ###############

Sertraline_plot <- function(){
  plot(frequencies$Sertraline ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "sertraline", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Sertraline_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Sertraline, pch = 16, cex = cex.points)
  box()
}



############### O.desmethylvenlafaxine ###############

O.desmethylvenlafaxine_plot <- function(){
  plot(frequencies$O.desmethylvenlafaxine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "O-desmethylvenlafaxine", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), O.desmethylvenlafaxine_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$O.desmethylvenlafaxine, pch = 16, cex = cex.points)
  box()
}


############### Hydroxybupropion ###############

Hydroxybupropion_plot <- function(){
  plot(frequencies$Hydroxybupropion ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Hydroxybupropion", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Hydroxybupropion_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Hydroxybupropion, pch = 16, cex = cex.points)
  box()
}




############### Bupropion ###############

Bupropion_plot <- function(){
  plot(frequencies$Bupropion ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Bupropion", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Bupropion_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Bupropion, pch = 16, cex = cex.points)
  box()
}



############### Venlafaxine ###############

Venlafaxine_plot <- function(){
  plot(frequencies$Venlafaxine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Venlafaxine", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Venlafaxine_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Venlafaxine, pch = 16, cex = cex.points)
  box()
}




############### Desmethylcitalopram ###############

Desmethylcitalopram_plot <- function(){
  plot(frequencies$Desmethylcitalopram ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Desmethylcitalopram", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Desmethylcitalopram_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Desmethylcitalopram, pch = 16, cex = cex.points)
  box()
}




############### Citalopram ###############

Citalopram_plot <- function(){
  plot(frequencies$Citalopram ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Citalopram", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Citalopram_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Citalopram, pch = 16, cex = cex.points)
  box()
}



############### Norfluoxetina ###############

#Norfluoxetina_model <- loess_func(frequencies$Norfluoxetina, x =Decil_unique,  span = span)

#Norfluoxetina_plot <- function(){
#  plot(frequencies$Norfluoxetina ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
#  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
#  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
#  title(ylab = "Norfluoxetine", line = 1.25, cex.lab = ylab.cex)
#  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
#  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
#  lines(c(0:100), Norfluoxetina_model, lty = 1, lwd  = line_lwd)
#  points(x = Decil_unique, y = frequencies$Norfluoxetina, pch = 16, cex = cex.points)
#  box()
#}



############### Fluoxetine ###############

Fluoxetine_plot <- function(){
  plot(frequencies$Fluoxetine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Fluoxetine", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Fluoxetine_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Fluoxetine, pch = 16, cex = cex.points)
  box()
}



############### Amitriptyline ###############

Amitriptyline_plot <- function(){
  plot(frequencies$Amitriptyline ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Amitriptyline", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Amitriptyline_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Amitriptyline, pch = 16, cex = cex.points)
  box()
}



############### Trimipramine ###############

Trimipramine_plot <- function(){
  plot(frequencies$Trimipramine ~ Decil_unique , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of", line = 2, cex.lab = ylab.cex)
  title(ylab = "Trimipramine", line = 1.25, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(1, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = tck); axis(2, tick = FALSE, cex.axis = axis.cex, line = line_axis, gap.axis = -10)
  lines(c(0:100), Trimipramine_model, lty = 1, lwd  = line_lwd)
  points(x = Decil_unique, y = frequencies$Trimipramine, pch = 16, cex = cex.points)
  box()
}
```

Plotting the curves:

``` r
ylab.cex <- 0.9
xlab.cex <- 0.9
axis.cex <- 0.75
line_labs <- 1.25
line_axis <- -0.75
tck <- -0.03
cex.points <- 0.7
line_lwd <- 1.25
cex.letter <- 1.25
x_letter = 2
y_letter = 98


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
#Norfluoxetina_plot()
```

<img src="frequency_molecules_loess_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Plotting the curves together:

``` r
lines <- list(O.desmethylvenlafaxine_model,
              Desmethylcitalopram_model,
              Hydroxybupropion_model,
              Amitriptyline_model,
              Sertraline_model,
              Venlafaxine_model,
              Fluoxetine_model,
              Citalopram_model,
              Bupropion_model,
              Trimipramine_model)
              #Norfluoxetina_model)

library(colorspace)
library(shape)


plot_loess <- function(letter = "", x_letter = 0, y_letter = 100){
  
  length(lines)
  Cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  Cols <- c(rep(c("#F0E442", "#D55E00", "#56B4E9", "#009E73", "#CC79A7"),2))
  bg <- c(rep("black", 5),rep("white", 5))
  names <- names(frequencies)
  
  margins<- par()$mar
  par(mar = c(7,3.5,2,0.5))
  
  plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0,1), xlim = c(0,100))
  title(xlab = "Urban cover (%)", line = line_labs, cex.lab = xlab.cex)
  title(ylab = "Frequency of API", line = line_labs, cex.lab = ylab.cex)
  axis(1, tick = TRUE, line = 0, labels = FALSE,  tck = -0.02); axis(1, tick = FALSE, cex.axis = 0.8, line = -0.7, gap.axis = -10)
  axis(2, tick = TRUE, line = 0, labels = FALSE,  tck = -0.02); axis(2, tick = FALSE, cex.axis = 0.8, line = -0.7, gap.axis = -10)
  
  
  for(i in 1:length(lines)){
    lines(c(0:100), lines[[i]], lty = 1, lwd  = 3.5, col = bg[i])
    lines(c(0:100), lines[[i]], lty = 1, lwd  = 2, col = Cols[i])
  }
  box()
  
  #plot(NA , type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xlim = c(0,100), ylim = c(0,100))
 
  
  par(new=TRUE, mar = c(0,0,0,0)) 
  plot(NA, type = "n", xlim = c(0,100), xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = F, ylim = c(0,100))
  
  legend(x = 2.5, y = 15,
         legend = names, bty = "n", col = bg, cex = 0.9, ncol = 3, lty = 1, lwd = 3.5, text.width   = c(32.5, 17.5))
  
  legend(x = 2.5, y = 15,
         legend = names, bty = "n", col = Cols, cex = 0.9, ncol = 3, lty = 1, lwd = 2, text.width   = c(32.5, 17.5))
  
  text(x = x_letter, y = y_letter, labels = letter, font = 2, cex = 1.5, adj = 0)
  par(mar = margins)
  
}


par(mfrow= c(1,1))
plot_loess()
```

<img src="frequency_molecules_loess_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />
