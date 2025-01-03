RDA
================
Rodolfo Pelinson
2024-11-14

``` r
source("C:/Users/rodol/OneDrive/repos/Antidepressants_sp/scripts/ajeitando_planilhas.R")

pred_mod_sel <- select(predictors, area, declividade, IC, discharge, imperm, casas_sem_saneamento, renda, habitantes)

pred_names <- c("Area", "Slope", "CI", "Discharge", "Impermeabilization", "HWNS", "Average Income", "Inhabitants")

colnames(pred_mod_sel) <- pred_names

colnames(data_farmacos) <-gsub("q_","", colnames(data_farmacos))
colnames(data_farmacos) <-gsub("_","-", colnames(data_farmacos))


pred_mod_sel_stand <- decostand(pred_mod_sel,  method = "stand")

rda <- rda(X = data_farmacos, Y =  pred_mod_sel_stand)

anova(rda)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(X = data_farmacos, Y = pred_mod_sel_stand)
    ##          Df Variance      F Pr(>F)  
    ## Model     8   5406.1 2.3209  0.011 *
    ## Residual 42  12228.8                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj (rda)
```

    ## $r.squared
    ## [1] 0.306558
    ## 
    ## $adj.r.squared
    ## [1] 0.1744738

``` r
#plot(rda, scaling = 0)
#plot(rda, scaling = 1)
#plot(rda, scaling = 2)
#plot(rda, scaling = 3)
```

``` r
library(shape)
library(vegan)
library(dplyr)
library(scales)
```

``` r
pal <- col_numeric(
  palette = c("white", "black"),
  domain = urb,
  na.color = "red",
  alpha = FALSE,
  reverse = FALSE)
colors_urb <- pal(urb)

col_classes <- data_farmacos_class
col_classes[col_classes == "SSRI"] <- "#0072B2"
col_classes[col_classes == "SNRI"] <- "#D55E00"
col_classes[col_classes == "Aminoketone"] <- "#E69F00"
col_classes[col_classes == "Tricyclic"] <- "#009E73"


col_classes <- c(col_classes, "black")


#####################################################

scores <- scores(rda, scaling  = 3)

loadings <- scores$biplot * attr(scores, "const")
rdas_bacias <- scores$sites
rdas_farmacos <- scores$species

eigenvalues <- rda$CCA$eig
importance <- eigenvalues/sum(eigenvalues)


RDA1_label <- paste("RDA1 (",round(importance[1]*100,2),"%)",sep = "")
RDA2_label <- paste("RDA2 (",round(importance[2]*100,2),"%)",sep = "")


xmin <- min(c(rdas_bacias[,1], loadings[,1], rdas_farmacos[,1]))*1.3
xmax <- max(c(rdas_bacias[,1], loadings[,1], rdas_farmacos[,1]))*1.3
ymin <- min(c(rdas_bacias[,2], loadings[,2], rdas_farmacos[,2]))*1.2
ymax <- max(c(rdas_bacias[,2], loadings[,2], rdas_farmacos[,2]))*1.2



#pdf("Resultados/Resultados Antidepressivos/Figuras/PDF/rda.pdf", height = 3.5* 1.1, width = (4 * 1.2), pointsize = 8)

par(bg = "white", mar = c(6,4,.1,5), bty = "o", cex = 1)
plot(rdas_bacias[,1], rdas_bacias[,2], xlim = c(xmin,xmax), ylim = c(ymin, ymax),
     type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")

abline(h = 0, v = 0, lty = 2)



#####################
shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}
#####################

shadowtext(x = rdas_farmacos[,1], y = rdas_farmacos[,2], labels = rownames(rdas_farmacos), col=col_classes, bg="white", cex = 0.9, r = .1)

points(rdas_bacias, pch = 21, cex = 1.2, bg = colors_urb)


Arrows(x0 <- rep(0, nrow(loadings)),
       y0 <- rep(0, nrow(loadings)),
       x1 <- loadings[,1],
       y1 <- loadings[,2], arr.type = "triangle", arr.length = 0.2, col = "grey30", lwd = 1.5)



#text(x = loadings[,1]+ 0.5, y = loadings[,2]+ 0.5, labels = rownames(loadings), col = col_classes, adj = c(0,1), cex = 0.9)
shadowtext(x = loadings[,1], y = loadings[,2], labels = rownames(loadings), col="grey30", bg="white", cex = 0.9, r = .1)



axis(1, cex.axis = 1)
axis(2, cex.axis = 1, las = 2)
title(xlab = RDA1_label, cex.lab = 1.1, line = 2.5)
title(ylab = RDA2_label, cex.lab = 1.1, line = 2.5)
#title(main = "Pesticides", line = 0.5, adj = 0, cex.main = 1.75)

par(new=TRUE, mar = c(0,0,0,0)) 
plot(NA, type = "n", xlim = c(0,100), xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = F, ylim = c(0,100))

legend(x = 0, y = 5,
       legend = unique(data_farmacos_class), bty = "n", cex = 1, ncol = 4, lty = 0, lwd = 3,text.width  = c(10, 20, 10), text.col = unique(col_classes), text.font = 2)

#library(grDevices)
legend_image <- as.raster(matrix(pal(0:100), ncol=1))
#plot(c(80,90),c(10,90),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 92, 22, 97,92)
axis(4, at=c(22, (22+92)/2, 92), pos=97, labels= F, tck=-.01, las = 1, cex.axis = 0.8) 
axis(4, at=c(22, (22+92)/2, 92), pos=96, labels= c("100","50","0"), tick = FALSE, las = 1, cex.axis = 0.8) 

text(x = c(94,94), y = c(100,95), labels = c("Urban", "cover"), adj = 0.5, cex.axis = 0.8)
```

<img src="RDA_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />
