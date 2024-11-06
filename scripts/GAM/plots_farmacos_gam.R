data_farmacos <- read.csv("antidepressivos_r.csv")
data_general <- read.csv("planilha_master_r.csv")

data_farmacos_reference <- data_farmacos[1,]
data_farmacos <- data_farmacos[-1,]

data_farmacos$urb <- data_general$urb[match(data_farmacos$bacia,data_general$microbacia)]

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)

data_farmacos$urb[data_farmacos$bacia == "EBB_08" | 
                    data_farmacos$bacia == "EBB_09" |
                    data_farmacos$bacia == "EBB_11" |
                    data_farmacos$bacia == "EBB_12"] <- 0 

data_farmacos$urb[is.na(data_farmacos$urb)] <- data_farmacos$Decil[is.na(data_farmacos$urb)] + 5

data.frame(data_farmacos$bacia, data_farmacos$urb, data_farmacos$Decil)



############### O-desmethylvenlafaxine
png("Plots/farmacos/O-desmethylvenlafaxine.png", width = 10, height = 10, units = "cm", pointsize = 10, res = 600)

par(mar = c(3.5,3.5,1,0.1))
plot(data_farmacos$O.desmethylvenlafaxine ~ data_farmacos$urb, xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = 1)
title(xlab = "Urbanização (%)", line = 2, cex.lab = 1.5)
title(ylab = "O-desmethylvenlafaxine (10 ng/L)", line = 2, cex.lab = 1.5)
axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

dev.off()


library(mgcv)
library(bbmle)


gam_O.desmethylvenlafaxine_0 <-   gam(O.desmethylvenlafaxine ~ 1, family = "tw", method = "REML", select = TRUE, data = data_farmacos)
gam_O.desmethylvenlafaxine_1 <-   gam(O.desmethylvenlafaxine ~ s(urb, k = 3, fx = FALSE, bs = "tp") + 0, family = "tw", method = "REML", select = TRUE, data = data_farmacos)
gam_O.desmethylvenlafaxine_2 <-   gam(O.desmethylvenlafaxine ~ s(urb, k = 6, fx = FALSE, bs = "tp") + 0, family = "tw", method = "REML", select = TRUE, data = data_farmacos)
gam_O.desmethylvenlafaxine_3 <-   gam(O.desmethylvenlafaxine ~ s(urb, k = 12, fx = FALSE, bs = "tp") + 0, family = "tw", method = "REML", select = TRUE, data = data_farmacos)
gam_O.desmethylvenlafaxine_4 <-   gam(O.desmethylvenlafaxine ~ s(urb, k = 24, fx = FALSE, bs = "tp") + 0, family = "tw", method = "REML", select = TRUE, data = data_farmacos)

summary(gam_O.desmethylvenlafaxine_1)
summary(gam_O.desmethylvenlafaxine_2)
summary(gam_O.desmethylvenlafaxine_3)
summary(gam_O.desmethylvenlafaxine_4)


AICctab(gam_O.desmethylvenlafaxine_0, gam_O.desmethylvenlafaxine_1, gam_O.desmethylvenlafaxine_2, gam_O.desmethylvenlafaxine_3, gam_O.desmethylvenlafaxine_4)


new_urbanization <- data.frame(urb = c(0:99))
pred_abundance_gam <- predict(gam_O.desmethylvenlafaxine_3, se = TRUE, type = "link", newdata = new_urbanization)

fit_gam <- pred_abundance_gam$fit
upper_gam <- exp(fit_gam + (1.96*pred_abundance_gam$se.fit))
lower_gam <- exp(fit_gam - (1.96*pred_abundance_gam$se.fit))
fit_gam <- exp(fit_gam)

plot(data_farmacos$O.desmethylvenlafaxine ~ data_farmacos$urb , type = "p", ylim = c(0,600), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
title(xlab = "Urban coverage (%)", line = 2, cex.lab = 1.5)
title(ylab = "O-desmethylvenlafaxine (ng/L)", line = 2, cex.lab = 1.5)
abline(h = 10, lty = 2, col = "brown")
axis(1, tick = TRUE, line = 0, labels = FALSE); axis(1, tick = FALSE, cex.axis = 1, line = -0.5)
axis(2, tick = TRUE, line = 0, labels = FALSE); axis(2, tick = FALSE, cex.axis = 1, line = -0.5)

lines(new_urbanization$urb, fit_gam, lty = 1)
lines(new_urbanization$urb, upper_gam, lty = 2)
lines(new_urbanization$urb, lower_gam, lty = 2)
