My_Imagine <- function (comm, col = c(0,1,"grey50"), order = TRUE, scores = 1, fill = TRUE,
                        xlab = "", ylab = "", yline = 0, xline = 0, sitenames = rownames(comm), cex.envlab = 1,
                        speciesnames = colnames(comm),   xlab_line = 3, cex.xlab = 2, ylab_line = 3, cex.ylab = 2, main = "", main_line = 3, cex.main = 2,
                        Env1 = NULL, Env2 = NULL, Env.col_1 = NULL, Env.label_1 = NULL, Env.col_2 = NULL, Env.label_2 = NULL,
                        cex.site = 1, cex.species = 1, top_margin = 2, left_margin = 2, bottom_margin = 3, right_margin = 1,EigenVal = F, box.lwd = 1,
                        Env.col_3 = NULL, Env.label_3 = NULL, Env3 = NULL,
                        Env.col_4 = NULL, Env.label_4 = NULL, Env4 = NULL,
                        Env.col_5 = NULL, Env.label_5 = NULL, Env5 = NULL,
                        gap.axis = 0, Empty = NULL, Empty_col = "grey75", sitesfont = 1, speciesfont = 1)
  
  
{
  
  require(scales)
  
  if(isFALSE(is.null(Env1))  & is.null(Env2) & is.null(Env3)& is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(0.95,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & is.null(Env3)& is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(0.9,0.05,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths = c(0.85,0.05,0.05,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & isFALSE(is.null(Env4)) & is.null(Env5)){
    layout(matrix(c(1,2,3,4,5), 1, 5, byrow = TRUE), widths = c(0.868,0.033,0.033,0.033,0.033))
  }
  if(isFALSE(is.null(Env1))  & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & isFALSE(is.null(Env4)) & isFALSE(is.null(Env5))){
    layout(matrix(c(1,2,3,4,5,6), 1, 6, byrow = TRUE), widths = c(0.835,0.033,0.033,0.033,0.033,0.033))
  }
  
  require(metacom)
  if (order == TRUE) {
    new_scores <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = T)
    if(is.null(sitenames)==FALSE){
      sitenames <- sitenames[order(new_scores$sitescores)]
    }
    if(is.null(speciesnames)==FALSE){
      speciesnames <- speciesnames[order(new_scores$speciesscores)]
    }
    if(is.null(Env1) == FALSE){
      Env1 <- Env1[order(new_scores$sitescores)]
    }
    if(is.null(Env2) == FALSE){
      Env2 <- Env2[order(new_scores$sitescores)]
    }
    
    if(is.null(Env3) == FALSE){
      Env3 <- Env3[order(new_scores$sitescores)]
    }
    
    if(is.null(Env4) == FALSE){
      Env4 <- Env4[order(new_scores$sitescores)]
    }
    
    if(is.null(Env5) == FALSE){
      Env5 <- Env5[order(new_scores$sitescores)]
    }
    
    
    
    EigenValues <- decorana(comm, ira = 1)$evals
    Relative_EigenValues <- EigenValues/sum(EigenValues)
    comm <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = FALSE)
    
  }
  
  
  #Filling embedded absences
  if (fill == TRUE) {
    temp1 = comm
    for (i in 1:dim(comm)[2]) {
      temp2 = comm[, i]
      if (sum(temp2) < 2) {
        comm[, i] = temp2
      }
      else {
        first <- min(which(temp2 > 0))
        last <- max(which(temp2 > 0))
        temp1[first:last, i] <- 1
      }
    }
    for (i in 1:dim(comm)[1]){
      for (j in 1:dim(comm)[2]){
        if(temp1[i,j] == 1 & comm[i,j] == 0){
          comm[i,j] <- 2
        }
      }
    }
  }
  
  #Reversing the order of rows?
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  
  
  #Column to be empty communities?
  if(is.null(Empty) == FALSE){
    new_empty <- comm[,Empty]
    new_empty[which(new_empty == 1)] <- 3
    comm[,Empty] <- new_empty
  }
  
  
  
  #ploting
  par(mar = c(bottom_margin,left_margin, top_margin, right_margin), cex = 1)
  if(is.null(Empty) == FALSE){
    image(1:dim(comm)[2], 1:dim(comm)[1], t(comm),  col = c(col, Empty_col),
          xlab = "", ylab = "", axes = FALSE)
  }else{
    image(1:dim(comm)[2], 1:dim(comm)[1], t(comm),  col = c(col),
          xlab = "", ylab = "", axes = FALSE)
  }
  
  title(xlab = xlab, line = xlab_line, cex.lab = cex.xlab, adj = 1)
  title(ylab = ylab, line = ylab_line, cex.lab = cex.ylab)
  title(main = main, line = main_line, cex.main = cex.main, adj = 0, font.main = 1)
  
  
  box(lwd = box.lwd)
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1,
         cex.axis = cex.site, tick = FALSE, line = yline, gap.axis = gap.axis, font.axis = sitesfont)
  }
  if (length(speciesnames) > 1) {
    if(is.null(Empty) == FALSE){
      axis(3, at = Empty, labels = speciesnames[Empty], las = 2,
           cex.axis = cex.species, tick = FALSE, line = xline, gap.axis = gap.axis, font.axis = 1)
      axis(3, at = (1:dim(comm)[2])[-Empty], labels = speciesnames[-Empty], las = 2,
           cex.axis = cex.species, tick = FALSE, line = xline, gap.axis = gap.axis, font.axis = speciesfont)
    }else{
      axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2,
           cex.axis = cex.species, tick = FALSE, line = xline, gap.axis = gap.axis, font.axis = speciesfont)
    }

  }
  
  if(is.null(Env1)==F){
    
    if(is.numeric(Env1)){
      pal <- col_numeric(
        palette = c(Env.col_1[1], Env.col_1[2]),
        domain = Env1,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_1 <-pal(Env1)
    }
    
    if(is.factor(Env1)){
      Env1<-as.numeric(Env1)
    }
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env1),x = 1, z = t(as.matrix(1:length(Env1))),
          col = Env.col_1,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_1, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
    
  }
  
  if(is.null(Env2)==F){
    
    if(is.numeric(Env2)){
      pal <- col_numeric(
        palette = c(Env.col_2[1], Env.col_2[2]),
        domain = Env2,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_2 <-pal(Env2)
    }
    
    if(is.factor(Env2)){
      Env2<-as.numeric(Env2)
    }
    
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env2),x = 1, z = t(as.matrix(1:length(Env2))),
          col = Env.col_2,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_2, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }
  
  if(is.null(Env3)==F){
    
    if(is.numeric(Env3)){
      pal <- col_numeric(
        palette = c(Env.col_3[1], Env.col_3[2]),
        domain = Env3,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_3 <-pal(Env3)
    }
    
    if(is.factor(Env3)){
      Env3<-as.numeric(Env3)
    }
    
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env3),x = 1, z = t(as.matrix(1:length(Env3))),
          col = Env.col_3,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_3, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }
  
  
  if(is.null(Env4)==F){
    
    if(is.numeric(Env4)){
      pal <- col_numeric(
        palette = c(Env.col_4[1], Env.col_4[2]),
        domain = Env4,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE
      )
      Env.col_4 <- pal(Env4)
    }
    
    if(is.factor(Env4)){
      Env4<-as.numeric(Env4)
    }
    
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env4),x = 1, z = t(as.matrix(1:length(Env4))),
          col = Env.col_4,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_4, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }
  
  
  if(is.null(Env5)==F){
    
    if(is.numeric(Env5)){
      pal <- col_numeric(
        palette = c(Env.col_5[1], Env.col_5[2]),
        domain = Env5,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_5 <-pal(Env5)
    }
    
    if(is.factor(Env5)){
      Env5<-as.numeric(Env5)
    }
    
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env5),x = 1, z = t(as.matrix(1:length(Env5))),
          col = Env.col_5,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_5, las = 2,
         cex.axis = cex.envlab, tick = FALSE,  line = xline)
    box(lwd = box.lwd)
  }
  
  
  
  if(isTRUE(EigenVal)){
    y <- -((dim(comm)[2]^2)/10000)
    text(x=(dim(comm)[2]),y=y, labels = paste(round(Relative_EigenValues[scores]*100,2),"%", sep = " "), xpd=NA, adj=c(1,1))
  }
  
  if(isTRUE(EigenVal)){
    y <- -((dim(comm)[2]^2)/10000)
    text(x=(dim(comm)[2]),y=y, labels = paste(round(Relative_EigenValues[scores]*100,2),"%", sep = " "), xpd=NA, adj=c(1,1))
  }
  
  layout(matrix(c(1), 1, 1))
  par(mar = c(4,4,4,4))
}


OrderGradient <- function (comm, gradient, scores = 1){
  new_scores <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = T)
  gradient <- gradient[order(new_scores$sitescores)]
  return(gradient)
}

OrderSpecies <- function (comm, speciesnames, scores = 1){
  new_scores <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = T)
  speciesnames <- speciesnames[order(new_scores$speciesscores)]
  return(speciesnames)
}

OrderSites <- function (comm, sitesnames, scores = 1){
  new_scores <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = T)
  sitesnames <- speciesnames[order(new_scores$sitescores)]
  return(speciesnames)
}


