



install.packages("renv")
renv::init()

#Snapshot?
renv::snapshot()
renv::restore()



#Packages
install.packages("usethis")


install.packages("devtools")
library(devtools)
#install.packages("Rcpp")
#devtools::install_github("JenniNiku/gllvm")
install.packages("vegan")
install.packages("glmmTMB", type = "source")
install.packages("emmeans")
install.packages("car")
install.packages("DHARMa")
install.packages("bbmle")
install.packages("Matrix", type = "source")
install.packages("TMB", type = "source")
install.packages("sdmTMB", type = "source")
install.packages("randomForest")
install.packages("mgcv")


#install.packages("spaa")
#install.packages("doParallel")
#install.packages("pbapply")


usethis::use_git()
usethis::use_github()
usethis::use_readme_rmd()
