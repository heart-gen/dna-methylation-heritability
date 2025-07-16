#!/usr/bin/env Rscript
if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("./f_0.75/Intervene_upset.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('hippocampus'=3395,'dlpfc'=4413,'dlpfc&hippocampus'=3710,'caudate'=9144,'caudate&hippocampus'=640,'caudate&dlpfc'=581,'caudate&dlpfc&hippocampus'=1220)
upset(fromExpression(expressionInput), nsets=3, nintersects=30, show.numbers="yes", main.bar.color="#36383B", sets.bar.color="#36383B", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
