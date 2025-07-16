#!/usr/bin/env Rscript
if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("./F_0.5/Intervene_upset.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('hippocampus'=3141,'dlpfc'=3716,'dlpfc&hippocampus'=4270,'caudate'=7835,'caudate&hippocampus'=863,'caudate&dlpfc'=709,'caudate&dlpfc&hippocampus'=2178)
upset(fromExpression(expressionInput), nsets=3, nintersects=30, show.numbers="yes", main.bar.color="#36383B", sets.bar.color="#36383B", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
