#!/usr/bin/env Rscript
if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("./f_0.25/Intervene_upset.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('hippocampus'=2404,'dlpfc'=3177,'dlpfc&hippocampus'=4063,'caudate'=7402,'caudate&hippocampus'=883,'caudate&dlpfc'=760,'caudate&dlpfc&hippocampus'=2540)
upset(fromExpression(expressionInput), nsets=3, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
