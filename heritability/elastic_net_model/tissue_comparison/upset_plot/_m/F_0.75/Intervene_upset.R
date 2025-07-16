#!/usr/bin/env Rscript
if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("./F_0.75/Intervene_upset.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('hippocampus'=4134,'dlpfc'=4585,'dlpfc&hippocampus'=4258,'caudate'=8367,'caudate&hippocampus'=791,'caudate&dlpfc'=657,'caudate&dlpfc&hippocampus'=1770)
upset(fromExpression(expressionInput), nsets=3, nintersects=30, show.numbers="yes", main.bar.color="#36383B", sets.bar.color="#36383B", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
