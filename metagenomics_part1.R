#Test your R code here

library(openxlsx)
library(dplyr)

#load data
annotatedtb <- read.xlsx("metagenomics1Annotation.xlsx", sheet = 4)
df <- data.frame(annotatedtb)

#create a table of summary values by phylum
numtable1 <- rowsum(df$ERR1250035, df$phylum, reorder = FALSE)
numtable1 <- na.omit(numtable1)
totalsum <- sum(numtable1)
numtable2 <- lapply(numtable1[,1], function(x) { x/totalsum })
numtable2 <- lapply(numtable2, function(x) { round(x*100, digits=2)})


#add labels and count length
labels <- row.names(numtable1)
colnum <- length(labels)
colorfill <- factor(c(labels))
numtable2 <- lapply(numtable2, function(x) { paste(colorfill[x], x, " %")})

#make a pie chart
par(mar=c(4,4,4,4))
par(xpd=TRUE)
pie(table(annotatedtb[,1]), col = rainbow(colnum), radius = 01, main="ERR1250035 by Phylum", c(numtable2), cex = 0.6)
legend(-2,1,labels,c(rainbow(colnum)), cex = 0.5)
par(mar=c(6,4,4,4))
barplot(table(annotatedtb[,1]), col = rainbow(colnum), main="ERR1250035 by Phylum", xlab = "Bacteria Type", las = 2, cex.names = 0.5)
