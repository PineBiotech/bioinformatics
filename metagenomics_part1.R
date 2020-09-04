#Test your R code here

library(openxlsx)

#load data
bacteria_order <- read.xlsx("metagenomics1Annotation.xlsx", sheet = 3)

#create a table of results
numtable <- table(bacteria_order$ERR1250035)

prop.table(numtable)
labels <- unique(as.vector(bacteria_order$byphylum))

print(labels)
colnum <- length(labels)

par(mar=c(2,2,2,2))
pie(numtable, col=heat.colors(colnum), main="Summary Pie Chart", c(labels))

