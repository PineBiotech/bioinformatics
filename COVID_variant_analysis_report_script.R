library(ggplot2)
library(reshape)
library(dplyr)
library(gt)
library(ggpubr)
library(gridExtra)

filenames <- list.files(pattern="*.txt")

#Strain Name
UKVOC <- "UK: 20I/501Y.V1"
SAVOC <- "SA: 20H/501Y.V2"

Report_Type <- c("UK","SA")

VOC <- c(UKVOC,SAVOC)
#Provide variant positions
UKpos <- factor(c(8782,14408,14676,15279,23063,23403,28881,28882))
SApos <- factor(c(1059,8782,14408,23403,25563,23063,23012))

make_report <- function(sample_name1,VOC){

  report_type <- ifelse(VOC==UKVOC,"UK","SA")
  
  if(report_type=="UK"){
    my_positions <- UKpos
  } else {
    my_positions <- SApos
  }
  
  report_type1 <- as.character(paste0("_",report_type,".pdf"))
  report_name <- toString(sub(".txt", report_type1, sample_name1))
  
  print(report_type)
  print(my_positions)
  print(UKpos)
  print(report_name)

#Read data having variant coverage information
mydata1 <- read.table(sample_name1, sep="\t", header=TRUE)

#remove insertion type column
mydata1 <- mydata1[,-10]

#Compute coverage
Coverage <- rowSums(mydata1[,4:9])

#Add coverage column to the data
mydata1 <- cbind(mydata1,Coverage)

#Extract selected position columns
mydata2 <- mydata1[mydata1$Pos %in% my_positions,-1]

#Remove extra column and make final data with vital info
cleanCoverage <- mydata2[,-9]
cleanCoverage

#Melt data
mydataMelt <- melt(cleanCoverage, id=c("Pos","Nucl"))
mydataMelt

#Remove coverage rows
mydataMelt1 <- subset(mydataMelt, !(variable == "Coverage"))

# Findvariant information
TopVariants <- mydataMelt1 %>% 
  group_by(Pos) %>% 
  mutate(max_score = max(value)) %>% 
  ungroup() %>% 
  filter(value == max_score)

#Print top variants
TopVariants

#Extract coverage
mydataMelt_cov <- mydataMelt %>% filter(variable == "Coverage")

#Complete data with coverage
cov <- mydataMelt_cov %>% 
  group_by(Pos) %>% 
  mutate(max_score = max(value)) %>% 
  ungroup() %>% 
  filter(value == max_score)

#Print mat
cov

#Create same order based on pos for variant and coverage 
TopVariants1 <- TopVariants[order(TopVariants$Pos),]
cov1 <- cov[order(cov$Pos),]

#Make final table with key information of variant and coverage
TopVariants_new <- cbind(TopVariants1[1:4], cov1$max_score)
TopVariants_new 

#Add Column names
colnames(TopVariants_new ) <- c("Pos", "Nucl", "Variant", "Var_reads", "Coverage")

#Compute %age for coverage
Fraction <- round((TopVariants_new$Var_reads/ TopVariants_new$Coverage) * 100, 1)
Fraction <- paste(Fraction,"%")

#Combine variant data with %age column
TopVariants_new1 <-cbind(TopVariants_new, Fraction)

#Compute whether it is a variant
TopVariants_new1$VarYN <- ifelse(as.character(TopVariants_new1$Variant) == as.character(TopVariants_new1$Nucl),"N","Y")

#make a table of top variants for report
TopVariants_new1 %>% gt()

#make a plotting theme
mytheme <- theme(
  text=element_text(family="Helvetica", size=10),
  axis.text.x = element_text(angle = 90,size = 5),
  axis.text.y = element_text(size = 5),
  panel.background = element_rect(fill = NA),
  panel.border = element_rect(fill = NA, color = "grey75"),
  axis.ticks = element_line(color = "grey85"),
  panel.grid.major = element_line(color = "grey95", size = 0.2),
  panel.grid.minor = element_line(color = "grey95", size = 0.2),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)

#summary table
variant_summary <- ggtexttable(TopVariants_new1, theme=ttheme("light"))
variant_summary <- variant_summary %>%
  tab_add_title(text = paste("Variant Summary for ",VOC), face = "plain", size = 10)

#plot coverage
p1 <-ggplot(data=mydata1, aes(x=Pos, y=Coverage, aspect_ratio = 2.5)) +
  geom_area(color="lightgray",fill="lightgray",alpha=0.5) +
  geom_point(aes(color = Coverage)) +
  scale_color_gradientn(colours=c('red','orange',"cyan","blue",'blueviolet')) +
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 50)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ggtitle(paste("Coverage for ",sample_name1)) +
  mytheme

#plot coverage per variant
mydata3 <- mydata2[,-10]
mydata4 <- cbind(mydata3,Var = TopVariants1$variable)
mydata4 <- melt(mydata4, id=c("Pos","Nucl", "Var"))

p2 <-ggplot(mydata4, aes(x=as.factor(Pos)), y=value) +
  geom_col(aes(y=value, fill=variable, group=Nucl), width = 0.7) +
  ggtitle(paste("Variant details for ",VOC)) +
  geom_text(aes(y = -100, label = Nucl), color = "black") +
  geom_text(aes(y = max(value)+200, label = Var), color = "red") +
  mytheme

grid <- grid.arrange(p1,arrangeGrob(p2,variant_summary,ncol=2), nrow=2)
ggsave(grid, file=report_name, height = 8 , width = 10)

}

for (i in 1:length(filenames)){
  make_report(sample_name1=filenames[i],VOC=UKVOC)
  make_report(sample_name1=filenames[i],VOC=SAVOC)
}


