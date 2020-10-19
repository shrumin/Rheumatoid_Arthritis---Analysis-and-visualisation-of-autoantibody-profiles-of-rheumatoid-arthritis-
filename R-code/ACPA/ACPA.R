setwd("~/Desktop/RA/ACPA")
table1=read.csv("ACPA_Data.csv",header = TRUE,row.names = 1)
RA_Patients=read.csv("final_experiment_info.csv",row.names = 1,header=TRUE)

#citrullinated data only
cit_data=read.csv("Cit_Data.csv",row.names = 1,header = TRUE)

#analysis using cutoff.csv
t1=read.csv("Cutoff.csv")
t1=t(t1)
cit_dataCopy=cit_data
for(i in 1:15){cit_dataCopy[,i]=cit_dataCopy[,i]-t1[i]}
count_function=function(x){length(which(x>0))}
c=apply(cit_dataCopy,1,count_function)
cit_dataCopy2=cit_data
cit_dataCopy2=merge(cit_dataCopy2,c,by=0)
row.names(cit_dataCopy2)=cit_dataCopy2$Row.names
cit_dataCopy2=cit_dataCopy2[,c(2:17)]
names(cit_dataCopy2)[names(cit_dataCopy2) == "y"] <- "Number_of_Positives"

# merging with RA_Patients 
merge1=merge(cit_dataCopy2,RA_Patients,by=0)
row.names(merge1)=merge1$Row.names
merge1=merge1[,c(2:18)]
write.csv(merge1,"Number_of_Positive-1.csv")
test=merge1[,c(16,17)]

#visualise in R
library("ggpubr")
test=merge1[,c(16,17)]
ggboxplot(test, x = "group", y = "Number_of_Positives", 
            color = "group",
            order = c("Negative","CCP", "CCP_Acet","CCP_Carb_Acet"),
             ylab = "Number_of_Positives", xlab = "group")
ggline(test, x = "group", y = "Number_of_Positives", 
           add=c("mean_se","jitter"),
           order = c("Negative","CCP", "CCP_Acet","CCP_Carb_Acet"),
           ylab = "Number_of_Positives", xlab = "group")

 
#Kruskal walis test
kruskal.test(Number_of_Positives ~ group, data = test)

# multiple pairwise comparison between groups
pairwise.wilcox.test(test$Number_of_Positives, test$group, p.adjust.method = "BH")

# dunn test
library(FSA)
PT = dunnTest(Number_of_Positives ~ group,data=test,method="bh")
PT
PT = PT$res
library(rcompanion)
cldList(comparison = PT$Comparison,p.value= PT$P.adj,threshold  = 0.05)

#########################################################################
#joint_narrowing with ACPA merging and analysis
narrowing_exp_info=read.csv("Joint_Narrowing_exp_info.csv",row.names = 1)
merge2=merge(cit_dataCopy2,narrowing_exp_info,by=0)
row.names(merge2)=merge2$Row.names
merge2=merge2[,c(2:18)]
write.csv(merge2,"Number_of_Positive-JointNarrowing.csv")
test2=merge2[,c(16,17)]

#visualising
library("ggpubr")
ggboxplot(test2, x = "group", y = "Number_of_Positives", 
          color = "group",
          order = c("Joint_Narrowing", "No_Joint_Narrowing"),
          ylab = "Number_of_Positives", xlab = "group")

ggline(test2, x = "group", y = "Number_of_Positives", 
       add=c("mean_se","jitter"),
       order = c("Joint_Narrowing","No_Joint_Narrowing"),
       ylab = "Number_of_Positives", xlab = "group")

#Kruskal walis test
kruskal.test(Number_of_Positives ~ group, data = test2)

#########################################################################
#progression with ACPA merging and analysis
progression_data=read.csv("Progression_exp_info.csv",row.names = 1)
merge3=merge(cit_dataCopy2,progression_data,by=0)
row.names(merge3)=merge3$Row.names
merge3=merge3[,c(2:18)]
write.csv(merge3,"Number_of_Positive-Progression.csv")
test3=merge3[,c(16,17)]

#Kruskal walis test
kruskal.test(Number_of_Positives ~ group, data = test3)

#########################################################################
#erosion with ACPA merging and analysis
erosion_data=read.csv("Erosion_exp_info.csv",row.names = 1)
merge4=merge(cit_dataCopy2,erosion_data,by=0)
row.names(merge4)=merge4$Row.names
merge4=merge4[,c(2:18)]
write.csv(merge4,"Number_of_Positive-Erosion.csv")
test4=merge4[,c(16,17)]

#visualise the results
library("ggpubr")
ggboxplot(test4, x = "group", y = "Number_of_Positives", 
          color = "group",
          order = c("Erosion", "No_Erosion"),
          ylab = "Number_of_Positives", xlab = "group")

ggline(test4, x = "group", y = "Number_of_Positives", 
       add=c("mean_se","jitter"),
       order = c("Erosion","No_Erosion"),
       ylab = "Number_of_Positives", xlab = "group")


#Kruskal walis test
kruskal.test(Number_of_Positives ~ group, data = test4)

########################################################################
#Factor analysis

#number of factors determined by eigen vectors
library(nFactors)
ev <- eigen(cor(cit_data))
ap <- parallel(subject=nrow(cit_data),var=ncol(cit_data), rep=100, cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

#principal component analysis
cit_data.pca=princomp(cit_data)
summary(cit_data.pca)
plot(cit_data.pca)
library(factoextra)
fviz_pca_biplot(cit_data.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
#exploratory factor analysis
fit_4 <- factanal(cit_data, factors=4, rotation="varimax")
print(fit_4, digits=2, cutoff=.3, sort=TRUE)
fit_4_reg <- factanal(cit_data, 4,scores = "regression", rotation="varimax")
head(fit_4_reg$scores)

#interpretation of factors
plot(fit_4$loadings[,1], 
     fit_4$loadings[,2],
     xlab = "Factor 1", 
     ylab = "Factor 2", 
     ylim = c(-1,1),
     xlim = c(-1,1),
     main = "Varimax rotation")
text(fit_4$loadings[,1]-0.08, 
     fit_4$loadings[,2]+0.08,
     colnames(cit_data),
     col="blue")
#########################################################################
##Vienn diagram

# number of peptide positive in each group
experimental_data=cit_dataCopy
experimental_data$Filaggrin_Positive=as.factor(experimental_data$CCP_1.cit>0)
experimental_data$Vimentin_positive=as.factor(experimental_data$Vim60_75.cit>0 | experimental_data$Vim2_17.cit>0)
experimental_data$Fibrinogen_positive=as.factor(experimental_data$Fib36_52.cit>0 | experimental_data$Fib573.cit>0 | 
experimental_data$Fib591.cit>0|experimental_data$Fib.Alpha36_50.cit>0|experimental_data$Fib.alpha621_635.cit>0|experimental_data$Fib.beta60_74.cit>0)
experimental_data$alpha_enolase_Positive=as.factor(experimental_data$CEP_1>0)
experimental_data$hnRNP_A3_Positive=as.factor(experimental_data$Pept.1>0|
experimental_data$Pept_Z1>0|experimental_data$Pept_Z2>0|experimental_data$Pept.5>0|experimental_data$Bla26>0)

#merging with RA patient profiles
experimental_data2=experimental_data[,c(16:20)]
experimental_merge=merge(experimental_data2,RA_Patients,by=0)
row.names(experimental_merge)=experimental_merge$Row.names
experimental_merge=experimental_merge[,c(2:7)]
write.csv(experimental_merge,"Positives_group_merged.csv")

#identifying  number of members in each combination groups

#Filaggrin
nrow(subset(experimental_merge,Filaggrin_Positive==TRUE&group=='CCP'))
nrow(subset(experimental_merge,Filaggrin_Positive==TRUE & group=='CCP_Acet'))
nrow(subset(experimental_merge,Filaggrin_Positive==TRUE & group=='CCP_Carb_Acet'))
nrow(subset(experimental_merge,Filaggrin_Positive==TRUE & group=='Negative'))

#Vimentin
nrow(subset(experimental_merge,Vimentin_positive==TRUE & group=='CCP'))
nrow(subset(experimental_merge,Vimentin_positive==TRUE & group=='CCP_Acet'))
nrow(subset(experimental_merge,Vimentin_positive==TRUE & group=='CCP_Carb_Acet'))
nrow(subset(experimental_merge,Vimentin_positive==TRUE & group=='Negative'))

#Fibrinogen
nrow(subset(experimental_merge,Fibrinogen_positive==TRUE & group=='CCP'))
nrow(subset(experimental_merge,Fibrinogen_positive==TRUE & group=='CCP_Acet'))
nrow(subset(experimental_merge,Fibrinogen_positive==TRUE & group=='CCP_Carb_Acet'))
nrow(subset(experimental_merge,Fibrinogen_positive==TRUE & group=='Negative'))

#alpha enolase
nrow(subset(experimental_merge,alpha_enolase_Positive==TRUE & group=='CCP'))
nrow(subset(experimental_merge,alpha_enolase_Positive==TRUE & group=='CCP_Acet'))
nrow(subset(experimental_merge,alpha_enolase_Positive==TRUE & group=='CCP_Carb_Acet'))
nrow(subset(experimental_merge,alpha_enolase_Positive==TRUE & group=='Negative'))

#hnRNP_A3
nrow(subset(experimental_merge,hnRNP_A3_Positive==TRUE & group=='CCP'))
nrow(subset(experimental_merge,hnRNP_A3_Positive==TRUE & group=='CCP_Acet'))
nrow(subset(experimental_merge,hnRNP_A3_Positive==TRUE & group=='CCP_Carb_Acet'))
nrow(subset(experimental_merge,hnRNP_A3_Positive==TRUE & group=='Negative'))

#collective overview table
data_collective=read.csv("PositiveData_Groupwise.csv",header = TRUE,row.names = 1)

#####################################################################################
#donut chart CCP 
data <- data.frame(
        category=c("Fillagrin", "Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        count=c(9,16,30,13,23))
data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$category, "\n value: ", data$count)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
        geom_rect() +
        geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
        scale_fill_brewer(palette=4) +
        coord_polar(theta="y") +
        xlim(c(1, 4)) +
        theme_void() +
        theme(legend.position = "none")

#donut chart CCP_Acet
data <- data.frame(
        category=c("Fillagrin", "Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        count=c(6,10,15,7,14))
data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$category, "\n value: ", data$count)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
        geom_rect() +
        geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
        scale_fill_brewer(palette=5) +
        coord_polar(theta="y") +
        xlim(c(1, 4)) +
        theme_void() +
        theme(legend.position = "none")

#donut chart Triple positive
data <- data.frame(
        category=c("Fillagrin", "Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        count=c(49,54,64,53,62))
data$fraction = data$count / sum(data$count)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$category, "\n value: ", data$count)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
        geom_rect() +
        geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
        scale_fill_brewer(palette=8) +
        coord_polar(theta="y") +
        xlim(c(1, 4)) +
        theme_void() +
        theme(legend.position = "none")

########################################################################################
##vein diagram

##CCP 
##checking combination 
nrow(subset(experimental_merge, Filaggrin_Positive==TRUE & Vimentin_positive==TRUE & Fibrinogen_positive==TRUE & 
                    alpha_enolase_Positive==TRUE & hnRNP_A3_Positive==TRUE & group=='CCP'))
library(grid)
library(VennDiagram)
grid.newpage()
venn.plot <- draw.quintuple.venn(
        area1 = 9,
        area2 = 16,
        area3 = 30,
        area4 = 13,
        area5 = 23,
        n12 = 4,
        n13 = 9,
        n14 = 5,
        n15 = 7,
        n23 = 16,
        n24 = 7,
        n25 = 14,
        n34 = 12,
        n35 = 23,
        n45 = 10,
        n123 = 4,
        n124 = 2,
        n125 = 4,
        n134 = 5,
        n135 = 7,
        n145 = 4,
        n234 = 7,
        n235 = 14,
        n245 = 7,
        n345 = 10,
        n1234 = 2,
        n1235 = 4,
        n1245 = 2,
        n1345 = 4,
        n2345 = 7,
        n12345 = 2,
        category = c("Filaggrin", "Vimentin", "Fibrinogen", "alpha_enolase", "hnRNP_A3"),
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 2,
        margin = 0.05,
        ind = TRUE
)


##Triple positive
grid.newpage()
venn.plot <- draw.quintuple.venn(
        area1 = 49,
        area2 = 54,
        area3 = 64,
        area4 = 53,
        area5 = 62,
        n12 = 42,
        n13 = 49,
        n14 = 42,
        n15 = 48,
        n23 = 54,
        n24 = 48,
        n25 = 52,
        n34 = 53,
        n35 = 61,
        n45 = 51,
        n123 = 42,
        n124 = 38,
        n125 = 41,
        n134 = 42,
        n135 = 48,
        n145 = 41,
        n234 = 48,
        n235 = 52,
        n245 = 46,
        n345 = 51,
        n1234 = 38,
        n1235 = 41,
        n1245 = 37,
        n1345 = 41,
        n2345 = 46,
        n12345 = 37,
        category = c("Filaggrin", "Vimentin", "Fibrinogen", "alpha_enolase", "hnRNP_A3"),
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 2,
        margin = 0.05,
        ind = TRUE
)

## CCP_Acet

grid.newpage()
venn.plot <- draw.quintuple.venn(
        area1 = 6,
        area2 = 10,
        area3 = 15,
        area4 = 7,
        area5 = 14,
        n12 = 5,
        n13 = 6,
        n14 = 5,
        n15 = 6,
        n23 = 9,
        n24 = 5,
        n25 = 9,
        n34 = 7,
        n35 = 13,
        n45 = 7,
        n123 = 5,
        n124 = 4,
        n125 = 5,
        n134 = 5,
        n135 = 6,
        n145 = 5,
        n234 = 5,
        n235 = 9,
        n245 = 5,
        n345 = 7,
        n1234 = 4,
        n1235 = 5,
        n1245 = 4,
        n1345 = 5,
        n2345 = 5,
        n12345 = 4,
        category = c("Filaggrin", "Vimentin", "Fibrinogen", "alpha_enolase", "hnRNP_A3"),
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 2,
        margin = 0.05,
        ind = TRUE
)

## allgroups together

library(grid)
library(VennDiagram)
grid.newpage()
venn.plot <- draw.quintuple.venn(
        area1 = 64,
        area2 = 81,
        area3 = 115,
        area4 = 74,
        area5 = 100,
        n12 =51,
        n13 = 64,
        n14 = 52,
        n15 = 61,
        n23 = 79,
        n24 = 60,
        n25 = 75,
        n34 = 73,
        n35 = 98,
        n45 = 68,
        n123 = 51,
        n124 = 44,
        n125 = 50,
        n134 = 52,
        n135 = 61,
        n145 = 50,
        n234 = 60,
        n235 = 75,
        n245 = 58,
        n345 = 68,
        n1234 = 44,
        n1235 = 50,
        n1245 = 43,
        n1345 = 50,
        n2345 = 58,
        n12345 = 43,
        category = c("Filaggrin", "Vimentin", "Fibrinogen", "alpha_enolase", "hnRNP_A3"),
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 2,
        margin = 0.05,
        ind = TRUE
)

######################################################################################
#pie-chart
library(ggplot2)
library(tidyverse)

# Create Data
data <- data.frame(
        group=c("Filaggrin","Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        value=c(9,16,30,13,23)
)
# the position of labels
data <- data %>% 
        arrange(desc(group)) %>%
        mutate(prop = value / sum(data$value) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )

#pie chart
ggplot(data, aes(x="", y=prop, fill=group)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
        theme_void() + 
        theme(legend.position="none") +
        
        geom_text(aes(y = ypos, label = group), color = "white", size=5) +
        scale_fill_brewer(palette="Set1")

##CCP_Acet
data <- data.frame(
        group=c("Filaggrin","Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        value=c(6,10,15,7,14)
)

# the position of labels
data <- data %>% 
        arrange(desc(group)) %>%
        mutate(prop = value / sum(data$value) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )

#pie chart
ggplot(data, aes(x="", y=prop, fill=group)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
        theme_void() + 
        theme(legend.position="none") +
        
        geom_text(aes(y = ypos, label = group), color = "white", size=5) +
        scale_fill_brewer(palette="Set1")

##triple positive pie chart
data <- data.frame(
        group=c("Filaggrin","Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        value=c(49,54,64,53,62)
)

# the position of labels
data <- data %>% 
        arrange(desc(group)) %>%
        mutate(prop = value / sum(data$value) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )

#pie chart
ggplot(data, aes(x="", y=prop, fill=group)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
        theme_void() + 
        theme(legend.position="none") +
        
        geom_text(aes(y = ypos, label = group), color = "white", size=5) +
        scale_fill_brewer(palette="Set1")

###############################################################################
# circular bar plot
library(tidyverse)

# Create dataset
data <- data.frame(
        individual=c("Filaggrin","Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        group=c( rep('CCP', 5), rep('CCP_Acet', 5), rep('CCP_Carb_Acet', 5), rep('Negative', 5)) ,
        value=c(9,16,30,13,23,6,10,15,7,14,49,54,64,53,62,0,1,6,1,1 )
)

#using percentage
data <- data.frame(
        individual=c("Filaggrin","Vimentin","Fibrinogen","alpha.enolase","hnRNP_A3"),
        group=c( rep('CCP', 5), rep('CCP_Acet', 5), rep('CCP_Carb_Acet', 5), rep('Negative', 5)) ,
        value=c(25,44.4,83.3,36.1,63.8,33.3,55.5,83.3,38.8,77.7,73.13,80.59,95.52,79.1,92.53,0,3.84,23.07,3.84,3.84 )
)
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
        group_by(group) %>% 
        summarize(start=min(id), end=max(id) - empty_bar) %>% 
        rowwise() %>% 
        mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +      
        
        geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
        
       
        geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        
        
        annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
        
        geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
        ylim(-100,120) +
        theme_minimal() +
        theme(
                legend.position = "none",
                axis.text = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                plot.margin = unit(rep(-1,4), "cm") 
        ) +
        coord_polar() + 
        geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
        
        # Add base line information
        geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
        geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p
###########################################################################################################
## odds ratio calculation
present <- c("No", "Yes")
peptide=c("Filaggrin","Vimentin","Fibrinogen","alpha_enolase","hnRNP_A3")
library(epitools)

##CCP
nrow(subset(experimental_merge,group=='CCP'))
dat <- matrix(c(27,9,20,16,6,30,23,13,13,23), nrow = 5, ncol = 2, byrow = TRUE)
dimnames(dat) <- list("ACPA" = peptide, "AE Present" = present)
or_fit <- oddsratio(dat)

##CCP_Acet
nrow(subset(experimental_merge,group=='CCP_Acet'))
dat <- matrix(c(12,6,8,10,3,15,11,7,4,14), nrow = 5, ncol = 2, byrow = TRUE)
dimnames(dat) <- list("ACPA" = peptide, "AE Present" = present)
or_fit <- oddsratio(dat)

##triple positive
nrow(subset(experimental_merge,group=='CCP_Carb_Acet'))
dat <- matrix(c(18,49,13,54,3,64,14,53,5,62), nrow = 5, ncol = 2, byrow = TRUE)
dimnames(dat) <- list("ACPA" = peptide, "AE Present" = present)
or_fit <- oddsratio(dat)

#########################################################################################################
##CLUSTERING
library(factoextra)

######### kmer mean clustering

df <-cit_dataCopy
df=t(df)
df <- scale(df)
head(df)
#distance matrix calculation and visualisation
distance <- get_dist(df)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#clustering
k2 <- kmeans(df, centers = 2, nstart = 25)
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)

#plotting the clusters
p1=fviz_cluster(k4, geom = "point", data = df) + ggtitle("k = 4")
p2=fviz_cluster(k5, geom = "point", data = df) + ggtitle("k = 5")
p3=fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p4=fviz_cluster(k3, geom = "point", data = df) + ggtitle("k = 3")
library(gridExtra)
grid.arrange( p3,p4,p1,p2, nrow = 2)

# elbow method
fviz_nbclust(df, kmeans, method = "wss") +
        geom_vline(xintercept = 4, linetype = 2)+
        labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")

#gap statistics method
set.seed(123)
library(cluster)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,K.max = 10, B = 500)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

#kmean clustering
set.seed(123)
final4 <- kmeans(df, 4, nstart = 25)
print(final4)
print(km.res)
fviz_cluster(final4, data = df)

################ Hierachial clustering
library(tidyverse)  
library(cluster)   
library(factoextra) 
library(dendextend) 

df <-cit_dataCopy
df=t(df)
df <- scale(df)
head(df)
# Dissimilarity matrix
d <- dist(df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

#######Agglomerative hierarchial clustering(AGNES)
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {agnes(df, method = x)$ac}
map_dbl(m, ac)
hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
cutree(as.hclust(hc3), k = 4)

##### divisive hierarchical clustering(DIANA)
hc4 <- diana(df)
# Divise coefficient; amount of clustering structure found
hc4$dc
# plot dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")
cutree(as.hclust(hc4), k = 4)
plot(hc4, cex = 0.6)
rect.hclust(hc4, k = 4, border = c(2:5))

##########
# Ward's method
hc5 <- hclust(d, method = "ward.D2" )

# Cut tree into 4 groups
sub_grp <- cutree(hc5, k = 4)

# Number of members in each cluster
table(sub_grp)
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 4, border = c(2:5))
fviz_cluster(list(data = df, cluster = sub_grp))

# Compute 2 hierarchical clusterings
hc1 <- hclust(d, method = "complete")
hc2 <- hclust(d, method = "ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2)

dend_list <- dendlist(dend1, dend2)

#calculating entanglement to see the quality in alignment in tree(entaglement scores comes as 0.2 so alignment is good )
tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, 
           common_subtrees_color_lines = FALSE, 
           common_subtrees_color_branches = TRUE,  
           main = paste("entanglement =", round(entanglement(dend_list), 2)))

##latent class analysis
library(poLCA)

## replacing all negatives with 2 and positives with 3 as 1 and 2 is throwing an error
d=cit_dataCopy
d1 <- apply(d, 2, function(x) ifelse(x > 0, 2, x))
d1 <- apply(d1, 2, function(x) ifelse(x < 0, 1, x))
data=d1+1;
data=as.data.frame(data)

#formula
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~1

#run a sequence of models with 1-10 classes and print out the model with the lowest BIC
max_II <- -100000
min_bic <- 100000
for(i in 2:10){
        lc <- poLCA(f, data, nclass=i, maxiter=3000, 
                    tol=1e-5, na.rm=FALSE,  
                    nrep=10, verbose=TRUE, calc.se=TRUE)
        if(lc$bic < min_bic){
                min_bic <- lc$bic
                LCA_best_model<-lc
        }
}    	
LCA_best_model
LCA_final_model=poLCA(f, data, nclass=2, maxiter=3000, 
                      tol=1e-5, na.rm=FALSE,  
                      nrep=10, verbose=TRUE, calc.se=TRUE,graphs = TRUE)

##some stats aboutour best fit model
table(LCA_final_model$predclass)
round(prop.table(table(LCA_final_model$predclass)),4)*100

# covariates
r1=read.csv("radio_merge.csv",header=TRUE,row.names = 1)
setwd("~/Desktop/RA")
t1=read.csv("Experiment_info.csv",header = TRUE,row.names = 1)
merged_final=merge(data,r1,by=0)
row.names(merged_final)=merged_final$Row.names
merged_final=merged_final[,c(2:26)]
merged_final=merge(merged_final,t1,by=0)
row.names(merged_final)=merged_final$Row.names
merged_final=merged_final[,c(2:30)]

#total Progression as covariate
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Progression_total.x
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 
nrow(subset(merged_final,Progression_total.x==0))

#only joint erosion as covariate
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Progression_erosions.x
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE)
nrow(subset(merged_final,Progression_erosions.x==0)) ## total 74 with no progression out of 92

#only joint progression
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Progression_joint.x
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 
nrow(subset(merged_final,Progression_joint.x==0))

#gender as covariate
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Sex
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 
nrow(subset(merged_final,Sex=='Female'))

#gender and progression
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Sex+Progression_total.x
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 

#AlcoholInatke and progression
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~AlcoholIntake+Progression_total.x
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 

# Gender and AlcoholInatke
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Sex+AlcoholIntake
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 

# Gender and AlcoholInatke and Progression
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~Sex+Progression_total.x+AlcoholIntake+AlcoholIntake
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 

# alcoholInatke
f <- cbind(CCP_1.cit,Vim60_75.cit,Vim2_17.cit,Fib36_52.cit,Fib573.cit,Fib591.cit,Fib.alpha621_635.cit,Fib.Alpha36_50.cit,Fib.beta60_74.cit,CEP_1,Pept.1,Pept_Z1,Pept_Z2,Pept.5,Bla26)~AlcoholIntake
ch2 <- poLCA(f,merged_final,nclass=2,graphs = TRUE) 

##############################################################################################################################
#MPlus Latent class analysis
setwd("~/Desktop/RA/ACPA/MPlus")
library(MplusAutomation)
prepareMplusData(experimental_data2, "demo.dat")
#createModels("t1.txt")
runModels(target = getwd())
output=readModels(target = getwd(), recursive = FALSE, what = c("tech1","tech8"), quiet = FALSE)


##############################################################################################################################
library(poLCA)
d=experimental_data2
data=as.data.frame(data)
f <- cbind(Filaggrin_Positive,Vimentin_positive,Fibrinogen_positive,alpha_enolase_Positive,hnRNP_A3_Positive)~1
max_II <- -100000
min_bic <- 100000
for(i in 2:10){
        lc <- poLCA(f, data, nclass=i, maxiter=3000, 
                    tol=1e-5, na.rm=FALSE,  
                    nrep=10, verbose=TRUE, calc.se=TRUE)
        if(lc$bic < min_bic){
                min_bic <- lc$bic
                LCA_best_model<-lc
        }
} 
LCA_final_model=poLCA(f, data, nclass=2, maxiter=3000, 
                      tol=1e-5, na.rm=FALSE,  
                      nrep=10, verbose=TRUE, calc.se=TRUE,graphs = TRUE)


merged_final=merge(data,r1,by=0)
row.names(merged_final)=merged_final$Row.names
merged_final=merged_final[,c(2:16)]
merged_final=merge(merged_final,t1,by=0)
row.names(merged_final)=merged_final$Row.names
merged_final=merged_final[,c(2:21)]

#gender as covariate
f <- cbind(Filaggrin_Positive,Vimentin_positive,Fibrinogen_positive,alpha_enolase_Positive,hnRNP_A3_Positive)~Sex


##############################################################################################################################################
setwd("~/Desktop/RA/ACPA/MPlus")
disease_activity=read.csv("experiment_info_updated_20200703.csv",header = TRUE,row.names = 1)
disease_activity1=disease_activity[,c(1,3,4,15,16)]
experimental_data3=merge(experimental_data2,disease_activity1,by=0)
row.names(experimental_data3)=experimental_data3$Row.names
experimental_data3=experimental_data3[,c(2:11)]

#MPlus
setwd("~/Desktop/RA/ACPA/MPlus")
library(MplusAutomation)
prepareMplusData(experimental_data3, "demo.dat")
runModels(target = getwd())

######################################################################################################
#zero inflated binomial regression
library(pscl) 
library(emmeans) 


#RA_Activity and 5 groups 
z1=read.csv("z1.csv")
z1=z1[,c(2:7)]
z_1 = zeroinfl(RA_Activity ~ Filaggrin_Positive +Vimentin_positive+Fibrinogen_positive+alpha_enolase_Positive+hnRNP_A3_Positive , dist="negbin", data=z1)
summary(z_1)

#Progression with 5 groups
z2=read.csv("z2.csv")
z2=z2[,c(2:7)]
z_2 = zeroinfl(group ~ Filaggrin_Positive +Vimentin_positive+Fibrinogen_positive+alpha_enolase_Positive+hnRNP_A3_Positive , dist="negbin", data=z2)
summary(z_2)

#########################################################################################
#RA_Activity and 2 classes (logistic regression)
class=read.csv("class_RA_Activity.csv",header=1)
z3=class[,c(6,9)]
z3[z3$RA_Activity==1,]$RA_Activity="H"
z3[z3$RA_Activity==2,]$RA_Activity="L"
z3[z3$RA_Activity==3,]$RA_Activity="M"
z3[z3$RA_Activity==4,]$RA_Activity="R"
xtabs(~RA_Activity+Latent.Class,data=z3)


