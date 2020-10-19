setwd("~/Desktop/RA/ACPA")
#citrullinated data only
cit_data=read.csv("cit_updated.csv",row.names = 1,header = TRUE)

#RA_patients
RA_Patients=read.csv("experiment_info_updated_20200723.csv",row.names = 1,header=TRUE)
RA_Patients=subset(RA_Patients,Diagnosis=='RA (ACR-EULAR 2010)')

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
test=merge1[,c(2:33)]

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

###########################################################################################
# Hierachial clustering
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

###########################################################################################
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

#Positive in each group
experimental_data=cit_dataCopy
experimental_data$Filaggrin_Positive=as.factor(experimental_data$CCP_1.cit>0)
experimental_data$Vimentin_positive=as.factor(experimental_data$Vim60_75.cit>0 | experimental_data$Vim2_17.cit>0)
experimental_data$Fibrinogen_positive=as.factor(experimental_data$Fib36_52.cit>0 | experimental_data$Fib573.cit>0 | 
                                                  experimental_data$Fib591.cit>0|experimental_data$Fib.Alpha36_50.cit>0|experimental_data$Fib.alpha621_635.cit>0|experimental_data$Fib.beta60_74.cit>0)
experimental_data$alpha_enolase_Positive=as.factor(experimental_data$CEP_1>0)
experimental_data$hnRNP_A3_Positive=as.factor(experimental_data$Pept.1>0|
                                                experimental_data$Pept_Z1>0|experimental_data$Pept_Z2>0|experimental_data$Pept.5>0|experimental_data$Bla26>0)
experimental_data=experimental_data[,c(16:20)]

#Preparing data for mplus
experimental_data2=merge(test,experimental_data,by=0)
row.names(experimental_data2)=experimental_data2$Row.names
experimental_data2=experimental_data2[,c(34:38,23:33)]
experimental_diseaseActivity=experimental_data2[,c(1:5,15,16)]
sapply(experimental_diseaseActivity, function(x) sum(is.na(x)))
experimental_diseaseActivity=na.omit(experimental_diseaseActivity)
experimental_dataProgression=experimental_data2[,c(1:5,14)]
sapply(experimental_dataProgression, function(x) sum(is.na(x)))
experimental_dataProgression=na.omit(experimental_dataProgression)

#Latent class analysis in MPLUS

#disease activity data
setwd("~/Desktop/RA/ACPA/MPlus")
library(MplusAutomation)
prepareMplusData(experimental_diseaseActivity, "demo.dat")
runModels(target = getwd())

#Progression
library(MplusAutomation)
prepareMplusData(experimental_dataProgression, "demo.dat")
runModels(target = getwd())
