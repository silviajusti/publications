#Libraries used
library(ape)
library(dplyr)
library(reshape2)
library(tidyverse)
library(data.table)
library(plyr)


#Read COI (or any other marker...) alignment and calculate pairwise distances

data <- read.FASTA("BOLD_noPERU_alignnment.fasta")
dist_data <- as.data.frame(dist.dna(data, model = "K80", 
                    as.matrix = TRUE, pairwise.deletion = T))
d<-dist_data #this "d" will be used later
write.csv(dist_data, file="dist_data.csv") #save pairwise distances to file

#Format the table to a 3 column pairwise comparison, and keep only query sample name and calculated distance value.
df<-reshape2::melt(dist_data)

#Sort the pairwise distance values select only the valuse between 0.01 and 0.05. Because 0.01 is a rough intraspecific distance and 0.05 is usually interspecific, for this marker for insects. 

sorted_dists<-as.data.table(df[order(df$value, decreasing = FALSE),])
dists_filtered<-sorted_dists[sorted_dists$value <= 0.05 ] 
dists_filtered2<-dists_filtered[dists_filtered$value > 0.01 ]

#Take a quick look ate the values filtered for the intra-inter specific distances region.

plot(dists_filtered2$value)

#Now we find the first barcode gap, by subtracting the distance value, by the previous distance value successively
dt <- as.data.table(dists_filtered2)
setkey(dt, value)
dt[, diff := value - shift(value, fill = first(value))]


#And we plot the differences to take a look

plot(dt$diff)

#Check the aproximate distance of the first observed gap, and round that number down. E.g. if the first observed gap is around y=0.0015, use 0.001

dt_high<-dt[dt$diff > 0.0010 ] 
dt_high #look at the valuse

#Now, find the max difference, which will be the correspondent to higher end of the barcode gap

max(dt_high$diff)

#Now, find the sample that generated the gap

t<-dt$diff < max(dt_high$diff)
min(which(t == FALSE))
row_id<-(min(which(t == FALSE)))
dt[(row_id-1):(row_id+1),]

#Find the gap using the values calculated above

intra<-as.matrix(dist_data)
inter<-as.matrix(dist_data)
intra[intra>min((dt[(row_id-1):(row_id+1),]) %>% select (2))] <- NA; intra  
inter[inter<max((dt[(row_id-1):(row_id+1),]) %>% select (2))] <- NA; inter  
df_final <- data.frame(melt(intra), melt(inter)) %>% select(1, 3,6)
names(df_final) <- c("lineage", "intra_ID", "inter_ID")


#Plot the gap with all samples

boxplot(df_final$inter_ID ~ df_final$lineage,  ylim=c(0,0.15),  
        ylab = "K2p pairwise distance", border= "gray31", xaxt = "n", xlab = "",
        col="gray47", pch=21, frame=F, cex=0.5)

boxplot(df_final$intra_ID ~ df_final$lineage, add=TRUE, ylim=c(0,0.15),  
        ylab = "K2p pairwise distance", border= "gray61", xaxt = "n", xlab = "",
        col="gray90", pch=21, frame=F, cex=0.5)
abline(h=min((dt[(row_id-1):(row_id+1),]) %>% select (2)), col="green",lty=3)
abline(h=max ((dt[(row_id-1):(row_id+1),]) %>% select (2)), col="red",lty=3)


#And the barcode gap is...

min((dt[(row_id-1):(row_id+1),]) %>% select (2))
max ((dt[(row_id-1):(row_id+1),]) %>% select (2))

#Write the files with the  comparisons
groups<-reshape2::melt(intra)
write.csv(groups, file="intra_groups.csv")
groups_inter<-reshape2::melt(inter)
write.csv(groups_inter, file="inter_groups.csv")

#Write the groups and plot the clusters

dist_data<-as.dist(dist_data)
hc<-hclust(dist_data,"complete")
plot(hclust(dist_data),cex=0.5)
abline(h= round(min((dt[(row_id-1):(row_id+1),]) %>% 
                      select (2)),digits=7), col="red",lty=3)
final_clusters<-as.data.frame(cutree(hc, h =  round(min((dt[(row_id-1):(row_id+1),]) %>% select(2)),digits=7)))
names(final_clusters) <- c("lineage")
final_clusters
write.csv(final_clusters,file="final_clusters.csv")

#Look at the csv. Find out how many cluseters and what you are anming them. In the example below, I have 4 clusters and I am naming them according to their taxonomy.

clusters<-with(data.frame(final_clusters), replace(final_clusters, final_clusters == 1, "Ch. bonnae"))
clusters<-with(data.frame(clusters), replace(clusters, clusters == 2, "Ch. avielae"))
clusters<-with(data.frame(clusters), replace(clusters, clusters == 3, "Ch. bathana"))
clusters<-with(data.frame(clusters), replace(clusters, clusters == 4, "Chagasia sp"))
clusters #check the clusters
```

#Now format the dataframe so that, instead of the name of the sample, we have the name of the cluster, so we can calculate diversity within and between clusters

names <- rownames(d)
rownames(d) <- NULL
data <- cbind(names,d)
dist_cluster<-cbind(clusters,data[,2:ncol(data)])
e <- dist_cluster
names <- colnames(e)
colnames(e) <- NULL
data <- t(rbind(names,e))
colnames(data) <- data[1,]
data = data[-1,]
dist_cluster<-cbind(clusters$lineage,data)
dist_cluster<-cbind(clusters$lineage,data[,-1])
dist_cluster2 <- dist_cluster[,-1]
rownames(dist_cluster2) <-dist_cluster[,1]
molten_clusters<-reshape2::melt(dist_cluster2)
comps<-molten_clusters


#Writing out the comparisons and then calculating means per group

write.table(molten_clusters,file = "molten_comps.txt", sep="\t")
comps<-read.table("molten_comps.txt", header=TRUE, sep="\t")
names(comps) <- c("var1", "var2", "dist")
SD_comps<-reshape2::dcast(comps, var1~var2,sd)
mean_comps<-reshape2::dcast(comps, var1~var2,mean)
mean_and_SD<-cbind((reshape2::melt(mean_comps)),(reshape2::melt(SD_comps))[,3])
names(mean_and_SD) <- c("Lineage 1", "Lineage 2", "Mean distance", "Standard deviation")
mean_and_SD
write.table(mean_and_SD,file = "mean_and_sd.txt", sep="\t")

