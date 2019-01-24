setwd("C:/Users/My_Data/")
library("flowCore")
library(cytofCore)
library(cytofkit)
library(cluster)
to_split_on<-"_np_"
set.seed(123456)
maindir<-"C:/Users/My_Data/results"
dir <- getwd()
files <- list.files(dir ,pattern='.csv$', full=F)

for (i in files) {                                                                             
datamatrix<-read.csv(i)
print(i)                                                       
setwd("./results")
name<-strsplit(i, ".csv")[[1]]
#name<-unlist(name)
                                                                            
cytofCore.write.FCS(as.matrix(datamatrix), filename=paste( "em_", name,".fcs", sep=""), what = "numeric")     
setwd("./..") } 

setwd("./results")
dir <- getwd()
filenames <- list.files(dir ,pattern='.fcs$', full=F)
combined_data_transformed <- cytof_exprsMerge(fcsFiles = filenames[1:3], comp=FALSE,                                             
                                              transformMethod = "none",
                                              mergeMethod = "all")
                                              
cytofCore.write.FCS(combined_data_transformed, filename="simple_combined_data_transformed.fcs", what = "numeric")


#remove unwanted markers 
data_transformed<-combined_data_transformed[,-c(1:6,13,21)] 


print("step tsne")
data_transformed_tsne <- cytof_dimReduction(data=data_transformed, method = "tsne")


 print("step 30")
cluster_PhenoGraph_30 <- Rphenograph(data_transformed, k=30 )
cluster_PhenoGraph_k30<-cluster_PhenoGraph_30$membership

data_all_k30 <- cbind(data_transformed, data_transformed_tsne, PhenoGraph = cluster_PhenoGraph_k30)


tsne_30<-cytof_clusterPlot(data=data_all_k30, cluster="PhenoGraph",  xlab="tsne_1", ylab="tsne_2")

pdf("cluster_tsne_phenograph_bexp_cct_my_test_30.pdf")

tsne_30
dev.off()

print("step 45")
cluster_PhenoGraph_45 <- Rphenograph(data_transformed, k=45 )
cluster_PhenoGraph_k45<-cluster_PhenoGraph_45$membership
data_all_k45 <- cbind(data_transformed, data_transformed_tsne, PhenoGraph = cluster_PhenoGraph_k45)
tsne_45<-cytof_clusterPlot(data=data_all_k45, cluster="PhenoGraph",  xlab="tsne_1", ylab="tsne_2")
pdf("cluster_tsne_phenograph_bexp_cct_my_test_45_2.pdf")
tsne_45
dev.off()


print("step 60")
cluster_PhenoGraph_60 <- Rphenograph(data_transformed, k=60 )
cluster_PhenoGraph_k60<-cluster_PhenoGraph_60$membership
data_all_k60 <- cbind(data_transformed, data_transformed_tsne, PhenoGraph = cluster_PhenoGraph_k60)

tsne_60<-cytof_clusterPlot(data=data_all_k60, cluster="PhenoGraph",  xlab="tsne_1", ylab="tsne_2")
 pdf("cluster_tsne_phenograph_bexp_cct_my_test_60.pdf")
tsne_60
dev.off()


print("step 75")
cluster_PhenoGraph_75 <- Rphenograph(data_transformed, k=75 )
cluster_PhenoGraph_k75<-cluster_PhenoGraph_75$membership
data_all_k75 <- cbind(data_transformed, data_transformed_tsne, PhenoGraph =cluster_PhenoGraph_k75)

tsne_75<-cytof_clusterPlot(data=data_all_k75, cluster="PhenoGraph",  xlab="tsne_1", ylab="tsne_2")
pdf("cluster_tsne_phenograph_bexp_cct_my_test_75.pdf")
tsne_75
dev.off()


#########################################Concatenate#######################################


env <-ls()
lista_data_all<- env[grep("data_all_k", env)]
for (b in lista_data_all){
data_all_chnls<-cbind(combined_data_transformed[,c(1:6,13,21)] ,get(b) )

#pdf(paste(b,"cluster_tsne_phenograph_bexp_cct_all_chnls.pdf",sep=""))
#cytof_clusterPlot(data=data_all_chnls, cluster="PhenoGraph",  xlab="tsne_1", ylab="tsne_2")
#dev.off()


dir.create(b)
setwd(file.path(maindir, b))

dir.create("FCS_per_cluster")

setwd("./FCS_per_cluster")

data_all_chnls_fr<-as.data.frame(data_all_chnls)
lista_cluster<-unique(data_all_chnls_fr$PhenoGraph)


for (i in lista_cluster) {
sub_data<- data_all_chnls_fr[data_all_chnls_fr$PhenoGraph==i,]
cytofCore.write.FCS(as.matrix(sub_data), filename=paste("clusterN_",i,".fcs",sep=""), what = "numeric")
}




cytofCore.write.FCS(data_all_chnls, filename="combined_data_transformed_wTnsePheno.fcs", what = "numeric")

dir.create(b)
setwd(file.path(maindir, b))

dir.create("FCS_per_sample")


setwd("./FCS_per_sample")
 
 lista_samples<-rownames(data_all_chnls_fr)
 
 
 sampleID<-NULL

for (h in lista_samples) {
sample<-unlist(strsplit(h, to_split_on))[1]
sampleID<-c(sampleID,sample)
print("check 1")
}


data_all_chnls_fr_samples<-cbind(data_all_chnls_fr,sampleID)


samples<-unique(sampleID)







Tot_counts<-c(1: length(lista_cluster))
Tot_percentages<-c(1: length(lista_cluster))
order<-as.character(1: length(lista_cluster))

for (v in samples) {
sub_data_samples<- data_all_chnls_fr[data_all_chnls_fr_samples$sampleID==v,]
print("check 2")
count_cluster<-as.matrix(table(sub_data_samples$PhenoGraph))
colnames(count_cluster)<-v

percentages<-prop.table(count_cluster)*100
colnames(percentages)<-v





cytofCore.write.FCS(as.matrix(sub_data_samples), filename=paste("sample_",v,".fcs",sep=""), what = "numeric")

#for(z in 1:length(rownames(percentages))){
#if(rownames(percentages)[z])
#}


Tot_percentages<-merge(Tot_percentages,percentages,  by="row.names",all.y=T,all.x=T , sort=F)

rownames(Tot_percentages)<-Tot_percentages[,1]
Tot_percentages_order<-Tot_percentages[order, ]
Tot_percentages<-Tot_percentages_order[,-1]                       
print("check 3")
Tot_counts<-merge(Tot_counts,count_cluster,  by="row.names", all.y=T,all.x=T, sort=T)
rownames(Tot_counts)<-Tot_counts[,1]
Tot_counts_order<-Tot_counts[order, ]
Tot_counts<-Tot_counts_order[,-1]
print("check 4")
#Tot_counts<-cbind(Tot_counts, count_cluster)
#Tot_percentages<-cbind(Tot_percentages,percentages)

}
write.table(Tot_percentages,"Tot_percentages.txt",sep="\t")
write.table(Tot_counts,"Tot_counts.txt",sep="\t")
setwd(maindir) }

 save(list=ls(),file="workspace_stepFinal.rda")
