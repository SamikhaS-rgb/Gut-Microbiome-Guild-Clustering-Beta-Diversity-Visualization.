install.packages("plyr")
install.packages("FSA")
install.packages("patchwork")
install.packages("reshape2")
install.packages("pheatmap")

library(ggplot2)
library(ggpubr)
library(vegan)
library(ape)
library(rcompanion)
library(multcompView)
library(multcomp)
library(plyr)
library(stats)
library(FSA)
library(patchwork)
library(readxl)
library(reshape2)
library(pheatmap)


##Permanova##
pairwise.adonis <- function(x, factors, p.adjust.m) {
  library(vegan)
  co <- as.matrix(combn(unique(factors), 2))
  pairs <- c()
  F.Model <- c()  # Fixed variable name here
  R2 <- c()
  p.value <- c()
  for (elem in 1:ncol(co)) {
    ind1 <- which(factors == co[1, elem])
    ind2 <- which(factors == co[2, elem])
    m <- x[c(ind1, ind2), c(ind1, ind2)]
    group <- c(factors[ind1], factors[ind2])
    m <- as.dist(m)
    group <- as.data.frame(group)
    set.seed(33)
    ad <- adonis2(m ~ group, data = group, permutations = 999)
    
    pairs <- c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
    F.Model <- c(F.Model, ad$F[1])  # Ensure correct variable is used
    R2 <- c(R2, ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }
  
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  return(pairw.res)
}

##repeated measurement correlation##
CROS=function(mydata)
{
  result=c()
  tmpfit=lm(var1~sub+var2,data=mydata)
  tmpanova=anova(tmpfit)
  if(!is.na(tmpfit$coefficients["var2"]))
  {
    if(tmpfit$coefficients["var2"]==0)
    {
      result[1]=0
    }else if(tmpfit$coefficients["var2"]<0)
    {
      result[1]=-sqrt(tmpanova[2,2]/(tmpanova[2,2]+tmpanova[3,2]))
    }else
    {
      result[1]=sqrt(tmpanova[2,2]/(tmpanova[2,2]+tmpanova[3,2]))
    }
    result[2]=tmpanova[2,5]
  }else{
    result[1]=0
    result[2]=1
  }
  return(result)
}



getwd() 
setwd("/Users/samikhasrinivasan/Desktop/Microbiome Data Analysis")
myallgroup=read.table('mydata.R.StorchLab.03292024.txt',header = T,row.names = 1,sep = '\t')
#add another column 
colnames(myallgroup)[1]="Genotypeold"
myallgroup$Genotype=plyr::mapvalues(myallgroup$Genotypeold,
                                    from=c("flox", "liver", "int"), to=c("WildType","LFABPliverKO","LFABPIntKO"))
myallgroup$group=paste(myallgroup$Genotype,myallgroup$Time_Point,sep = "_")
myallgroup$Time_Point=factor(myallgroup$Time_Point, levels = c("T1","T2"))
myallgroup_time_color=c("T1"="#327DBE","T2"="#D76464")
myallgroup$Genotype=factor(myallgroup$Genotype,levels = c("WildType","LFABPliverKO","LFABPIntKO"))
myallgroup_genotype_color=c("WildType"="#030303", "LFABPliverKO"="#FF8C00","LFABPIntKO"="#0509fa")
myallgroup$group=factor(myallgroup$group,
                        levels=c("WildType_T1","WildType_T2","LFABPliverKO_T1","LFABPliverKO_T2","LFABPIntKO_T1","LFABPIntKO_T2"))
#myallgroup_genotypetp_color=c("#030303","#030303","#FF8C00","#FF8C00","#0509fa","#0509fa")
myallgroup_genotypetp_color = c("#030303", "#D0D0D0", "#0509FA", "#B0E0E6", "#FF8C00", "#ffebcc")

names(myallgroup_genotypetp_color)=levels(myallgroup$group)
head(myallgroup)

myallbraycurtis=read.table('bray_curtis_distance_matrix.txt',header = T,row.names = 1,sep = '\t', check.names=F)
myallbraycurtis=myallbraycurtis[rownames(myallgroup), rownames(myallgroup)]
myalljaccard=read.table('jaccard_distance_matrix.txt',header = T,row.names = 1,sep = '\t',check.names=F)
myalljaccard=myalljaccard[rownames(myallgroup), rownames(myallgroup)]
myallunweiunifrac=read.table('unweighted_unifrac_distance_matrix.txt',header = T,row.names = 1,sep = '\t',check.names=F)
myallunweiunifrac=myallunweiunifrac[rownames(myallgroup), rownames(myallgroup)]
myallweiunifrac=read.table('weighted_unifrac_distance_matrix.txt',header = T,row.names = 1,sep = '\t',check.names=F)
myallweiunifrac=myallweiunifrac[rownames(myallgroup), rownames(myallgroup)]
rownames(myallbraycurtis)==rownames(myallgroup)

myallabun=data.frame(read_xlsx('Abundance.xlsx', sheet='13000depth',range='B2:BU1180'),check.names=F, row.names=1)
#change the row and column
myallabun=data.frame(t(myallabun),check.names=F)
apply(myallabun,1,sum)
#change the order to match with myallgroup#
myallabun=myallabun[rownames(myallgroup),]
rownames(myallabun)==rownames(myallgroup)

#2 transformations: relative#
myallabun_rel=decostand(myallabun,method="total")*100
apply(myallabun_rel,1,sum)
#robust central log ratio transformation# the total is fixed, thus when one bacteria goes down, the other will change accordingly. but some bacteria they all increase. the ratio stays the same no matter we use relative abundance or the absolute number (which is impossible to get)
#(no NA value)the above transformation created some NA value which should be avoided#
#(no NA value)check which column or row is NA#
myallabun_rclr=decostand(myallabun,method="rclr")
#based on the rclr distance = generate a new distance#
myallabun_rclr_euclidean=as.matrix(vegdist(myallabun_rclr,method="euclidean"))

#alpha diversity. #change this so that x axis is abundance and y axis is the group of guild (abundance across guilds across groups to see which guild changed between the groups)
myallalpha=myallgroup[,c("shannon","observed.ASVs","faith_pd","pielou_e")]
myallalpha_name=c("Shannon Index","ASV number", "Faith's phylogenetic diversity","Evenness")
myresult=data.frame(matrix(data=NA, nrow=4, ncol=3+4+4))
mytime=c("T1","T2")
mygenotype=unique(myallgroup$Genotype)

for(i in 1:4)
{
  for(j in 1:2){
    ind_tmp=which(myallgroup$Time_Point==mytime[j])
    mygroup=myallgroup[ind_tmp,]
    myalpha=myallalpha[ind_tmp,]
    tmpdata=data.frame(group=mygroup$Genotype, value=myalpha[,i])
    mytest=kruskal.test(x=tmpdata$value,g=tmpdata$group)
    myduntest=dunnTest(x=tmpdata$value,g=tmpdata$group)
    myresult[i,(4*(j-1)+1):(4*j)]=c(mytest$p.value,myduntest$res$P.unadj)
  }
  for(k in 1:length(mygenotype)){
    ind_tmp=which(myallgroup$Genotype==mygenotype[k])
    mygroup=myallgroup[ind_tmp,]
    myalpha=myallalpha[ind_tmp,]
    subject_count=table(mygroup$MouseID)
    subject_2=names(subject_count)[which(subject_count==2)]
    mygroup=mygroup[mygroup$MouseID%in%subject_2,]
    myalpha=myalpha[rownames(mygroup),]
    tmpdata=data.frame(group=mygroup$Time_Point,mouseid=mygroup$MouseID,value=myalpha[,i])
    tmpdata=tmpdata[order(tmpdata$mouseid),]
    mytest=wilcox.test(x=tmpdata$value[tmpdata$group=="T1"],y=tmpdata$value[tmpdata$group=="T2"], paired=T)
    myresult[i,9+k-1]=mytest$p.value #bc we have 3 genotypes
  }
}
rownames(myresult)=colnames(myallalpha)
colnames(myresult)=c("T1_Kruskal","T1_Intvsliver","T1_IntvsWildtype","T1_livervsWildtype",
                     "T2_Kruskal","T2_intvsliver","T2_IntvsWildtype","T2_livervsWildtype",
                     "Wildtype_T1vsT2","Int_T1vsT2","liver_T1vsT2") #same order from the Dunn-test
write.table(myresult,"t.txt",sep="\t", quote=F)


#Plots
myallalpha_name=c("Shannon Index","ASV number", "Faith's phylogenetic diversity","Evenness")
mypic=list()
for(i in 1:length(myallalpha_name)){
  drawdata=data.frame(group=myallgroup$group,value=myallalpha[,i], time=myallgroup$Time_Point,genotype=myallgroup$Genotype)
  mypic[[i]]=ggplot(data=drawdata,aes(x=time,y=value,color=genotype))+
    geom_boxplot(aes(color=group), stat = "boxplot",position = "dodge",outlier.size=2,outlier.shape=NA)+
    geom_jitter(aes(color=group),size=0.5)+
    labs(y=myallalpha_name[i],x="")+
    theme_bw()+
    theme(text=element_text(size=18))+scale_color_manual(values=myallgroup_genotypetp_color)
}
wrap_plots(mypic)+plot_annotation(tag_levels = "a")+plot_layout(guides="collect")

#for(i in 1:length(myallalpha_name)){
#  drawdata=data.frame(group=myallgroup$group,value=myallalpha[,i], time=myallgroup$Time_Point,genotype=myallgroup$Genotype)
#  mypic[[i]]=ggplot(data=drawdata,aes(x=time,y=value,fill=genotype))+
#    geom_boxplot(aes(fill=group), stat = "boxplot",position = "dodge",outlier.size=2,outlier.shape=NA)+
#    geom_jitter(size=1)+
#    labs(y=myallalpha_name[i],x="")+
#    theme_bw()+
#    theme(text=element_text(size=18))+scale_fill_manual(values=myallgroup_genotypetp_color)
#}
#wrap_plots(mypic)+plot_annotation(tag_levels = "a")+plot_layout(guides="collect")



#for(i in 1:length(myallalpha_name)){
#  drawdata=data.frame(group=myallgroup$group,value=myallalpha[,i])
#  mypic[[i]]=ggplot(data=drawdata,aes(x=group,y=value))+
#    geom_boxplot(aes(fill=group), stat = "boxplot",position = "identity",outlier.size=2,outlier.shape=NA)+
#    geom_jitter(size=1)+
#    labs(y=myallalpha_name[i],x="")+
#    coord_flip()+
#    theme_bw()+
#    theme(text=element_text(size=18), legend.position="none")+scale_fill_manual(values=myallgroup_genotypetp_color)
#}
#wrap_plots(mypic)+plot_annotation(tag_levels = "a")+plot_layout(guides="collect")



#beta diversity
#put all four into one variable, the list class, the last one is generated
BetaDistance=list(myallbraycurtis,myalljaccard,myallweiunifrac,myallunweiunifrac,myallabun_rclr_euclidean)
BetaDistance_name=c("Bray-Curtis Distance","Jaccard Distance","Weighted UniFrac Distance","Unweighted UniFrac Distance","Euclidean Distance on rclr")

#loop
ind=1:nrow(myallgroup)
mygroup=myallgroup[ind,]
#to store the result
pic=list()
pic_permanova=list()
mypermanova_result=data.frame()
#start the loop
for(i in 1:length(BetaDistance_name))
{
  tmppcoa=pcoa(BetaDistance[[i]])
  a=diag(t(tmppcoa$vectors)%*%tmppcoa$vectors)
  explain=a/sum(a)*100
  explain=round(explain,2)
  #combine
  drawdata=data.frame(mygroup,tmppcoa$vectors)
  
  # pic[[i]]=ggplot(drawdata,aes(x=Axis.1,y=Axis.2))+
  #   geom_point(aes(color=group, shape=Sex),size=5)+
  #   labs(x=paste("PC1(",explain[1],"%)"),y=paste("PC2(",explain[2],"%)"), title=BetaDistance_name[i])+
  #   stat_ellipse(aes(color=group))+
  #   labs(x=paste("PC1 (", explain[1],"%)"),y=paste("PC2 (", explain[2], "%)"), color ="", title=BetaDistance_name[i])+
  #   theme_bw()+
  #   theme(text=element_text(size=18),plot.title=element_text(hjust=0.5,size=18),legend.title=element_blank())+
  #   scale_color_manual(values=myallgroup_genotypetp_color)
  
  drawdata$newgroup=paste(drawdata$group, drawdata$Sex, sep="__")
  mydot=tmppcoa$vectors
  mydot=as.data.frame(mydot)
  mydot$group=as.factor(drawdata$newgroup)
  mydot_mean=aggregate(.~group, data=mydot,FUN=mean)
  mydot_sem=aggregate(.~group,data=mydot, FUN=function(x){sd(x)/sqrt(length(x))})
  newdrawdata=data.frame(mydot_mean, mydot_sem)
  b=data.frame(matrix(unlist(strsplit(as.character(newdrawdata$group), "__")), ncol=2,byrow=T))
  colnames(b) =c("Group", "Sex")
  newdrawdata=cbind(newdrawdata,b)
  dim(newdrawdata)
  
  pic[[i]]=ggplot(newdrawdata,aes(x=Axis.1,y=Axis.2))+
    geom_point(aes(color=Group), size=6)+
    geom_errorbarh(aes(xmin=Axis.1-Axis.1.1,xmax=Axis.1+Axis.1.1))+
    geom_errorbar(aes(ymin=Axis.2-Axis.2.1,ymax=Axis.2+Axis.2.1))+
    labs(x=paste("PC1 (", explain[1],"%)"),y=paste("PC2 (", explain[2], "%)"), color ="", title=BetaDistance_name[i])+
    theme_bw()+
    theme(text=element_text(size=18), plot.title=element_text(hjust = 0.5, size = 18), legend.title = element_blank()) + 
    scale_color_manual(values=myallgroup_genotypetp_color)
  
  
  tmppermanova=pairwise.adonis(x=BetaDistance[[i]],factors=paste(mygroup$group, mygroup$Sex, sep = "_"),p.adjust.m="BH")
  tmppermanova$Distance=BetaDistance_name[i]
  mypermanova_result=rbind(mypermanova_result, tmppermanova)
}
wrap_plots(pic)+plot_annotation(tag_levels = "a")+plot_layout(guides = "collect")
write.table(mypermanova_result, "guild_b.txt", sep = "\t", quote = F)
  #figure
  #remove everything after vs using gsub
 # tmppermanova$P1=gsub(" vs.*","",tmppermanova$pairs)
  #remove everything before vs using gsub
 # tmppermanova$P2=gsub(".*vs ","",tmppermanova$pairs)
  #change the order to make the smaller value be in the P1
 # for(ii in 1:nrow(tmppermanova))
 # {
 #   if(tmppermanova$P1[ii]<tmppermanova$P2[ii]){
 #     a=tmppermanova$P1[ii]
 #     tmppermanova$P1[ii]=tmppermanova$P2[ii]
 #     tmppermanova$P2[ii]=a
 #   }
 # }
  #add another column to add the symbol for P
 # tmppermanova$symbol=cut(tmppermanova$p.value,breaks=c(0,0.001,0.01,0.05,0.1),labels = c("****","***","**","*"))
  #figure
 # tmppermanova$P1=factor(tmppermanova$P1,levels=unique(mygroup$group))
 # tmppermanova$P2=factor(tmppermanova$P2,levels=unique(mygroup$group))
 # a=tmppermanova
 # a$P1=tmppermanova$P2
  #a$P2=tmppermanova$P1
 # b=rbind(tmppermanova,a)
 # pic_permanova[[i]]=ggplot(data=b,aes(P1,P2,fill=R2))+
  #  geom_tile()+geom_text(aes(P1,P2,label=symbol),size=10)+
 #   theme_classic()+
 #   scale_fill_gradient(low='#EDE387',high='#868b9f')+
 #   labs(fill=bquote(R^2),x="",y="",title=BetaDistance_name[i])+
 #   theme(text=element_text(size=4),plot.title=element_text(hjust = 0.5,size=18),
 #         legend.title=element_text(size=12),legend.text=element_text(size=12),axis.text=element_text(color="black"),
 #         axis.text.x=element_text(angle=90,size=14),axis.text.y=element_text(size=14))
 # tmppermanova$Distance=BetaDistance_name[i]
tmppermanova=pairwise.adonis(x=BetaDistance[[i]],factors=paste(mygroup$group, mygroup$Sex, sep = "_"),p.adjust.m="BH")
  mypermanova_result=rbind(mypermanova_result,tmppermanova)
}
wrap_plots(pic)+plot_annotation(tag_levels="a")+plot_layout(guides="collect")
wrap_plots(pic_permanova)
write.table(mypermanova_result,"t.txt",sep="\t",quote=F)


##Guild-base analysis##
#prevalence cutoff# Prevalence vs abundance# if a bacteria has a low prevalence, it will show a false positive correlation across the sample 
ind=1:nrow(myallabun_rel)
myabun_rel=myallabun_rel[ind,]
myabun_rel=myallabun_rel[,apply(myabun_rel,1,sum)>0]
#only choose the column with sum bigger than 0, thus no ASV will be 0
mygroup=myallgroup[ind,]
pre=apply(myabun_rel,2,function(x)(length(which(x>0))))
#myabun_rel is matrix, 2 means do it by column, the function that counts how many rows are not zero
#length of pre will be the column number of the abundance matrix, which is unique ASV
myresult=matrix(rep(0,21*3),ncol=3)
#the first column is the prevalence cutoff, the second columnw ill be the number of ASVs that meet the prevalence, the thirs the total abundance of ASVs meets the cutoff
for(i in 0:20)
{
  cutoff=i*5
  indtmp=which((pre/nrow(myabun_rel)*100)>=cutoff)
  myresult[i+1,1]=cutoff
  myresult[i+1,2]=length(indtmp)
  if(length(indtmp)>1)
    #at least 2 #
  {
    myresult[i+1,3]=mean(apply(myabun_rel[,indtmp],1,sum))
    #"1" means by row #
    #mean: the summed abundance of ASVs meets the prevalence cutoff, averaged across all the samples #
  }
  else
  {myresult[i+1,3]=mean(myabun_rel[,indtmp])}
  ## because we are using "apply" thus we need if and else
}
myresult=data.frame(myresult)
write.table(myresult, "t.txt", sep="\t",quote=F)
ggplot(myresult,aes(x=X1))+
  geom_line(aes(y=X2,color="ASVs"),size=1)+
  geom_line(aes(y=X3*20,color="Relative Abundance"),size=1)+
  geom_point(aes(y=X2,color="ASVs"),size=2)+
  geom_point(aes(y=X3*20,color="Relative Abundance"),size=2)+
  scale_y_continuous(sec.axis = sec_axis(~./20,name = "Relative Abundance"))+
  labs(y="# of ASVs",x="Prevalence Cutoff (%)",color="Item")+
  theme(legend.position = c(0.8,0.9),text=element_text(size=30))

#30% prevalence cutoff, based on the slope of both ASVs and relative abundance, the number of samples should be over 20#
#Bc the samples are paired, instead of Fastpar, we use repeated measurement correlation
cutoff=0.3
myabun_rel=myallabun_rel[,which(pre>=cutoff*dim(myallabun_rel)[1])]
#dim is a different method compare to nrow#
dim(myabun_rel)
myabun_rclr=myallabun_rclr[,colnames(myabun_rel)]
dim(myabun_rclr) 
#R value
mygenomecros=matrix(rep(1,dim(myabun_rclr)[2]^2),nrow=dim(myabun_rclr)[2])
#P value
mygenomecros_p=matrix(rep(0,dim(myabun_rclr)[2]^2),nrow=dim(myabun_rclr)[2])

##double check if myabun=

rownames(myabun_rclr)==rownames(myallgroup)

for(i in 2:dim(myabun_rclr)[2])
{
  for(j in 1:(i-1))
    #compare the 1:3 in i matrix to 1:1 and 1:2 in j matrix
  {tmpdata=data.frame(var1=myabun_rclr[,i],var2=myabun_rclr[,j],sub=as.factor(myallgroup$MouseID))
  tmpresult=CROS(tmpdata)
  mygenomecros[i,j]=tmpresult[1]
  mygenomecros[j,i]=tmpresult[1]
  mygenomecros_p[i,j]=tmpresult[2]
  mygenomecros_p[j,i]=tmpresult[2]
  }
}

rownames(mygenomecros)=colnames(myabun_rel)
colnames(mygenomecros)=colnames(myabun_rel)
rownames(mygenomecros_p)=colnames(myabun_rel)
colnames(mygenomecros_p)=colnames(myabun_rel)

write.table(mygenomecros,"t.txt",sep="\t",quote=F)
write.table(mygenomecros_p,"t.txt",sep="\t",quote=F)

mydist=1-mygenomecros #convert CROS correlation to distance, R matrix to distance matrix
rownames(mydist)=colnames(myabun_rel)
colnames(mydist)=colnames(myabun_rel)
mycluster=hclust(as.dist(mydist), "ward.D2") #hclust, the function to create the clustering tree
b=plot(mycluster,hang=-1, cex = 0.8) #plotting

write.table(mycluster$labels[mycluster$order],"hcluster.txt",sep="\t",quote=F)

#reorder to make the order matching with the tree
mydist_order=mydist[mycluster$labels[mycluster$order],mycluster$labels[mycluster$order]]
clade1_start=228
clade1_end=230
clade2_start=231
clade2_end=236
tmpgroup=data.frame(group=c(rep("g1",clade1_end-clade1_start+1), 
                            rep("g2",clade2_end-clade2_start+1)))

#bc we are comparing the largest groups so we are comparing all the data
tmpdist=mydist_order[clade1_start:clade2_end,clade1_start:clade2_end]
set.seed(315)
adonis2(tmpdist~group,tmpgroup,permutations = 9999) #P<0.001

#load guild result
myguildcluster = data.frame(read_xlsx("hcluster.xlsx", sheet = "Guilds"))

#guild level abundance 
myabun_rel = myallabun_rel[,myguildcluster$Item] #71 samples x 236 ASVs in guilds
tmp = t(myabun_rel) #236 ASVs x 71 samples
tmp = aggregate(tmp, by = list(myguildcluster$Result), sum) #[can use rowsum() as well]
rownames(tmp) = tmp$Group.1
tmp$Group.1 = NULL
print(dim(tmp))
tmp = data.frame(t(tmp))
print(dim(tmp))
colnames(tmp)
tmp = tmp[,order(as.numeric(gsub("CAG", "", colnames(tmp))))]
colnames(tmp)
summary(apply(myabun_rel, 1, sum))
summary(apply(tmp, 1, sum))
myallguild_rel = tmp

#### compare ASV level and guild level, betadiversity results #compare ASV and guild matrix to decide how similar the 2 matrices are 
myasv_bray=as.matrix(vegdist(myallabun_rel, method = "bray"))
myguild_bray=as.matrix(vegdist(myallguild_rel,method ="bray"))
set.seed(315)
mantel(myasv_bray, myguild_bray)
myasv_bray_pcoa = pcoa(myasv_bray) # using the PCOA plots of PC1 and PC2 from ASV and Guild, compare the clusters and see if there's a similarity
myguild_bray_pcoa = pcoa(myguild_bray)
set.seed(315)
protest(X=myasv_bray_pcoa$vectors[,1:2], Y = myguild_bray_pcoa$vectors[,1:2], choices = c(1,2)) #Correlation : the higher, the better (0.9901) in this case, #also provides sum of squares(m^2) 

#### betadiversity with sex information 

#beta diversity
#put all four into one variable, the list class, the last one is generated
BetaDistance=list(myallbraycurtis,myalljaccard,myallweiunifrac,myallunweiunifrac,myallabun_rclr_euclidean)
BetaDistance_name=c("Bray-Curtis Distance","Jaccard Distance","Weighted UniFrac Distance","Unweighted UniFrac Distance","Euclidean Distance on rclr")

#loop
ind=1:nrow(myallgroup)
mygroup=myallgroup[ind,]
#to store the result
pic=list()
pic_permanova=list()
mypermanova_result=data.frame()
#start the loop
for(i in 1:length(BetaDistance_name))
{
  tmppcoa=pcoa(BetaDistance[[i]])
  a=diag(t(tmppcoa$vectors)%*%tmppcoa$vectors)
  explain=a/sum(a)*100
  explain=round(explain,2)
  #combine
  drawdata=data.frame(mygroup,tmppcoa$vectors)
  
  pic[[i]]=ggplot(drawdata,aes(x=Axis.1,y=Axis.2))+
    geom_point(aes(color=group, shape=Sex),size=5)+
    labs(x=paste("PC1(",explain[1],"%)"),y=paste("PC2(",explain[2],"%)"), title=BetaDistance_name[i])+
    stat_ellipse(aes(color=group))+
    
    theme_bw()+
    theme(text=element_text(size=18),plot.title=element_text(hjust=0.5,size=18),legend.title=element_blank())+
    scale_color_manual(values=myallgroup_genotypetp_color)
  
  tmppermanova=pairwise.adonis(x=BetaDistance[[i]],factors=paste(mygroup$group, my group$Sex, sep = "_")p.adjust.m="BH")
  #figure
  #remove everything after vs using gsub
  tmppermanova$P1=gsub(" vs.*","",tmppermanova$pairs)
  #remove everything before vs using gsub
  tmppermanova$P2=gsub(".*vs ","",tmppermanova$pairs)
  #change the order to make the smaller value be in the P1
  for(ii in 1:nrow(tmppermanova))
  {
    if(tmppermanova$P1[ii]<tmppermanova$P2[ii]){
      a=tmppermanova$P1[ii]
      tmppermanova$P1[ii]=tmppermanova$P2[ii]
      tmppermanova$P2[ii]=a
    }
  }
  #add another column to add the symbol for P
  tmppermanova$symbol=cut(tmppermanova$p.value,breaks=c(0,0.001,0.01,0.05,0.1),labels = c("****","***","**","*"))
  #figure
  tmppermanova$P1=factor(tmppermanova$P1,levels=unique(mygroup$group))
  tmppermanova$P2=factor(tmppermanova$P2,levels=unique(mygroup$group))
  a=tmppermanova
  a$P1=tmppermanova$P2
  a$P2=tmppermanova$P1
  b=rbind(tmppermanova,a)
  pic_permanova[[i]]=ggplot(data=b,aes(P1,P2,fill=R2))+
    geom_tile()+geom_text(aes(P1,P2,label=symbol),size=10)+
    theme_classic()+
    scale_fill_gradient(low='#EDE387',high='#868b9f')+
    labs(fill=bquote(R^2),x="",y="",title=BetaDistance_name[i])+
    theme(text=element_text(size=4),plot.title=element_text(hjust = 0.5,size=18),
          legend.title=element_text(size=12),legend.text=element_text(size=12),axis.text=element_text(color="black"),
          axis.text.x=element_text(angle=90,size=14),axis.text.y=element_text(size=14))
  tmppermanova$Distance=BetaDistance_name[i]
  mypermanova_result=rbind(mypermanova_result,tmppermanova)
}
wrap_plots(pic)+plot_annotation(tag_levels="a")+plot_layout(guides="collect")
wrap_plots(pic_permanova)
write.table(mypermanova_result,"t.txt",sep="\t",quote=F)



#enter save.image('20250116_StorchLab.Rdata')
#Each ASV is one unique bacteria 
#See if bray curtis of reduced data matches with the ASV one
#do PCOA plot based on guild level, you already have PCOA plot based on ASV plot
#Next steps : Visualizing what we have done, etc, and   #focus only only relevant microbiomes and remove low prevalence bacteria

# PCOA plot based on guild level
tmppcoa=pcoa(myguild_bray)
a=diag(t(tmppcoa$vectors)%*%tmppcoa$vectors)
explain=a/sum(a)*100
explain=round(explain,2)
drawdata=data.frame(mygroup,tmppcoa$vectors)

guild_pcoa_plot=ggplot(drawdata,aes(x=Axis.1,y=Axis.2))+
  geom_point(aes(color=group, shape=Time_Point),size=5)+
  labs(x=paste("PC1(",explain[1],"%)"),y=paste("PC2(",explain[2],"%)"), title="Guild Level Bray-Curtis Distance")+
  stat_ellipse(aes(color=group))+
  theme_bw()+
  theme(text=element_text(size=18),plot.title=element_text(hjust=0.5,size=18),legend.title=element_blank())+
  scale_color_manual(values=myallgroup_genotypetp_color)

print(guild_pcoa_plot)



# Ensure myallguild_rel is a data frame
myallguild_rel <- as.data.frame(myallguild_rel)

# Add metadata (group, Time_Point, Genotype) to myallguild_rel
myallguild_rel$group <- myallgroup$group
myallguild_rel$Time_Point <- myallgroup$Time_Point
myallguild_rel$Genotype <- myallgroup$Genotype

# Melt the data for plotting
melted_guild <- melt(myallguild_rel, id.vars = c("group", "Time_Point", "Genotype"), variable.name = "Guild", value.name = "Abundance")

# Plot abundance of each guild across groups
guild_plots <- list()
for (guild in unique(melted_guild$Guild)) {
  guild_data <- subset(melted_guild, Guild == guild)
  
  guild_plot <- ggplot(guild_data, aes(x = group, y = Abundance, fill = Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    labs(title = guild, x = "Group", y = "Relative Abundance") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = myallgroup_genotypetp_color)
  
  guild_plots[[guild]] <- guild_plot
}

# Combine plots using patchwork
combined_plots <- wrap_plots(guild_plots, ncol = 3) + plot_annotation(tag_levels = "A")
print(combined_plots)

##Guild level abundance of RCLR range mean; rclr euclidean## #Compare rclr ASVs level vs Guild level#### 3/13/25
#extract data for the rclr
tmpabunrclr=myallabun_rclr[,myguildcluster$Item]
tmpabunrclr_range=decostand(tmpabunrclr,method="range") #range scale# #transform the value to 0-1 and then put them in a relative order
tmpabun=aggregate(t(tmpabunrclr_range),by=list(myguildcluster$Result),mean) #(transform the data, get the guild information, average the abundance of all ASVs from the same guild in each sample)#
rownames(tmpabun)=tmpabun$Group.1
tmpabun$Group.1=NULL
tmpabun=data.frame(t(tmpabun))
myallguild_rclr_range_mean=tmpabun #the second guild abundance system vs the sum#

myasv_rclr_euclidean=as.matrix(vegdist(myallabun_rclr,method="euclidean")) #the vegdist output is not a matrix
myguild_rclr_euclidean=as.matrix(vegdist(myallguild_rclr_range_mean,method="euclidean"))

set.seed(315)
mantel(myasv_rclr_euclidean,myguild_rclr_euclidean)
myasv_rclr_euclidean_pcoa=pcoa(myasv_rclr_euclidean)
myguild_rclr_euclidean_pcoa=pcoa(myguild_rclr_euclidean)

set.seed(315)
protest(X=myasv_rclr_euclidean_pcoa$vectors[,1:2],Y=myguild_rclr_euclidean_pcoa$vectors[,1:2],choices=c(1,2))

#dev.off()
#Generating the Heatmap #8x20 to export full size heatmap
#Create 2 heatmaps separate my timepoints
rownames(myallgroup)==rownames(myallguild_rel)
mygroup=myallgroup[myallgroup$Time_Point=="T1",]
mygroup=mygroup[order(mygroup$Genotype,mygroup$Time_Point,mygroup$Sex),]
myguild_rel=myallguild_rel[rownames(mygroup),]
tmpannotation_col=data.frame(Sex=mygroup$Sex, Genotype=mygroup$Genotype)
rownames(tmpannotation_col)=rownames(myguild_rel)
mysex_color=c("F"="#000000", "M"="#FFFFFF")
tmpannotation_color=list(Sex=mysex_color,Time=myallgroup_time_color,Genotype=myallgroup_genotype_color)
pheatmap(t(log10(myguild_rel+1)),cluster_rows=F, cluster_cols=F, show_colnames=F,
         color=colorRampPalette(c("#0D1740","white","#EDE387","#F17C67"))(50),cellheight = 12,cellwidth = 12,
         #gaps_col = c(6,12,18,24,30),fontsize = 10, #T1
         gaps_col = c(6,12,18,23,29),fontsize = 10,  #T2
         annotation_col = tmpannotation_col,annotation_colors = tmpannotation_color)

mygroup=myallgroup[myallgroup$Time_Point=="T1",]. #Replace T1 and T2
mygroup=mygroup[order(mygroup$Genotype,mygroup$Time_Point,mygroup$Sex),]
myguild_rel=myallguild_rel[rownames(mygroup),]
myresult=data.frame(matrix(data=NA,nrow=ncol(myguild_rel),ncol=12))
mysex=unique(mygroup$Sex)
for(i in 1:ncol(myguild_rel))
{
  for(j in 1:2){
    ind_tmp=which(myallgroup$Sex==mysex[j])
    tmpmygroup=mygroup[ind_tmp,]
    myalpha=myguild_rel[ind_tmp,]
    tmpdata=data.frame(group=tmpmygroup$Genotype,value=myalpha[,i])
    mytest=kruskal.test(x=tmpdata$value,g=tmpdata$group)
    myduntest=dunnTest(x=tmpdata$value,g=tmpdata$group,method = "bh")
    tmpadjp=myduntest$res[,4]
    names(tmpadjp)=gsub(" ","",myduntest$res$Comparison)
    tmpletter=multcompLetters(tmpadjp)$Letters
    myresult[i,(3*(j-1)+1):(3*j)]=tmpletter
  }
  for(k in 1:length(mygenotype)) {
    ind_tmp=which(mygroup$Genotype==mygenotype[k])
    tmpmygroup=mygroup[ind_tmp,]
    myalpha=myguild_rel[ind_tmp,]
    tmpdata=data.frame(group=tmpmygroup$Sex,value=myalpha[,i])
    mytest=wilcox.test(value~group,tmpdata)
    myresult[i,6+k]=mytest$p.value
    
  }
  myresult[i,10:12]=p.adjust(myresult[i,7:9],method="BH")
}
rownames(myresult)=colnames(myallguild_rel)
colnames(myresult)=c(paste(mysex[1],names(tmpletter),sep = "_"),
                     paste(mysex[2],names(tmpletter),sep = "_"),
                     paste("RawP",mygenotype,sep = "_"),
                     paste("BHadjp",mygenotype,sep = "_"))
write.table(myresult,"t.txt",sep="\t", quote=F)

#Save results in excel and separate which ones at T1 and which ones are T2


#each column in one gender and on each gender is the genotypes (LEFT) 
