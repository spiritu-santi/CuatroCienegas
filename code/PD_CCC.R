library(tidyverse)
library(picante)
setwd("~/Documents/1.PROYECTOS/13.PD_CCC/") #### Cambia tu directorio aqui
### There is a stupid transposing of tables....... Correct!!!
nullmodels.diversity <- function(input.tree,input.table,output.base,root=FALSE, rarify=TRUE) {
tree <- ape::read.tree(input.tree)
if(rooted) tree <- ape::root.phylo(tree,"ASV_1",resolve.root = T)
tabla <- read.table(input.table,sep="\t",header=T)
if(rarify){  
raremax <- min(colSums(tabla))
tabla <- t(vegan::rrarefy(t(tabla), sample=raremax))
}
print(colSums(tabla))
readline("Press enter if OK (if not abort)")
print(summary(rowSums(tabla)))
empirical <- tibble(Site=colnames(tabla),PSC_total=NA, SR_total=NA, PD_total=NA,SRp_total=NA, PSC.rare=NA,SR.rare=NA,PSC.abund=NA,SR.abund=NA,
                          PD.rare=NA,SRp.rare=NA,PD.abund=NA,SRp.abund=NA,.rows=dim(tabla)[2])
#### Estimados de Clustering por sitio para la muestra total (más cercano a cero == más clustering)
empirical[,2:3] <-  picante::psc(t(tabla),tree,scale.vcv = TRUE)
#### Estimados de Faith's PD por sitio para la muestra total
empirical[,4:5]<- picante::pd(t(tabla),tree,include.root = TRUE) 
tabla_na <- tabla
tabla_na <- apply(tabla_na,2,function(x) as.numeric(sub(0,NA,x)))
d <- apply(tabla_na,2,summary,na.rm=T)
rares.site <- ceiling(d[seq(2,70,7)])
abund.site <- floor(d[seq(5,70,7)])

for (i in 1:dim(tabla)[2]) {
  cat("Calculating empirical diversity indices for site",empirical[[i,1]],"\n")
  empirical[i,6:7] <- picante::psc(t(tabla)[i,which(t(tabla)[i,] < rares.site[i])], tree, scale.vcv = F)
  empirical[i,8:9] <- picante::psc(t(tabla)[i,which(t(tabla)[i,] > abund.site[i])], tree, scale.vcv = F)
                  uno <- t(tabla[which(tabla[,i] < rares.site[i]),i])
                  colnames(uno) <- rownames(tabla[which(tabla[,i] < rares.site[i]),])
                  empirical[i,10:11] <- picante::pd(as.matrix(uno), tree, include.root = TRUE)
                  uno <- t(tabla[which(tabla[,i] > abund.site[i]),i])
                  colnames(uno) <- rownames(tabla[which(tabla[,i] > abund.site[i]),])
                  empirical[i,12:13] <- picante::pd(as.matrix(uno), tree, include.root = TRUE)
                  }



### Standardized Effect Sizes (1000 réplicas)
##### Regular for loop.....
# system.time({ 
#   null_model <- list()
#   for (k in 1:2){ 
#     cat(toupper("Randomization"),k,"\n")
#     tabla_null <- randomizeMatrix(tabla,null.model = "frequency") ### the matrix is flipped, so this is in fact "richness"
#     indices <- tibble(Site=colnames(tabla_null),PSC_total=NA, SR_total=NA, PD_total=NA,SRp_total=NA, PSC.rare=NA,SR.rare=NA,PSC.abund=NA,SR.abund=NA,
#                       PD.rare=NA,SRp.rare=NA,PD.abund=NA,SRp.abund=NA,.rows=10)
#     indices[,2:3] <-  picante::psc(t(tabla_null),tree,scale.vcv = TRUE)
#     indices[,4:5]<- picante::pd(t(tabla_null),tree,include.root = TRUE) 
#     tabla_na <- tabla_null
#     tabla_na <- apply(tabla_na,2,function(x) as.numeric(sub(0,NA,x)))
#     d <- apply(tabla_na,2,summary,na.rm=T)
#     rares.site <- ceiling(d[seq(2,70,7)])
#     abund.site <- floor(d[seq(5,70,7)])
#     for (i in 1:10) {
#       cat("Calculating for site",indices[[i,1]],"\n")
#       indices[i,6:7] <- picante::psc(t(tabla_null)[i,which(t(tabla_null)[i,] < rares.site[i])], tree, scale.vcv = TRUE)
#       indices[i,8:9] <- picante::psc(t(tabla_null)[i,which(t(tabla_null)[i,] > abund.site[i])], tree, scale.vcv = TRUE)
#       uno <- t(tabla_null[which(tabla_null[,i] < rares.site[i]),i])
#       colnames(uno) <- rownames(tabla_null[which(tabla_null[,i] < rares.site[i]),])
#       indices[i,10:11] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
#       uno <- t(tabla_null[which(tabla_null[,i] > abund.site[i]),i])
#       colnames(uno) <- rownames(tabla_null[which(tabla_null[,i] > abund.site[i]),])
#       indices[i,12:13] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
#     }
#     null_model[[k]] <- indices
#   }
#   empirical
# })

##### Try with futures
library(furrr)
tabla_null <- list()
cat("Creating null models (richness)")
for (i in 1:100) tabla_null[[i]] <- randomizeMatrix(tabla, null.model = "frequency") ### the matrix is flipped, so this is in fact "richness": mantains site's richness.
plan(multisession)
#system.time({ 
null_model <- tabla_null %>% furrr::future_map(function(x){
  indices <- tibble(Site=colnames(x),PSC_total=NA, SR_total=NA, PD_total=NA,SRp_total=NA, PSC.rare=NA,SR.rare=NA,PSC.abund=NA,SR.abund=NA,
                    PD.rare=NA,SRp.rare=NA,PD.abund=NA,SRp.abund=NA,.rows=dim(tabla)[2])
  indices[,2:3] <-  picante::psc(t(x),tree,scale.vcv = TRUE)
  indices[,4:5]<- picante::pd(t(x),tree,include.root = rooted) 
  tabla_na <- x
  tabla_na <- apply(tabla_na,2,function(xxx) as.numeric(sub(0,NA,xxx)))
  d <- apply(tabla_na,2,summary,na.rm=T)
  rares.site <- ceiling(d[seq(2,70,7)])
  abund.site <- floor(d[seq(5,70,7)])
  for (i in 1:dim(tabla)[2]) {
    indices[i,6:7] <- picante::psc(t(x)[i,which(t(x)[i,] < rares.site[i])], tree, scale.vcv = TRUE)
    indices[i,8:9] <- picante::psc(t(x)[i,which(t(x)[i,] > abund.site[i])], tree, scale.vcv = TRUE)
    uno <- t(x[which(x[,i] < rares.site[i]),i])
    colnames(uno) <- rownames(x[which(x[,i] < rares.site[i]),])
    indices[i,10:11] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
    uno <- t(x[which(x[,i] > abund.site[i]),i])
    colnames(uno) <- rownames(x[which(x[,i] > abund.site[i]),])
    indices[i,12:13] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
  }
  return(indices)
},.progress=TRUE)
#})
arr <- array( unlist(lapply(null_model,"[",-1)) , c(dim(tabla)[2],12,length(lapply(null_model,"[",-1))))
mean_indices <- apply( arr , 1:2 , function(x) round(mean(x),5)) %>% as_tibble(.) %>% bind_cols(Site = empirical$Site,.)
colnames(mean_indices) <- colnames(empirical)
sd_indices <- apply( arr , 1:2 , function(x) round(sd(x),5)) %>% as_tibble(.) %>% bind_cols(Site=empirical$Site,.)
colnames(sd_indices) <- colnames(empirical)
SES_richness <- ( empirical[,-1] - mean_indices[,-1] ) / sd_indices[,-1]
colSums(tabla)==colSums(tabla_null[[1]])

tabla_null <- list()
cat("Creating null models (richness)")
for (i in 1:100) tabla_null[[i]] <- randomizeMatrix(tabla, null.model = "richness") ### the matrix is flipped, so this is in fact "frequency": mantains species abundances.
#system.time({ 
  null_model <- tabla_null %>% furrr::future_map(function(x){
    indices <- tibble(Site=colnames(x),PSC_total=NA, SR_total=NA, PD_total=NA,SRp_total=NA, PSC.rare=NA,SR.rare=NA,PSC.abund=NA,SR.abund=NA,
                      PD.rare=NA,SRp.rare=NA,PD.abund=NA,SRp.abund=NA,.rows=dim(tabla)[2])
    indices[,2:3] <-  picante::psc(t(x),tree,scale.vcv = TRUE)
    indices[,4:5]<- picante::pd(t(x),tree,include.root = rooted) 
    tabla_na <- x
    tabla_na <- apply(tabla_na,2,function(xxx) as.numeric(sub(0,NA,xxx)))
    d <- apply(tabla_na,2,summary,na.rm=T)
    rares.site <- ceiling(d[seq(2,70,7)])
    abund.site <- floor(d[seq(5,70,7)])
    for (i in 1:dim(tabla)[2]) {
      indices[i,6:7] <- picante::psc(t(x)[i,which(t(x)[i,] < rares.site[i])], tree, scale.vcv = TRUE)
      indices[i,8:9] <- picante::psc(t(x)[i,which(t(x)[i,] > abund.site[i])], tree, scale.vcv = TRUE)
      uno <- t(x[which(x[,i] < rares.site[i]),i])
      colnames(uno) <- rownames(x[which(x[,i] < rares.site[i]),])
      indices[i,10:11] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
      uno <- t(x[which(x[,i] > abund.site[i]),i])
      colnames(uno) <- rownames(x[which(x[,i] > abund.site[i]),])
      indices[i,12:13] <- picante::pd(as.matrix(uno), tree,include.root = TRUE)
    }
    return(indices)
  },.progress=TRUE)
#})
plan(sequential)
arr <- array( unlist(lapply(null_model,"[",-1)) , c(dim(tabla)[2],12,length(lapply(null_model,"[",-1))))
mean_indices <- apply( arr , 1:2 , function(x) round(mean(x),5)) %>% as_tibble(.) %>% bind_cols(Site=empirical$Site,.)
colnames(mean_indices) <- colnames(empirical)
sd_indices <- apply( arr , 1:2 , function(x) round(sd(x),5)) %>% as_tibble(.) %>% bind_cols(Site=empirical$Site,.)
colnames(sd_indices) <- colnames(empirical)

SES_frequency <- ( empirical[,-1] - mean_indices[,-1] ) / sd_indices[,-1]
rowSums(tabla)==rowSums(tabla_null[[1]])
SES_frequency <- SES_frequency %>% as_tibble(.) %>% bind_cols(Site=empirical$Site,.)
SES_richness <- SES_richness %>% as_tibble(.) %>% bind_cols(Site=empirical$Site,.)
print(empirical)
print(SES_richness)
print(SES_frequency)


write.table(SES_richness, paste("output/SESrich_",output.base,".csv",sep=""),sep=",",row.names = F, col.names = T,quote = F)
write.table(empirical, paste("output/Empirical_",output.base,".csv",sep=""),sep=",",row.names = F, col.names = T,quote = F)
write.table(SES_frequency, paste("output/SESfreq_",output.base,".csv",sep=""),sep=",",row.names = F, col.names = T,quote = F)
}

nullmodels.diversity(input.tree = "data/Tree_16S.tre",input.table = "data/Table_16S.txt",output.base = "16S",root = FALSE,rarify = FALSE)









plot(indices$PD_total,indices$PSC_total,pch=21,bg="gray50",ylim=c(min(indices$PSC_total-0.01),max(indices$PSC_total+0.01)),
     xlab="Faith's PD", ylab="Phylogenetic Species Clustering")
text(indices$PD_total,indices$PSC_total,labels=rownames(pd_tabla),pos=3,offset=1)
abline(a=0,b=1,lty=3,col="gray50",lwd=2)

plot(pd_tabla$PD,pd_rare$PD,pch=21,bg="gray50",xlim=c(0,1000),ylim=c(0,1000),
     xlab="Faith's PD", ylab="Faith's PD (most rare OTUs)")
text(pd_tabla$PD,pd_rare$PD,labels=rownames(pd_tabla),pos=3,offset=1)
abline(a=0,b=1,lty=3,col="gray50",lwd=2)


plot(pd_tabla$PD,pd_abund$PD,pch=21,bg="gray50",xlim=c(0,1000),ylim=c(0,1000),
     xlab="Faith's PD", ylab="Faith's PD (most abundant OTUs)")
text(pd_tabla$PD,pd_abund$PD,labels=rownames(pd_tabla),pos=3,offset=1)
abline(a=0,b=1,lty=3,col="gray50",lwd=2)

plot(pd_rare$PD,pd_abund$PD,pch=21,bg="gray50",xlim=c(0,400),ylim=c(0,400),
     xlab="Faith's PD (most rare OTUs)", ylab="Faith's PD (most abundant OTUs)")
text(pd_rare$PD,pd_abund$PD,labels=rownames(pd_rare),pos=3,offset=1)
abline(a=0,b=1,lty=3,col="gray50",lwd=2)



plot(psc_rare$PSCs,psc_abund$PSCs,pch=21,bg="gray50",xlim=c(0,1),ylim=c(0,1),
     xlab="Clustering (most rare OTUs)", ylab="Clustering (most abundant OTUs)")
text(psc_rare$PSCs,psc_abund$PSCs,labels=rownames(psc_rare),pos=3,offset=1)
abline(a=0,b=1,lty=3,col="gray50",lwd=2)

