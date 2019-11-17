## Bolivia dung beetle study - data from 2011 and 2015

## Laura Allen

## 17 November 2019

rm(list=ls()) ##Clear memory

## packages----
# for ordinations
library(vegan)
library(ape)

# + for beta diversity, LCBD etc and plotting
library(adespatial)
require(vegan)  
require(adegraphics)
library(ade4)
require(RColorBrewer)

# alpha divversity hill number
library(iNEXT)
library(ggplot2)


#Ordination

########################
# Import data ----
# species x site
dung  <- read.csv("C:/Data/Bolivia-dungbeetles/Rawdata/species_sites_2019.csv",row.names=1) #column 1 has row names
summary(dung) # check smmary makes sense
tdung <- t(dung) #transpode data
rownames(tdung) # check site names are recognised

# community composition ----
# PCA to compare sites 

#custom pca plot
dung.transf = decostand(tdung,"hel")
rda.out <- rda(dung.transf) 
# pc1=0.24119 pc2=0.08143
scrs <- scores(rda.out,display = c("sites", "species"), scaling =1) #extract sepecies and site scores from rda summary
summary(rda.out)

biplot(rda.out, scaling=1,display=c("sites","species"))

## show the most stongly associated species in each direction (not all species, to reduce crowding)
spc1 <- scrs$species[order(scrs$species[,1]),]
spc2 <- scrs$species[order(scrs$species[,2]),]
sppc <- rbind(spc1[c(1,2,33,34),],spc2[c(1,2,33,34),])
#bw friendly version
tiff(file="C:/Data/Bolivia-dungbeetles/Outputs/PCA_commcomp_bw.tiff",width=190,height=190,units="mm",res=1000,pointsize=12)  
xr <- range(scrs$species[,1], scrs$sites[,1])
yr <- range(scrs$species[,2], scrs$sites[,2])
xr[1] <- xr[1]-0.2; xr[2] <- xr[2]+0.2; yr[1] <- yr[1]-0.2; yr[2] <- yr[2]+0.2
par(mar=c(1,1,1,1))
plot.new()
plot.window(xlim = xr, ylim = yr, asp = 1)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(scrs$sites, col="black", bg= "grey",pch=c(0,15,4,6,2,17,24,1,5),
                 cex=3,lwd=2)
text(sppc,labels = row.names(sppc),col="red")
legend("bottomleft",bty="n",pch=c(0,15,4,6,2,17,24,1,5),pt.cex=1.2,cex=1,col="black", pt.bg="grey",legend=c(row.names(scrs$sites)))
text(-0.03,-1.1,"PC2",srt=90,font=2)
text(-1.1,-0.03,"PC1",font=2)
text(c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),-0.01,c("-1","-0.8","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.8","1"))
text(-0.02,c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1),c("-1","-0.8","-0.6","-0.4","-0.2","0.2","0.4","0.6","0.8","1"))
dev.off()

## hills alpha diversity
head(dung)
db.mat <- as.matrix(dung)
max(colSums(db.mat))
iN0 <- iNEXT(db.mat, q=0, datatype="abundance", size=NULL, endpoint = 1000, se=TRUE, conf=0.95, nboot=50)

iN123 <- iNEXT(db.mat, q=c(0,0.5,1,2,Inf), datatype="abundance", size=NULL, endpoint = 1000, se=TRUE, conf=0.95, nboot=50)
ggiNEXT(iN0, type=1, se=TRUE, grey=FALSE)  
