#!/bin/bash 

entry=$1
# cutoff=$2 for those cases that there is a mismatch between the atoms in the original protein and the trajectory

if [ $# -eq 0 ]; then
        echo
        echo -e "Usage: MD_analysis_HMDB.sh \t entry_name"
	echo -e "there will be few cases where a cutoff to fix residue number will be needed"
        echo
        echo -e "\t This script writes out the R script to run the analysis on trajectories"
        echo -e "\t & generates all possible pdf files and some pdb files" 
        echo
        exit 0
fi

#echo -e "/usr/local/bin/Rscript --vanilla MD.R"

cat > MD.R << EOF 
library(bio3d)
library(cluster)
library(grid)
library(lattice)

plot.pca2 <-
function(x, pch=16, col=par("col"), cex=0.8, mar=c(4, 4, 1, 1), ...) {

  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  par(mfrow=c(2, 2), cex=cex, mar=mar)
  par(pty="s")
  plot(x\$z[,1],x\$z[,2], type="p", pch=pch, xlab="PC1", ylab="PC2", col=col, ...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x\$z[,3],x\$z[,2], type="p", pch=pch, xlab="PC3", ylab="PC2", col=col,...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x\$z[,1],x\$z[,3], type="p", pch=pch, xlab="PC1", ylab="PC3", col=col,... )
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot.pca.scree(x\$L,ylim=c(0,30))
}


##-- Read a trajectory file
trj <- read.ncdf("$entry.12ns_trj.nc")
#trj <- trj[1:1000,]  ### FOR TESTING!

##-- Read a PDB file with the CA from the minimized structure
##CA<-read.pdb("ErbB3.gas.pdb")
pdb <- read.pdb("$entry.gas.pdb")
# Trajectory Frame Superposition
ca.inds <- atom.select(pdb, elety="CA")
len<- length(ca.inds\$atom)
ca_gapped <- seq(1,len,10)
inds<-atom.select(pdb, "calpha")\$atom

residues<-aa321(pdb\$atom[inds,'resid'])
resnumber<-as.numeric(pdb\$atom[inds,'resno'])
names <- paste(residues,resnumber)

xyz <- fit.xyz(fixed=pdb\$xyz, mobile=trj, 
	fixed.inds=ca.inds\$xyz, 
	mobile.inds=ca.inds\$xyz)

# Root Mean Square Deviation (RMSD)
rd <- rmsd(xyz[1,ca.inds\$xyz], xyz[,ca.inds\$xyz])

pdf("rmsd.pdf")
par(mfrow=c(3,1))#,cex=1.2)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD") 
lines(density(rd), col="gray", lwd=3)
time <- seq(1,dim(trj)[1])*0.5
plot.bio3d(rd, time, typ="l", ylab="RMSD from first frame",xlab="time (ps)")#,ylim=c(0,2.5))
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
#points(lowess(rd[2001:6839]), typ="l", col="blue", lty=2, lwd=2)
dev.off()

summary(rd)

# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[,ca.inds\$xyz])
fluc <- rmsf(trj)
write.pdb(pdb, b=fluc, file="fluct.pdb")

pdf("rmsf.pdf", height=4, width=12, onefile=TRUE)
par(mfrow=c(1,1))
par(mar=c(5,4,4,4) + 0.1, cex=1.0)
#plot.bio3d(fluc, sse=sse, typ="l", ylim=c(0,8), ylab="RMSF (A)", xlab="Residues",helix.col = "red", sheet.col = "yellow", sse.border = FALSE, axes=FALSE)
plot.bio3d(rf, typ="l", ylim=c(0,10), ylab="RMSF (A)", xlab="Residues", axes=FALSE)
axis(2)
axis(1,at=ca_gapped,labels=names[ca_gapped],las=2)
#axis(3,at=ca_gapped,labels=names[ca_gapped],las=2)

#par(new=TRUE)
# if there is a pdb file
#plot(pdb\$atom[pdb\$calpha,"b"],type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n')
####plot(pdb\$atom[pdb\$calpha,"b"],type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n')
#axis(4, labels=FALSE)
#at = axTicks(4)
#mtext(side = 4, text = at, at = at, col = "blue", line = 1) 
#mtext('B-Factor',side=4,line=3,col='blue')
box()

dev.off()
#plot(rf, ylab="RMSF", xlab="Residue Position",typ="l")

# Principal Component Analysis
trj.pca <- pca.xyz(xyz[,ca.inds\$xyz])
pdf("pca.pdf")
plot.pca2(trj.pca,col=rainbow(length(trj.pca\$z[,1])))
dev.off()

#plot(pc, col=bwr.colors(nrow(xyz)) )

## Plot atom wise loadings
# it is included within pca.loadings
pdf("pca_rmsf.pdf")
x <- c(1:len)
par(mfrow=c(3,1),mar=c(3,4,1,1))
plot.bio3d(trj.pca\$au[, 1],x, typ="l", ylim=c(0,0.15),ylab="PC1 (A)", xlab="Residue Position")
plot.bio3d(trj.pca\$au[, 2],x, typ="l", ylim=c(0,0.15),ylab="PC2 (A)", xlab="Residue Position")
plot.bio3d(trj.pca\$au[, 3],x, typ="l", ylim=c(0,0.15),ylab="PC3 (A)", xlab="Residue Position")
dev.off()

#plot.bio3d(pc\$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l") 
#points(pc\$au[,2], typ="l", col="blue")

# Conformer Clustering in PC space. Clustering of structures from the trajectory in the PC1-3 planes
pdf("pca_cluster.pdf")
hc <- hclust(dist(trj.pca\$z[, 1:3]))
grps <- cutree(hc, k=3)
#png("pc_clustering.png",res=300, pointsize=2, width=1200, height=1200)
plot.pca2(trj.pca, col = c("red", "black", "blue")[grps], cex=1.0)
###plot(hclust(dist(trj.pca\$z[, 1:3])))
dev.off()

# Write PC trajectory
p1 <- mktrj.pca(trj.pca, pc=1,b=trj.pca\$au[,1], file="pc1.pdb",resno = resnumber, resid = residues) 
p2 <- mktrj.pca(trj.pca, pc=2,b=trj.pca\$au[,2], file="pc2.pdb",resno = resnumber, resid = residues)
p3 <- mktrj.pca(trj.pca, pc=3,b=trj.pca\$au[,3], file="pc3.pdb",resno = resnumber, resid = residues)

write.ncdf(p1, "trj_pc1.nc")

# Cross-Correlation Analysis 
cij <- dccm(xyz[,ca.inds\$xyz])
pdf("dccm.pdf", height=12, width=12, onefile=TRUE)

x <- y <- c(1:len)
ca_gapped <- seq(1,len,5)
inds<-atom.select(pdb, "calpha")\$atom

levels=c(-1, -0.7, -0.4, 0.4, 0.7, 1)
#colors=c('darkblue','skyblue','white','pink','red2')
colors=c('red2','pink','white','skyblue','darkblue')

# Color Scheme 1
filled.contour(x, y, cij,col=colors,levels=levels,xlab="Residue Number", ylab="Residue Number",main="DCCM: dynamic cross-correlation map",zlim=c(-1,1),plot.axes={axis(1,at=ca_gapped,labels=names[ca_gapped],las=2); axis(2,at=ca_gapped,labels=names[ca_gapped],las=2);box()})
dev.off()

pdf("dccm2.pdf")
myColour = colorRampPalette(c("red", "white","blue"))
colkey <- colorRampPalette(c('darkred','white','darkblue'))(50) #32)

# This part has to be optimized to include those regions of interest
##### Lim you need to find your regions

regions <-c(125,250,413,564)
regions.name <-c("Domain I","Domain II","Domain III ","Domain IV")

scales <- list(at=regions,
                  labels=regions.name,alternating=2,ticks=F)

ca_gapped <- seq(1,len,10)

plot.dccm(cij, 
        row.values = seq_len(nrow(cij)),
        column.values = seq_len(ncol(cij)), 
        #sse = sse, helix.col = "gray20", sheet.col = "gray80", inner.box = TRUE, outer.box = FALSE,
        scales=scales,
        at=c(-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1), 
        col.regions=(col=colkey),
        axes=FALSE, 
                #plot.axes={axis(1,at=ca_gapped,labels=names[ca_gapped],las=2); axis(2,at=ca_gapped,labels=names[ca_gapped],las=2);  box()},
    xlab="Residue No. $entry",
    ylab="Residue No. $entry",
    main="DCCM: dynamic cross-correlation map")

dev.off()
# View the correlations in pymol
#view.dccm(cij, pdb, launch=TRUE)
view.dccm(cij, pdb, launch=F)
EOF 


/usr/local/bin/Rscript --vanilla MD.R
echo -e "/usr/local/bin/Rscript --vanilla MD.R"
r --vanilla < MD.R
#rm -f MD.R
