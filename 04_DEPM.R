library("edgeR")
# create DGE list
# count is refers to the normalized expression matrix
# submap refers to the metadata file of samples

d = DGEList(counts=count, group=sub_map$Subtype)
d = calcNormFactors(d)

map.mat = model.matrix(~ 0 + d$samples$Subtype)
colnames(map.mat)=levels(sub_map$Subtype)
d2 = estimateGLMCommonDisp(d, map.mat)
d2 = estimateGLMTagwiseDisp(d2, map.mat)
fit = glmFit(d2, map.mat)

 BvsA <- makeContrasts(contrasts = "LBC-TNBC", levels=map.mat)

 lrt = glmLRT(fit,contrast=BvsA)
 # FDR test, less than 0.5%
 de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
 x=lrt$table
 x$sig=de_lrt
 J_H_enriched = row.names(subset(x,sig==1))
 J_H_depleted = row.names(subset(x,sig==-1))
 
 #热图展示差异OTU
 pair_group = subset(sub_map2,Subtype %in% c("LBC","TNBC"))
 DE = c(J_H_enriched,J_H_depleted)
 sub_norm = as.matrix(norm[DE,rownames(pair_group)])
 pdf(file= paste("heat_otu_BvsL_sig.pdf",sep = ""),height= 8,width= 8)
 heatmap.2(sub_norm,scale = "row", Rowv=FALSE, Colv=FALSE, dendrogram='row', trace='none', xlab =  col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),cexCol = 0.65,keysize=0.8,cexRow=0.8,density.info = "none",main = NULL)
 dev.off()
pdf(file= paste("heat_otu_BvsL_sig.pdf",sep = ""),height= 8,width= 12)
heatmap.2(sub_norm,col=cm.colors(255),scale = "col", Rowv=FALSE, Colv=FALSE, dendrogram='none', trace='none', xlab = sub_map$Subtype,main = "LBC-TNBC DEPM") 
dev.off()
