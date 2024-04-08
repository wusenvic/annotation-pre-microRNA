#all the script comes from https://www.biowolf.cn/
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("digest")
#install.packages("GOplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GOplot")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

pvalueFilter=0.05     #p值过滤条件
qvalueFilter=1        #矫正后的p值过滤条件

#定义图形的颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("D:\\AI\\trans\\breast\\combined\\mRNA")       #设置工作目录
rt=read.csv("gene03.txt", header=T, sep=",", check.names=F)     #读取输入文件

#获取核心基因的名称, 将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(GO, file="gene03_GO.txt", sep="\t", quote=F, row.names = F)

#柱状图
pdf(file="gene03_go_barplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#气泡图
pdf(file="gene03_go_bubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#install.packages("circlize")
#install.packages("RColorBrewer")


#引用包
library(circlize)
library(RColorBrewer)

input="down-GO.txt"      #输入文件
outpdf="down-GO2.circos.pdf"     #输出的结果文件
#setwd("C:\\Users\\Dell\\Desktop\\GOcircos.high")     #设置工作目录

#读取输入文件
data=read.table(input, header=T, sep="\t", check.names=F)
data=data[order(data$pvalue),]     #根据pvalue对数据进行排序
data=head(data, n=8)      #展示富集最显著的前八个GO
pvalue=round(-log(data$pvalue,10),2)
#获取展示基因的列表
genelist=unlist(strsplit(data$geneID,"/"))
genetable=table(genelist)

#定义图形的颜色
GO_col = colorRampPalette(brewer.pal(9, "Paired"))(nrow(data))
genes_col = rainbow(length(names(genetable)),s=0.7,v=0.7)
pvalue_col_fun <- colorRamp2(
	breaks = c(min(pvalue), mean(pvalue), max(pvalue)), 
	colors = c("#FFFFCC", "#FF9966", "#FF0000")
)

#获取基因和功能的对应关系
genedict = list()
for(i in names(genetable)) genedict[[i]] = 0
GOdict = list()
for(i in data$Description) GOdict[[i]] = 0
# link  gene--GO
gene_GO_list = list()
n = 0
for(i in 1:nrow(data)){
  genei = strsplit(data$geneID[i],"/")[[1]] # multi
  GOi = data$Description[i]
  for(j in 1:length(genei)){
    n = n+1
    genej = genei[j] # single
    gene_GO_list[[n]] = data.frame(gene=genej,GO=data$Description[i],start=genedict[[genej]],
                                        end=genedict[[genej]]+1,start2 = GOdict[[GOi]],end2=GOdict[[GOi]]+1,
                                        pvalue=pvalue[i])
    genedict[[genej]] = genedict[[genej]]+1
    GOdict[[GOi]] = GOdict[[GOi]]+1
  }
}
#gene \t GO \t start \t end \t pvalue
gene_GO_data = as.data.frame(do.call('rbind',gene_GO_list))
gene_GO_data$linkcol = GO_col[as.numeric(as.factor(gene_GO_data$GO))]
#right GO
data3 = data.frame(id=data$Description,start = 0, end = data$Count)
# left top
data1 = data.frame(id=names(genetable),start=0,end = as.numeric(genetable))
# main chrom
df = as.data.frame(rbind(data3,data1))

#定义显著性的标记
get_sig = function(p){
  ifelse(p> -log(0.001,10),"***",ifelse(p> -log(0.01,10),'**','*'))
}

bed3 = data.frame(data3,yt=0,yb=1,col=GO_col[as.numeric(as.factor(data3$id))],p=0,text='')
bed1 = data.frame(data1,yt=0.5,yb=1,col=genes_col[as.numeric(as.factor(data1$id))],p=0,text='')
bed2 = data.frame(id=gene_GO_data$gene,start=gene_GO_data$start,
                  end=gene_GO_data$end,yt=0,yb=0.5,
                  col=pvalue_col_fun(gene_GO_data$pvalue),p=gene_GO_data$pvalue,
                  text=get_sig(gene_GO_data$pvalue))

bed = as.data.frame(rbind(bed1,bed2,bed3))

#绘制图形
pdf(file=outpdf, width=12, height=7)
layout(mat=matrix(c(1,1,1,0,2,3),nc=2),width=c(6.5,3.5),height=c(2,2,7))
#初始化图形
circos.par(track.margin=c(0.01,0.01), start.degree=90)
circos.genomicInitialize(df,plotType="none")

#展示基因的名称
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  if(!any(data3$id%in%sector.index)){
  circos.text(mean(xlim), mean(ylim), sector.index, facing = "reverse.clockwise",cex = 1, niceFacing =TRUE)
  }
}, track.height = 0.08, bg.border = NA,bg.col = NA)
#展示图形的圆圈
circos.genomicTrack(bed, ylim = c(0, 1),track.height = 0.15,bg.border=NA,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = value[,1], ybottom = value[,2], col = value[,3],
                                         border = "black", ...)
                      for(j in 1:nrow(value)){
                        if(value[j,4]!=0){
                          circos.genomicText(region[j,], value[j,], y = 0.25, labels = value[j,5], adj=0.5,cex=1,...)
                        }
                      }
                    })

#展示功能和基因的联系
for(i in 1:nrow(gene_GO_data)){
  genei = gene_GO_data$gene[i]
  GOi = gene_GO_data$GO[i]
  circos.link(genei, c(gene_GO_data$start[i], gene_GO_data$end[i]), 
              GOi, c(gene_GO_data$start2[i], gene_GO_data$end2[i]), 
              col = gene_GO_data$linkcol[i])
}
circos.clear()

#绘制右边pvalue的图例
par(mar=c(3,0,3,10))
barplot(rep(1,100),col=pvalue_col_fun(seq(min(pvalue),max(pvalue),length=100)),
        space=0,border=NA,xaxt="n",yaxt="n",main="-log10(pvalue)")
axis(1,c(1,50,100),c(min(pvalue),mean(pvalue),max(pvalue)),tick=F)

#绘制GO名称的图例
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=bed3$id,col=bed3$col,pch=15,pt.cex=3,cex=1.2,bty="n",title="Gene Ontology")
dev.off()

pvalueFilter=0.05     #p值过滤条件
qvalueFilter=1        #矫正后的p值过滤条件

#定义图形的颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
#setwd("D:\\AI\\trans\\breast\\combined\\km_plot\\GO.high")      #设置工作目录
rt=read.csv("verified_down.txt", header=T, sep=",", check.names=F)     #读取输入文件

#获取核心基因的名称, 将基因名字转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="down-KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="down-kegg-barplot.pdf", width=8, height=7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format=100, color = colorSel)
dev.off()

#气泡图
pdf(file="down-kegg-bubble.pdf", width=8, height=7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format=100, color = colorSel)
dev.off()

input="down-KEGG.txt"      #输入文件
outpdf="down_KEGG.circos.pdf"     #输出结果文件
#setwd("C:\\Users\\Dell\\Desktop\\18.KEGGcircos")     #设置工作目录

#读取输入文件
data=read.table(input, header=T, sep="\t", check.names=F)
data=data[order(data$pvalue),]     #根据pvalue对数据进行排序
data=head(data, n=8)      #展示通路的数目(展示富集最显著的前八个通路)
pvalue=round(-log(data$pvalue,10),2)
#获取展示基因的列表
genelist=unlist(strsplit(data$geneID,"/"))
genetable=table(genelist)

#定义图形的颜色
pathway_col = colorRampPalette(brewer.pal(9, "Paired"))(nrow(data))
genes_col = rainbow(length(names(genetable)),s=0.7,v=0.7)
pvalue_col_fun <- colorRamp2(
	breaks = c(min(pvalue), mean(pvalue), max(pvalue)), 
	colors = c("#FFFFCC", "#FF9966", "#FF0000")
)

#获取基因和通路的对应关系
genedict = list()
for(i in names(genetable)) genedict[[i]] = 0
pathwaydict = list()
for(i in data$Description) pathwaydict[[i]] = 0
# link  gene--pathway
gene_pathway_list = list()
n = 0
for(i in 1:nrow(data)){
  genei = strsplit(data$geneID[i],"/")[[1]] # multi
  pathwayi = data$Description[i]
  for(j in 1:length(genei)){
    n = n+1
    genej = genei[j] # single
    gene_pathway_list[[n]] = data.frame(gene=genej,pathway=data$Description[i],start=genedict[[genej]],
                                        end=genedict[[genej]]+1,start2 = pathwaydict[[pathwayi]],end2=pathwaydict[[pathwayi]]+1,
                                        pvalue=pvalue[i])
    genedict[[genej]] = genedict[[genej]]+1
    pathwaydict[[pathwayi]] = pathwaydict[[pathwayi]]+1
  }
}
#gene \t pathway \t start \t end \t pvalue
gene_pathway_data = as.data.frame(do.call('rbind',gene_pathway_list))
gene_pathway_data$linkcol = pathway_col[as.numeric(as.factor(gene_pathway_data$pathway))]
#right pathway
data3 = data.frame(id=data$Description,start = 0, end = data$Count)
# left top
data1 = data.frame(id=names(genetable),start=0,end = as.numeric(genetable))
# main chrom
df = as.data.frame(rbind(data3,data1))

#定义显著性的标记
get_sig = function(p){
  ifelse(p> -log(0.001,10),"***",ifelse(p> -log(0.01,10),'**','*'))
}

bed3 = data.frame(data3,yt=0,yb=1,col=pathway_col[as.numeric(as.factor(data3$id))],p=0,text='')
bed1 = data.frame(data1,yt=0.5,yb=1,col=genes_col[as.numeric(as.factor(data1$id))],p=0,text='')
bed2 = data.frame(id=gene_pathway_data$gene,start=gene_pathway_data$start,
                  end=gene_pathway_data$end,yt=0,yb=0.5,
                  col=pvalue_col_fun(gene_pathway_data$pvalue),p=gene_pathway_data$pvalue,
                  text=get_sig(gene_pathway_data$pvalue))

bed = as.data.frame(rbind(bed1,bed2,bed3))

#绘制图形
pdf(file=outpdf, width=12, height=7)
layout(mat=matrix(c(1,1,1,0,2,3),nc=2),width=c(6.5,3.5),height=c(2,2,7))
#初始化圈图
circos.par(track.margin=c(0.01,0.01), start.degree=90)
circos.genomicInitialize(df,plotType="none")

#在图形中展示基因的名称
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  if(!any(data3$id%in%sector.index)){
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "reverse.clockwise", niceFacing = TRUE)
  }
}, track.height = 0.08, bg.border = NA,bg.col = NA)
#展示图形的圆圈(包括左半圆的基因和右半圆的通路)
circos.genomicTrack(bed, ylim = c(0, 1),track.height = 0.15,bg.border=NA,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = value[,1], ybottom = value[,2], col = value[,3],
                                         border = "black", ...)
                      for(j in 1:nrow(value)){
                        if(value[j,4]!=0){
                          circos.genomicText(region[j,], value[j,], y = 0.25, labels = value[j,5], adj=0.5,cex=1,...)
                        }
                      }
                    })

#展示基因和通路的联系
for(i in 1:nrow(gene_pathway_data)){
  genei = gene_pathway_data$gene[i]
  pathwayi = gene_pathway_data$pathway[i]
  circos.link(genei, c(gene_pathway_data$start[i], gene_pathway_data$end[i]), 
              pathwayi, c(gene_pathway_data$start2[i], gene_pathway_data$end2[i]), 
              col = gene_pathway_data$linkcol[i])
}
circos.clear()

#绘制右边pvalue的图例
par(mar=c(3,0,3,10))
barplot(rep(1,100),col=pvalue_col_fun(seq(min(pvalue),max(pvalue),length=100)),
        space=0,border=NA,xaxt="n",yaxt="n",main="-log10(pvalue)")
axis(1,c(1,50,100),c(min(pvalue),mean(pvalue),max(pvalue)),tick=F)

#绘制通路名称的图例
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=bed3$id,col=bed3$col,pch=15,pt.cex=3,cex=1.2,bty="n",title="KEGG pathway")
dev.off()

