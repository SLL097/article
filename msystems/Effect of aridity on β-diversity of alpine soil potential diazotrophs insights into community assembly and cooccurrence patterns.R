Computational code
###Sparcc network Analysis Code
####linux code
### otu_data.txt: otu table 
conda activate Fastspar##Activate Environment
fastspar --threshold 0.3 --otu_table otu_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv# Calculate correlation matrix

mkdir bootstrap_counts
fastspar_bootstrap --otu_table otu_data.txt --number 1000 --prefix bootstrap_counts/fake_data# Generate Random Matrix

mkdir bootstrap_correlation
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 10 ::: bootstrap_counts/*# Calculate the correlation of random matrices 
  fastspar_pvalues --otu_table otu_data.txt --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv #Calculate P-value

#### R Language code
cor_sparcc <- read.delim('median_correlation.tsv', row.names = 1, sep = '\t', check.names = FALSE)
cor_sparcc1 <- cor_sparcc

pvals <- read.delim('pvalues.tsv', row.names = 1, sep = '\t', check.names = FALSE)

#Preserve: Correlation  ≥ 0.3 and Values p<0.01
cor_sparcc[abs(cor_sparcc)<=0.3] <- 0


pvals[pvals>=0.05] <- -1
pvals[pvals<0.05 & pvals>=0] <- 1
pvals[pvals==-1] <- 0


adj <- as.matrix(cor_sparcc) %*% as.matrix(pvals) 

write.table(data.frame(adj, check.names = FALSE), 'neetwork.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

library(igraph)
neetwork_adj <- adj

neetwork_adj <- read.delim('neetwork.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)

g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    

E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)


adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))

#Graphml format, can be opened and visually edited using gephi software
write.graph(g, 'network.graphml', format = 'graphml')



####icamp community assembly 
#This is a modified version of the scripts: https://github.com/DaliangNing/iCAMP1
#OTU_ Table: OTU/ASV or group abundance table, with rows being OTUID and columns being sample names
#Taxonomy_ Table: Species classification table, row is OTUID, column is boundary to species name
#Sample_ Info: Sample grouping table, where rows are sample names and columns are processing/grouping names
#Env: Environmental factor table, where the row represents the sample name and the column represents each environmental factor

rm(list=ls())
com.file="OTU_table.txt"
tree.file="OTU_phylo.tre"
clas.file="Taxonomy_table.txt"
treat.file="Sample_info.txt"
env.file="Env.txt"
save.wd="~/total"
if(!dir.exists(save.wd)){dir.create(save.wd)}
nworker=50
memory.G=500 
prefix="Test" 
rand.time=999 

library(iCAMP)
library(ape)
setwd(wd)
comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
tree=read.tree(file = tree.file)
clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table(treat.file, header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)
env=read.table(env.file, header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE) 

sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
env=sampid.check$env 
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
}else{
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}


setwd(save.wd)
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)
ds = 0.25
bin.size.limit = 48 
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)
sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=5 
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)



bin.size.limit = 48
sig.index="Confidence"
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
icbin=icamp.bins(icamp.detail = icres$detail,treat = treat,
                 clas=clas,silent=FALSE, boot = TRUE,
                 rand.time = rand.time,between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) 
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)
i=1
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)


### neutral model
This is a modified version of the scripts: https://github.com/Weidong-Chen-Microbial-Ecology/Stochastic-assembly-of-river-microeukaryotes/blob/master/Neutral%20community%20model.r

rm(list=ls())
setwd("D:\\data\\ neutral model ")

library(Hmisc)
library(minpack.lm)
library(stats4)


#spp: Abundance table of species or taxa, with rows representing taxa and columns representing samples
spp<-read.delim('otu.txt',head=T,stringsAsFactors=F,row.names=1,sep = "\t")

spp<-t(spp)

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr 

#output
write.csv(p, file = "otu_p.csv")
write.csv(freq, file = "otu_freq.csv")
write.csv(freq.pred, file = "otufreq.pred.csv")

#p: mean relative abundance
#freq : observed occurrence frequency
#freq.pred :predict occurrence frequency
