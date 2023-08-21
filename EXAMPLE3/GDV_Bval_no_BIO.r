## 1) read PDB file, convert 3D atom model in to graph
## 2) calculate Graphlet Degree Vector
## 3) predict B-values

#INPUT: 6SK0.pdb
#OUTPUT: 6SK0.dataGDV, predicted.pdb, rescaled.pdb

pdb<-read.pdb('6SK0.pdb')
indsP <- atom.select(pdb, "protein") ## all chains, also "9", sym
indsH <- atom.select(pdb, "h") 
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb<-trim.pdb(pdb, inds = indsTRIM)
inds <- atom.select(pdb)
#
BF<-pdb$atom[inds$atom,"b"]
my.atoms<-pdb$atom[inds$atom,c("x", "y", "z")]
N<-dim(my.atoms)[1]
# calculate distnces between all atom pairs  
DD<-dist(my.atoms)
# set treshold for graph construction >> 5.0A contacts (or 7.0A)
#threshold<-5 #Upper treshold=5.0Å
threshold<-7 #Upper treshold=7.0Å
DD[DD>threshold]<-0
DD[DD>1.0E-5]<-1     
adj <- matrix(0,nrow = N, ncol = N) 
adj[lower.tri(adj, diag=FALSE)] <- DD
adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)] # upper matrix 
g<-graph.adjacency(adj,mode="undirected",weighted=NULL)
go<-count_orbits_per_node(g,max_graphlet_size=4) # count orbits
rm(g) 
# gcc graph is used only for smoothing >> cavalent contacts (not for GDV!)
#threshold<-2.0 #treshold=2.0Å, close/covalent bond contact
threshold<-2.15 #treshold=2.15Å, close/covalent bond contact (S-S < 2.1A)
DD[DD>threshold]<-0
DD[DD>1.0E-5]<-1  
adj[lower.tri(adj, diag=FALSE)] <- DD
adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)] # upper matrix           
gcc<-graph.adjacency(adj,mode="undirected",weighted=NULL)
rm(adj,DD) 
# do smoothing, only one round, for atoms closer than 2.0A
goSM<-go*0
for (n in 1:N) {
       t<-unlist(adjacent_vertices(gcc, n))
       if (length(t)>1){goSM[n,]<-(go[n,]+colMeans(go[t,]))/2 }
       else if (length(t)==1) { goSM[n,]<-(go[n,]+go[t,])/2 }
       else if (length(t)==0) { goSM[n,]<-go[n,] } 
}
go<-goSM # Graphlet Degree Vector after smoothing             

#write data frame in to file >> 6dnl.dataGDV
df<-data.frame(BF,go)
myfile<-("6SK0.dataGDV")
write.table(df, file=myfile, sep=",",
            eol="\n", append=FALSE,
            row.names=FALSE, col.names=FALSE)

# df is of size NxM, where N is numer of atoms and M is 16
# first column is B-values from pdb file, columns from
# 2 to 16, are degree of orbits (O0, O1, ..., O14)

## Predict the B values
## The predicted B-values are written to two files:
##  -predicted.pdb (it contains the predicted normalized B-values)
##  -rescaled.pdb (it contains the predicted rescaled B-values)
## Note that rescaled.pdb is only created if the standard deviation and
## mean B-value of the query structure are in the range [2, 100].

data<-scale(go) # Graphlet Degree Vector for each atom
#BF<-temp[,1] 
beta0 <- 0 # intercept is zero
# this are regression coeffcients, for details see equation (3) in
# following publication 
beta <- c(3.320838e-01,-2.475571e+00,2.981908e-01,-1.300847e+00,
          -7.121716e-01,1.173306e+00,3.489543e-01,-2.750241e-01,
          1.046642e-01,4.966259e-01,6.795167e-01,-8.367222e-02,
          4.285069e-01,3.710838e-02,2.792701e-01)
data<-as.matrix(data)
BFPR <- beta0 + data %*% beta # predicted B-values (linear model)
CC<-cor(BF,BFPR)
cat('\n')
cat('------------------------------\n')
cat('Correlation (predicted vs pdb): ',round(CC,2),'\n')
cat('------------------------------\n')
cat('\n')
### write predicted BF's in to PDB file
pdb1<-pdb
pdb1$atom$b<-BFPR
write.pdb(pdb1, file = "predicted.pdb")

## RESCALE B-values and write predicted-rescaled BF's in to PDB file
## first check mean B val and standard deviation of B
if (mean(BF)>2 & mean(BF)<100 & sd(BF)>2 & sd(BF)<100) {
    #cat('Rescale B-values')
    BFrescaled<-BFPR*sd(BF)+mean(BF)
    pdb1$atom$b<-BFrescaled
    write.pdb(pdb1, file = "rescaled.pdb")
}


cat('\n')
cat('**************** Prediction is complete. (see file predicted.pdb and rescaled.pdb) ******************* \n')
cat('\n')
