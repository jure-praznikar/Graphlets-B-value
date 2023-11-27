## Calculate the Graphlet Degree Vector (GDV)
## This file reads two pdb files: the PDB query file (6dnl.pdb) and the "renamePDB.pdb" file.
## The output is the file "GDV.rds" which contains the degree of orbits for each atom.

#INPUT: renamePDB.pdb, 6dnl.pdb
#OUTPUT: 6dnl.dataGDV

file <-("renamePDB.pdb")
cat('>>> ',file,'\n')
pdb<-read.pdb(file)
indsP <- atom.select(pdb, "protein") ## all chains, also "9" >> crystal contacts
indsH <- atom.select(pdb, "h") 
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb<-trim.pdb(pdb, inds = indsTRIM)
inds <- atom.select(pdb)
my.atoms<-pdb$atom[inds$atom,c("x", "y", "z")]
N<-dim(my.atoms)[1]
# Calculate distances between all atom pairs  
DD<-dist(my.atoms)
# set threshold for graph construction >> 5.0A contacts (or 7.0A)
#threshold<-5 #Upper treshold=5.0Å
threshold<-7 #Upper treshold=7.0Å
DD[DD>threshold]<-0
DD[DD>1.0E-5]<-1     
adj <- matrix(0,nrow = N, ncol = N) 
adj[lower.tri(adj, diag=FALSE)] <- DD
adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)] # upper triangle - matrix 
g<-graph.adjacency(adj,mode="undirected",weighted=NULL)
rm(adj,DD)
go<-count_orbits_per_node(g,max_graphlet_size=4) # count orbits
rm(g)

## All atoms (including crystal contacts) are used to construct the graph and
## determine the degree of orbits.
## Then a new graph is created containing only nodes (atoms) from
## the PDB file (no crystal contacts).
## So pdb+crsytall contacts are used to determine the degree of orbits.
## The smoothing process and prediction only include atoms from
## the PDB file (no crystal contacts). Note that crystal contacts increase
## the degree of orbits belonging to the outer atoms.


# select only chains from BIO 1  
fileBIO<-('6dnl.pdb') ## this is query "protein", protein of interest
BIO<-biounit(read.pdb(fileBIO))  
CHAIN<-unique(BIO[[1]]$atom$chain) # extract BIO 1 chain's
cat(CHAIN,'\n') 
indsTRIM <- atom.select(pdb,chain=CHAIN)
go<-go[indsTRIM$atom,]     
pdb<-trim.pdb(pdb, inds = indsTRIM)
inds <- atom.select(pdb) # re-index
   
BF<-pdb$atom[inds$atom,"b"]
 
my.atoms<-pdb$atom[inds$atom,c("x", "y", "z")]        
DD<-dist(my.atoms)        
N<-dim(my.atoms)[1]
cat('>>>>>> Protein size: ',N,' <<<<<<<\n')
adj <- matrix(0,nrow = N, ncol = N)

# gcc graph is used only for smoothing >> covalent contacts (not for GDV!)
#threshold<-2.0 #threshold=2.0Å, close/covalent bond contact
threshold<-2.15 #treshold=2.15Å, close/covalent bond contact (S-S < 2.1A)
DD[DD>threshold]<-0
DD[DD>1.0E-5]<-1  
adj[lower.tri(adj, diag=FALSE)] <- DD
adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)] # upper triangle - matrix           
gcc<-graph.adjacency(adj,mode="undirected",weighted=NULL)
rm(adj,DD)

# do smoothing, only one cycle, for atoms closer than 2.1A
goSM<-go*0
for (n in 1:N) {
       t<-unlist(adjacent_vertices(gcc, n))
       if (length(t)>1){goSM[n,]<-(go[n,]+colMeans(go[t,]))/2 }
       else if (length(t)==1) { goSM[n,]<-(go[n,]+go[t,])/2 }
       else if (length(t)==0) { goSM[n,]<-go[n,] } 
}
go<-goSM             
#write data frame in to file >> 6dnl.dataGDV
df<-data.frame(BF,go)
myfile<-("6dnl.dataGDV")
write.table(df, file=myfile, sep=",",
            eol="\n", append=FALSE,
            row.names=FALSE, col.names=FALSE)

# df is of size NxM, where N is the number of atoms and M is 16
# the first column is B-values from pdb file, columns from
# 2 to 16, are degrees of orbits (O0, O1, ..., O14)


cat('\n')
cat('**************** Calculation of GDV is complete. (see file 6dnl.dataGDV) ******************* \n')
cat('\n')


