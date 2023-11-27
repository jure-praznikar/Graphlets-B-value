# Add shell of symmetry-related residues
# This script adds crystal contacts - residues to the original pdb file.
# The original PDB file you downloaded from https://www.rcsb.org/
# must contain information about the symmetry and cell size,
# see lines below - example from PDB header.
# The default cutoff size for contacts is 5Å plus 3Å (=approx. 2*VdW radii = 3Å),
# so the total cut-off distance is 8Å.

# REMARK 290 CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS                            
# REMARK 290 THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM             
# REMARK 290 RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY                
# REMARK 290 RELATED MOLECULES.                                                   
# REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
# REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
# REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
# REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000        0.00000            
# REMARK 290   SMTRY2   2  0.000000  1.000000  0.000000       14.03000            
# REMARK 290   SMTRY3   2  0.000000  0.000000 -1.000000        0.00000  
          
#   .
#   .
#   .
#   .
#   .
# CRYST1   43.250   28.060   45.800  90.00 100.96  90.00 P 1 21 1      2 

# INPUT: 6dnl.pdb
# OUTPUT: 6dnl_CRYST.pdb

rm(list=ls())
#
pdbfile<-'6dnl.pdb' # input
outputfile<-'6dnl_CRYST.pdb' #output

#pdb<-read.pdb(pdbfile)
invisible(capture.output(pdb<-read.pdb(pdbfile))) # invisible - no messages 
# select protein atoms
idx<-atom.select(pdb,"protein")
pdb<-trim.pdb(pdb,inds=idx)
pdbORG<-pdb #store for later use to write final pdb file > write.pdb
# CLEAN it > to have a unique residue number overall chains
pdb<-clean.pdb(pdb, consecutive = TRUE,force.renumber = TRUE)
pdb$atom$chain[1:length(pdb$atom$chain)]<-"9" # set chain name to "9" >> discriminate from non-symm resi.
# select CB atoms - the first stage of searching close contacts
idxCB<-c(which(pdb$atom$elety=="CB"),which(pdb$atom$resid=="GLY" & pdb$atom$elety=="CA"))
# geometric center
gc<-colMeans(pdb$atom[,c("x","y","z")])
## max deviation in x, y, z
dxyz<-numeric(3)
dxyz[1]<-max(abs(gc[1]-pdb$atom[,c("x")]))
dxyz[2]<-max(abs(gc[2]-pdb$atom[,c("y")]))
dxyz[3]<-max(abs(gc[3]-pdb$atom[,c("z")]))
# d is max distance between the geometric center and the most distant atom
d<-sqrt(sum(dxyz^2))
#
VdW<-3 # approx. two times Van der Waals radii (3Å)
cutoff<-5+VdW # distance  = 5Å
# if the distance between the centers of the two molecules is
# more than 2*d+cutoff, then are not in contact
dc<-2*d+cutoff
#
xyz0 <- rbind(matrix(pdb$xyz, nrow=3), 1) # initial coordinates as vector
xyz30<-t(matrix(pdb$xyz,nrow=3)) # initial coordinates as Nx3 matrix
# extract symmetry data and Cell size (the data must be PDB header)
temp<-readLines(pdbfile,n=1000)
k<-grep('CRYST1',temp)
A<-unlist(strsplit(temp[k],split = ' '))
A<-A[A != ""]
CELL<-as.numeric(A[2:4])
ANGLE<-as.numeric(A[5:7])
k<-grep('SMTRY',temp)
temp<-temp[k]
s1<-seq(1,length(temp),3)
# create crystallographic  matrices (rotation and translation) - data from PDB header
M<-list()
for (i in 1:(length(temp)/3)){
     mat0<-matrix(0,ncol=4,nrow=3)
     k<-i*3-2
     #1
     A<-unlist(strsplit(temp[k],split = ' '))
     A<-A[A != ""]
     mat0[1,]<-as.numeric(A[5:8])
     #2
     A<-unlist(strsplit(temp[k+1],split = ' '))
     A<-A[A != ""]
     mat0[2,]<-as.numeric(A[5:8])
     #3
     A<-unlist(strsplit(temp[k+2],split = ' '))
     A<-A[A != ""]
     mat0[3,]<-as.numeric(A[5:8])
     # append
     M<-append(M,list(mat0))
}
#

deg2rad <- function(deg) {(deg * pi) / (180)}
a<-deg2rad(ANGLE[1]) # alpha
b<-deg2rad(ANGLE[2]) # beta
g<-deg2rad(ANGLE[3]) # gamma
temp<-(cos(b)*cos(g)-cos(a))/(sin(b)*sin(g))
a_st<-acos(temp)
#convert CELL to orthogonal 
#https://www.ucl.ac.uk/~rmhajc0/frorth.pdf
# equation 14,15, 16, from #https://www.ucl.ac.uk/~rmhajc0/frorth.pdf
#X <- xi*CELL[1] + yi*CELL[2]*cos(g) + zi*CELL[3]*cos(b)
#Y <- yi*CELL[2]*sin(g) - zi*CELL[3]*sin(b)*cos(a_st)
#Z <- zi*CELL[3]*sin(b)*sin(a_st)

mat<-matrix(0,ncol=4,nrow=3)
PDBlist<-list() # store contacts
# -2,-1,0,1,2 all combinations in x, y,z
myshift<-seq(-2,2) # shift along each axis (only translation)

#start_time <- Sys.time()
for (i in 1:length(M)) {
     mat<-M[[i]]
     for (xi in myshift) {
       for (yi in myshift) {
         for (zi in myshift) {
           #to orthogonal   
           mat[1,4] <- M[[i]][1,4] + ( xi*CELL[1] + yi*CELL[2]*cos(g) + zi*CELL[3]*cos(b) )
           mat[2,4] <- M[[i]][2,4] + ( yi*CELL[2]*sin(g) - zi*CELL[3]*sin(b)*cos(a_st) )
           mat[3,4] <- M[[i]][3,4] + ( zi*CELL[3]*sin(b)*sin(a_st) )        
           #init pdb
           pdb1<-pdb
           pdb2<-pdb                     
           #xyz <- matrix(mat %*% xyz0, nrow = 1)
           pdb2$xyz<-matrix(mat %*% xyz0, nrow = 1)
           pdb2ALL<-pdb2
           xyzT<-t(matrix(pdb2$xyz,nrow=3))
           XYZgc<-colMeans(xyzT)
           test<-dist(rbind(gc,XYZgc))
           #if the distance between the centers of the two molecules is more than 2*d+cutoff >> no contact
           if ( test>0.5 & test<dc ) {
              MAT <- as.matrix(pdist( xyz30[idxCB,], xyzT[idxCB,] ))
              v<-which(MAT<=(cutoff+7), arr.ind = TRUE) # CB atoms cutoff+7Å
              v1<-unique(v[,1])
              v2<-unique(v[,2])
              if (length(v2)>=1){
                  pdb1<-trim.pdb(pdb1,resno=pdb1$atom$resno[idxCB[v1]])
                  pdb2<-trim.pdb(pdb2,resno=pdb2$atom$resno[idxCB[v2]])
                  idx2<-atom.select(pdb2)
                  # On selected residues perform the pairwise calculation for all atoms (a subset of residues)
                  MAT <- as.matrix(pdist( t(matrix(pdb1$xyz,nrow=3)), t(matrix(pdb2$xyz,nrow=3)) ))
                  v<-which(MAT<=cutoff, arr.ind = TRUE) # cutoff
                  v<-unique(v[,2])
                  if (length(v2)>=1){
                      #cat(paste('OUT_',xi,yi,zi,'sym',i-1,'.pdb',sep=''),'\n')
                      cat(paste('> ',xi,yi,zi,' sym ',i-1,sep=''),'\n')
                      #write.pdb(pdb2ALL,file=paste('OUT_',xi,yi,zi,'sym',i-1,'.pdb',sep=''))
                      pdb3<-trim.pdb(pdb2,resno=unique(pdb2$atom$resno[idx2$atom[v]]))
                      PDBlist<-append(PDBlist,list(pdb3))
                  }
              }
           }       
       }   
      }
     }         
}    
#end_time <- Sys.time()
#cat(end_time - start_time,'\n')

#### Write pdb file
for ( f in 1:length(PDBlist)) {      
     #pdbORG <- cat.pdb(pdbORG,PDBlist[[f]],rechain=FALSE)
     suppressWarnings(pdbORG <- cat.pdb(pdbORG,PDBlist[[f]],rechain=FALSE))
#    invisible(capture.output(pdbORG <- cat.pdb(pdbORG,PDBlist[[f]],rechain=FALSE)))
    }
write.pdb(pdbORG,file=outputfile)
