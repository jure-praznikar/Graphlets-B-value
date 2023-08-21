## This script renames symmetry chains to the name "9".
## Symmetry atoms/chains are used only for graphlet degree vector calculation,
## not for "smoothing" and prediction.
## The renaming is necessary to distinguish between the atoms from the PDB file
## and the symmetry atoms/chains.
## The usual names for the chain ID are letters,
## so the name "9" (number) was used instead of letters.

# INPUT: 6dnl.pdb, symtry.pdb
# OUTPUT: renamePDB.pdb

# read pdb file 6dnl.pdb
file<-'6dnl.pdb'
pdb1<-read.pdb(file)
# select "heavy" protein atoms (H atoms excluded)
indsP <- atom.select(pdb1, "protein")
indsH <- atom.select(pdb1, "h") 
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb1<-trim.pdb(pdb1, inds = indsTRIM)
inds1 <- atom.select(pdb1)
# coordinates of "heavy" protein atoms
my.atoms1<-pdb1$atom[inds1$atom,c("x", "y", "z")]  
#
# read pdb file symtry.pdb (of 6dnl.pdb) - crystal contacts!
SYMfile<-'symtry.pdb'
pdb2<-read.pdb(SYMfile)
# select "heavy" protein atoms - crystal contacts!
indsP <- atom.select(pdb2, "protein")
indsH <- atom.select(pdb2, "h") 
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb2<-trim.pdb(pdb2, inds = indsTRIM)
inds2 <- atom.select(pdb2) # indeces in symtry.pdb
# coordinates of "heavy" protein atoms - crystal contacts!
my.atoms2<-pdb2$atom[inds2$atom,c("x", "y", "z")]

#find coordinates that exist in PDB (6dnl.pdb), but are not in symtry.pdb
idx<-numeric(0)
for (j in 1:dim(my.atoms1)[1]){
      X<-abs(my.atoms1[j,1]-my.atoms2[,1])
      Y<-abs(my.atoms1[j,2]-my.atoms2[,2])
      Z<-abs(my.atoms1[j,3]-my.atoms2[,3])
      n<-which(X<1.0E-3 & Y<1.0E-3 & Z<1.0E-3)
      idx<-c(idx,n) # this is in PDB
}
    
inds3<-inds2
# xyz are write in 1D array: (idx-1)*3+1=x, (idx-1)*3+2=y,  (idx-1)*3+3)=z 
inds3$xyz<-inds3$xyz[-c((idx-1)*3+1,(idx-1)*3+2,(idx-1)*3+3)] # xyz in symtry.pdb but not in 6dnl.pdb
inds3$atom<-inds3$atom[-idx] # index which are in symtry.pdb but not in 6dnl.pdb
pdb3<-trim.pdb(pdb2, inds = inds3)
pdb3$atom$chain="9" # rename chain to discriminate atoms/chains (chain ID="9" indicate symmetry atoms )
pdb4<-cat.pdb(pdb1,pdb3,rechain=FALSE) # merge PDB and renamed symetry atoms
# write file
FIXEDfile<-('renamePDB.pdb')
write.pdb(pdb4, file=FIXEDfile)
cat('\n')
cat('**************** Renaming is complete. (see file renamePDB.pdb) ******************* \n')
cat('\n')

#In the renamePDB.pdb file, the atoms belonging to the crystal contacts have chain id="9".
#The numeric name "9" was used because the original pdb file may already
#contain A,B,C, ..., X,Y,Z Id's. See below, for an example.

#   ATOM    912  CA  ARG A 115       2.394 -12.301   5.698  1.00 61.06           C  
#   ATOM    913  C   ARG A 115       1.883 -10.936   6.111  1.00 67.52           C  
#   ATOM    914  O   ARG A 115       2.607  -9.949   5.996  1.00 75.09           O  
#   ATOM    915  CB  ARG A 115       1.347 -13.383   6.004  1.00 62.38           C  
#   ATOM    916  CG  ARG A 115       1.159 -13.717   7.484  1.00 59.57           C  
#   ATOM    917  CD  ARG A 115      -0.078 -14.585   7.721  1.00 58.02           C  
#   ATOM    918  NE  ARG A 115      -1.270 -13.969   7.150  1.00 58.19           N  
#   ATOM    919  CZ  ARG A 115      -1.811 -14.334   5.995  1.00 61.13           C  
#   ATOM    920  NH1 ARG A 115      -1.271 -15.329   5.300  1.00 61.13           N  
#   ATOM    921  NH2 ARG A 115      -2.890 -13.711   5.544  1.00 66.95           N  
#   ATOM    922  N   LEU 9  16       6.877  11.599  -9.921  1.00 16.52           N  
#   ATOM    923  CA  LEU 9  16       7.024  10.147  -9.936  1.00 17.03           C  
#   ATOM    924  C   LEU 9  16       8.457   9.721  -9.644  1.00 15.73           C  
#   ATOM    925  O   LEU 9  16       8.678   8.676  -9.014  1.00 22.18           O  
#   ATOM    926  CB  LEU 9  16       6.564   9.597 -11.288  1.00 12.10           C 

