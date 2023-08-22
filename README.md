# **Introduction**  
The components of the graphlet degree vector, which describes the complexity<br> 
of the wiring of a given atom, can be used in a multiple linear regression model <br>
to predict atomic displacement parameters in protein structures. <br>
For details see: (TO DO!)<br>
By default, Biological unit 1 is selected from the pdb file.  <br>
Otherwise, all chains in the asymmetric unit are selected. <br> 

### **Prerequisites for running R scripts**  

**install packages**  
* install.packages("igraph")  
* install.packages("bio3d")  
* install.packages("remotes")
* remotes::install_github("alan-turing-institute/network-comparison")

**load libraries**  
* library(igraph)  
* library(bio3d)  
* library(netdist)  

## **Example 1 – Biological unit 1 and crystal contacts**  

By default, the analysis is performed using the Biological Assembly 1 (BIOMOLECULE: 1).  
The data about CRYSTALLOGRAPHIC SYMMETRY must be represented in the PDB file.  

**_STEP 1: Crystal contacts_**  

The EXAMPLE1 folder contains the pdb file “6dnl.pdb”.  
Go to https://swift.cmbi.umcn.nl/servers/html/symshel2.html and upload the 6dnl.pdb file (from the EXAMPLE1 folder).  

Wait for the server to finish its work and save the output file "symtry.pdb" in the EXAMPLE1 folder.   
There should now be two pdb files in your folder: 6dnl.pdb and symtry.pdb.   
The symtry.pdb file contains the crystallographic contact atoms.  

**_STEP 2: Rename symmetry chains_**     
Go to the EXAMPLE1 folder and run the R script "renameSYM.r"  
>source("renameSYM.r")    

This script renames symmetry chains to the name "9".   
Symmetry atoms/chains are used only for graphlet degree vector calculation, not for "smoothing" and prediction.  
The renaming is necessary to distinguish between the atoms from the PDB file and the symmetry atoms/chains.  
Usually letters are used for the chain ID's, so the (number) "9" was used instead of letters.    

The R script "renameSYM.r" creates the file "renamePDB.pdb".   
This file contains the renamed chains and will be used further, see step (4).  

**_STEP 3: Calculate the Graphlet Degree Vecotr (GDV): run the R script "GDV.r"_**    
>source("GDV.r")  

This file reads two pdb files: the PDB query file (6dnl.pdb) and the "renamePDB.pdb" file.  
The output is the file "GDV.rds" which contains the degree of orbits for each atom.  

**_STEP 4: Predict the B values by running the R script "predict.r"_**  
>source("predict.r")  

The correlation between the predicted and PDB B-values for this particular example should be 0.71.<br>

The predicted B-values are written to two files:  
* predicted.pdb (it contains the predicted normalized B-values)  
* rescaled.pdb (it contains the predicted rescaled B-values)  

Note that rescaled.pdb is only created if the standard deviation and <br>
mean B-value of the query structure are in the range [2, 100].

## **Example 2 – Biological unit 1, no crystal contacts (i.e., cryo-EM model)**
This example does not consider symmetry contacts (i.e., the cryo- EM model).  
By default, the analysis is performed using the Biological Assembly 1 (BIOMOLECULE: 1).  

**_STEP 1: Predict the B values by running the R script "GDV_Bval.r"_**  

> source("GDV_Bval.r")

or  

> Rscript GDV_Bval.r 

The script “GDV_Bval.r” reads the file 6SK0.pdb, which is already in the EXAMPLE2 folder.  
Then it converts the 3D protein model into a graph and calculates the Graphlet Degree Vector for each atom.  
The predicted B-values are written to two files:  
* predicted.pdb (it contains the predicted normalized B-values)  
* rescaled.pdb (it contains the predicted rescaled B-values)

The correlation between the predicted and PDB B-values for this particular example should be 0.61.<br>
   
Note that rescaled.pdb is only created if the standard deviation and  
mean B-value of the query structure are in the range [2, 100].  

## **Example 3 -  Asymmetric unit, no crystal contacts**
This script is a slightly modified version of “Example 2”.  
The only difference is that the script reads all the atoms of the protein   
in the PDB file (asymmetric unit) rather than in Biological Assembly 1.   

> source("GDV_Bval_no_BIO.r")

or  

> Rscript GDV_Bval_no_BIO.r

The correlation between the predicted and PDB B-values for this particular example should be 0.61.<br>
