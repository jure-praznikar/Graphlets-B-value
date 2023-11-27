# **Introduction**  
The components of the graphlet degree vector, which describes the complexity<br> 
of the wiring of a given atom, can be used in a multiple linear regression model <br>
to predict atomic displacement parameters in protein structures. <br>
For details see:<br>
J. Praznikar (2023). Acta Cryst. 79, https://doi.org/10.1107/S2059798323009142 <br>
By default, Biological unit 1 is selected from the pdb file.  <br>
Otherwise, all chains in the asymmetric unit are selected. <br> 

### **Prerequisites for running R scripts**  

**install packages**  
* install.packages("igraph")  
* install.packages("bio3d")
* install.packages("pdist")
* install.packages("remotes")
* remotes::install_github("alan-turing-institute/network-comparison")

**load libraries**  
* library(igraph)  
* library(bio3d)
* library(pdist)
* library(netdist)  

## **Example 1 – Biological unit 1 and crystal contacts**  

By default, the analysis is performed using the Biological Assembly 1 (BIOMOLECULE: 1).  
The data about CRYSTALLOGRAPHIC SYMMETRY must be represented in the PDB file.  
The EXAMPLE1 folder contains the pdb file “6dnl.pdb”.

**_STEP 1: Add Crystal contacts_**  
>source("add_symm_residues.r")

This R script reads 6dnl.pdb and adds crystall contacts (cutoff is 8Å). The output<br>
is written in "6dnl_CRYST.pdb".

**_STEP 2: Calculate the Graphlet Degree Vecotr (GDV): run the R script "GDV.r"_**    
>source("GDV.r")

This R script reads two pdb files: the PDB query file (6dnl.pdb) and the "6dnl_CRYST.pdb" file.  
The output is the file "GDV.rds" which contains the degree of orbits for each atom.  

**_STEP 3: Predict the B values by running the R script "predict.r"_**  
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
R script reads file 6SK0_noREMARK.pdb (it does not contain info about    
Biological Assembly and crystallographic symmetry).

> source("GDV_Bval_no_BIO.r")

or  

> Rscript GDV_Bval_no_BIO.r

The correlation between the predicted and PDB B-values for this particular example should be 0.61.<br>
