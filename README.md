# Graphlets-B-value-prediction
Using graphlet degree vector to predict atomic displacement parameters in protein structures


# List of libraries to check
libraries <- c("quantmod", "ggplot2", "tidyr", "lubridate")

# Check if each library is installed
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    cat("\n", lib, "library is about to get installed!\n")
    
    # If the library is not installed, install it
    install.packages(lib)
    # Load the library after installation
    library(lib, character.only = TRUE)
  }
  else{
    cat("\n", lib, "library is already installed!\n")
  }
}


## library needed, check if exists!!

Save the query PDB file in folder EXAMPLE1. For example go to
https://www.rcsb.org/structure/6DNL and download odb file (6dnl.pdb)

1) Open next link / website. Select pdb file (6dnl.pdb) from EXAMPLE1 folder.
https://swift.cmbi.umcn.nl/servers/html/symshel2.html

2) Wait for the server to finish the job and save the output file "symtry.pdb"
in to folder EXAMPLE1. In your folder are now should two pdb files. 
Query file (6dnl.pdb) and symtry.pdb

3) run R script
Go to folder EXAMPLE1 and run R script  "symetry.r"
>Rscript symetry.r

This script renames symmetry chains to name "9". Symmetry atoms/chains
are used only for the calculation of graphlet degree vectors, 
but not for "smoothing" and prediction.

The R script "symmetry.r" produces the file "fixed.pdb", this file
is further used by the next R script, see step (4).

4) run R script "work_GDV_for_fixed_sym.r"
> Rscriprt work_GDV_for_fixed_sym.r 

This file creates reads two pdb files: (i) query PDB file (6dnl.pdb) and
file that contains symmetry atoms "fixed.pdb".  
The output is "GDV.rds" file, which contains normalized orbits for
each atom. This part is used to predict normalized and
rescaled B-values. Rescaled B-values are calculated if
Bvalues in PDB file are available. Also, constant B-values can not be used for rescaling

5) predict, run the following script, "predict.r"
> Rscript ....  

This file predicts normalized B-values, which are written in the file
predicted.pdb.
Also, the correlation between the deposited (PDB file) and a predicted
file is calculated. The predicted PDB file does not contain hydrogen atoms,
HETERO atoms, solvent atoms, and only protein atoms are included.

>> Lahko bi dodal RESCALED > da se vzame povprecje iz PDB fajla in
tudi standardno deviacijo, ter se izracuna B-values!!!
