## Predict the B values
## The predicted B-values are written to two files:
##  -predicted.pdb (it contains the predicted normalized B-values)
##  -rescaled.pdb (it contains the predicted rescaled B-values)
## Note that rescaled.pdb is only created if the standard deviation and
## mean B-value of the query structure is in the range [2, 100].

# INPUT: 6dnl.pdb, 6dnl.dataGDV
# OUTPUT: predicted.pdb, rescaled.pdb

# read graphlet degree vector GDV (for example 6dnl.dataGDV)
temp<-read.table('6dnl.dataGDV', header = FALSE, sep = ",") 

data<-scale(temp[,2:16]) 
BF<-temp[,1] # B-values from pdb file 
beta0 <- 0 # intercept is zero
# These are regression coefficients, for details see equation (5) in
# following publication: J. Praznikar (2023). Acta Cryst. 79, https://doi.org/10.1107/S2059798323009142
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
## read pdb file (6dnl.pdb)
pdb<-read.pdb('6dnl.pdb') # query file
BIO <- tryCatch( { biounit(pdb) }, error = function(e){} )
pdb<-BIO[[1]]
cat(paste(attributes(BIO)),'\n')
indsP <- atom.select(pdb, "protein")
indsH <- atom.select(pdb, "h") 
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb<-trim.pdb(pdb, inds = indsTRIM)

### Write predicted B-values into PDB file
pdb1<-pdb
pdb1$atom$b<-BFPR
write.pdb(pdb1, file = "predicted.pdb")

## RESCALE B-values and write predicted-rescaled BFs into PDB file
## First check mean B val and standard deviation of B
if (mean(BF)>2 & mean(BF)<100 & sd(BF)>2 & sd(BF)<100) {
    #cat('Rescale B-values')
    BFrescaled<-BFPR*sd(BF)+mean(BF)
    pdb1$atom$b<-BFrescaled
    write.pdb(pdb1, file = "rescaled.pdb")
    }

cat('\n')
cat('**************** Prediction is complete. (see file predicted.pdb and rescaled.pdb) ******************* \n')
cat('\n')


