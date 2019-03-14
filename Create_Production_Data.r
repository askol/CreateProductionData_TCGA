## THIS FUNCTION CREATES MRNA AND MIRNA DATASETS THAT ARE PROCESSED TO:
## 1. EXCLUDE NORMAL TISSUE SAMPLES
## 2. AVERAGE SAMPLES THAT ARE FROM THE SAME PHYSICAL SAMPLE
## 3. REMOVE COUNT INFORMATION FROM GENES/MIRNA THAT HAVE TOO FEW READS TO BE
##    USEFUL TO ANALYZE
##
## THIS FUNCTION IS CALLED BY RUN_CREATE_PRODUCTION_DATA.R
##
##
## THE RESULTING DATASETS ARE THEN PRINTED IN THE DIRECTORIES THAT THE ORIGINAL
## DATA ORIGINATED AND GIVEN AN _USE.TXT SUFFIX
##
## RESULTING DATASETS ARE USED BY EXPLORE_EXPRESSION_RESULTS.R AND
## COMPARE_LIMMA_VOOM.R
##

source("/gpfs/data/stranger-lab/askol/TCGA2/Coexpression/Code/Create_Production_Data_funcs.r")

## ------ MAIN --------##
args <- commandArgs(TRUE)
project = args[1]
ExpDir <- args[2]
MetaDir = args[3]
DataDir = args[4]
Create_Production_Data(project, ExpDir, MetaDir, DataDir)

print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))
