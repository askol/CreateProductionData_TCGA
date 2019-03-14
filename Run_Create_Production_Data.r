job.files = c()
setwd("/gpfs/data/stranger-lab/askol/TCGA2/Data/Production")

projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt", header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 

pbsOutDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/Production/"
CodeDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Coexpression/Code/"
ExpDir <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Expression/"
MetaDir <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_MetaData/"
DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/"


for (project in projects){

    job.file = paste0(project,".pbs")
    job.files = c(job.files, job.file)

    write(file = job.file, "#!/bin/bash")
    write(file = job.file, paste("#PBS -l nodes=1:ppn=1,mem=8gb"), append=TRUE)
    write(file=job.file, "#PBS -l walltime=36:00:00", append=TRUE)
    write(file = job.file,
          paste0("#PBS -o ", pbsOutDir,  project,"_pbs.out"), append=TRUE)
    write(file = job.file, "#PBS -j oe", append=TRUE)
    write(file = job.file, paste0("#PBS -N" ,project), append=TRUE)
    write(file = job.file, paste0("R CMD BATCH  '--args ",project, " ", ExpDir, " ",
              MetaDir, " ",  DataDir,  " ' ",CodeDir,
              "Create_Production_Data.r ",
              pbsOutDir, project,".out"),  append=TRUE)    
    system(paste("chmod +x", job.file))
}

for (job in job.files){
    system(paste("qsub",job))
}
