system("module load htslib")
system("module load bcftools")
system("module load bedtools")

library(edgeR)
library(limma)
library(stringr)
library(readr)
library(SmartSVA)
library(biomaRt)
library(dplyr)
library(annotables)

ExpDir <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Expression/"
MetaDir <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_MetaData/"
DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/"

Create_Production_Data <- function(project, ExpDir, MetaDir, DataDir){

    make.expr.files(project, ExpDir, MetaDir, DataDir)
##    make.mirna.files(project)
    
}

make.expr.files<- function(project,
                           ExpDir, MetaDir, DataDir,
                           min.reads=5,
                           min.samples=.02,
                           use.existing.gene.info = TRUE)
{
    
    ## MUCH OF WHAT IS BELOW IS DERIVED FROM:
    ## https://www.bioconductor.org/help/workflows/RNAseq123/
    ## setwd("/gpfs/data/stranger-lab/askol/TCGA/TCGA_Expression/")

    ## GET GENE INFO ##
    gene.info.file <- paste0(DataDir, "GeneInfo/Gene.Info.txt")
    ## GET GENE INFORMATION FROM BIOMART ##
    if (use.existing.gene.info==TRUE){
        gene.info <- read.table(file = gene.info.file, as.is=T, header=T, sep="\t", quote="")
        print("Using existing gene.info file")
    }else{
        gene.info <- Get_Gene_Info()
        write.table(file = gene.info.file, gene.info, quote=F, row.names=F,
                    col.names=T, sep="\t")
    }  
    
    ## GET SAMPLE INFORMATION AND COVARIATES (CASE.ID, GENDER)
    covs <- get.covs(project, MetaDir)

    ## GET MRNA COUNT AND FPKM DATA ##
    countFile <- paste0(ExpDir,"/Count/",project,"/",project, "-Count.txt")

    out <- get.m.data(countFile, covs, gene.info, min.reads=5, min.samples=.02)
    
    CountData <- out$CountData
    covs <- out$covs       
    
    ## REMOVE SAMPLES WITH MISSING VALUES OF THE COV ##
    
    ind.na <- which(is.na(covs$gender))
    birthyear.na <- is.na(covs$BirthYear)
    if (mean(birthyear.na) < 0.10){
        ind.na <- unique(c(ind.na, which(birthyear.na)))
    }
    ids.na <- covs$cases.0.case_id[ind.na]
    
    print(paste0("Removing ",length(ids.na), " cases with no gender or BirthYear information"))
    
    if (length(ind.na) > 0){
        covs <- covs[-ind.na,]
        ind.na = which(gsub("A","", names(CountData)) %in% ids.na)
        CountData = CountData[,-ind.na]
    }
    covs$case.id <- paste0("A",covs$cases.0.case_id)

    ## KEEP ONLY NEEDED COVS
    covs <- covs[, c("case.id", "gender", "BirthYear")]
    
    gender <- covs[,c("case.id", "gender")]
    rm.dupes <- which(duplicated(gender$case.id))
    if (length(rm.dupes) > 0){
        gender <- gender[-rm.dupes,]
    }
    group <- data.frame(case.id = names(CountData)[-1])
    group <- merge(group, gender, by.x=1,
                   by.y="case.id", all.x = TRUE, all.y=FALSE, sort=F)
    rownames(CountData) <- CountData$ensmbl 

    ## MATCH UP GENE.INFO AND COUNTFILE BY ENSEMBL ID ##
    mtch <- match(CountData$ensmbl, gene.info$ENSEMBL)
    gene.info <- gene.info[mtch,]
    
    d <- DGEList(counts = CountData[ , -1], group=group$gender,
                 genes = gene.info)

    ## REMOVE GENES IF 20% OF SAMPLES HAVE 6 OR FEWER READS ##
    prop.gt.minread <- apply(d$counts, 1, function(x) mean(x > min.reads))
    keep <- which(prop.gt.minread >= min.samples)
    d <- d[keep,]
    
    d <- calcNormFactors(d, method="TMM")

    d$samples$ID <- rownames(d$samples)
   
    svs <- generate.svs(d, covs=c("group")) 
    covs <- data.frame(covs[,-grep("gender",names(covs))], svs)

    ## ADD SEX AND SVS TO SAMPLES ##
    d$samples <- merge(d$samples, covs, by.x="ID", by.y="case.id", 
                       all.x = T, all.y=F, sort=F)
  
    ## Return edgeR normalized/rescaled CPM (counts per million)
    lcpm <- cpm(d, log=TRUE)
    keep.ind <- which(rowSums(is.na(d$samples[, grep("SV|BirthYear", names(d$samples))])) == 0)
    lcpmNorm <- removeBatchEffect(lcpm,
                                  covariates=d$samples[,
                                      grep("SV|BirthYear", names(d$samples))])

    ## REMOVE X CHROMOSOME GENES ##
    xy.genes <- gene.info$ENSEMBL[gene.info$chr %in% c("X","Y")]    
    auto.ind <- which(!rownames(lcpmNorm) %in% xy.genes)
    lcpmNormAuto <- lcpmNorm[auto.ind,]       
    
    saveRDS(d, file = paste0(DataDir , "Expression/", project, ".RDS"))
            
    ## MAKE ARACNE FILES ##
    make.aracne.files(lcpmNormAuto, d, project, paste0(DataDir,"Expression/Aracne/"))

    ## MAKE WGCNA FILES ##
    make.wgcna.files(lcpmNormAuto, d, project, paste0(DataDir,"Expression/WGCNA/"))
}



Create_miRNA_Data <- function(project, min.samp.size=40){
        
    ROOT <- "/gpfs/data/stranger-lab/askol/TCGA/"
    MIRDIR <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_microRNA/"
    setwd(MIRDIR)

    map <- make.mi.map(project)

    ## DETERMINE NUMBER OF TUMOR SAMPLES##
    no.samples <- length(unique(map$cases.0.case_id[-grep("Normal", map$sample_type)]))
    if (no.samples < min.samp.size){
        print(paste0("Only ",no.samples, " found. Skip making miRNA files."))
        return(0)
    }
    mir.info <- get.mir.info()
    mir.info$mir <- gsub("-",".",mir.info$mir)
    mir.info$chr <- as.character(gsub("chr","",mir.info$chr))
    
    mi.file <- paste0(MIRDIR,project,"/",project,"_miR_count.txt")
    
    tmp <- get.mi.data(mi.file, map)
    mi <- tmp$CountData
    map <- tmp$map
    
    mi$ID <- gsub("-","\\.",mi$ID)
    rownames(mi) <- mi$ID
    
    ## REMOVE SAMPLES WITH MISSING VALUES OF THE COV ##
    ind.na <- which(is.na(map$cases.0.demographic.gender))
    ids.na <- map$cases.0.case_id[ind.na]

    print(paste0("Removing ",length(ids.na), " cases with no gender information"))
    
    if (length(ind.na) > 0){
        map <- map[-ind.na,]
        ind.na = which(gsub("A","", names(mi)) %in% ids.na)
        mi = mi[,-ind.na]
    }

    ## ### --------- ######
    map$case.id <- paste0("A",map$cases.0.case_id)
    
    group <- data.frame(case.id = names(mi)[-1])
    group <- merge(group, map[,c("case.id", "gender")],by.x=1,
                   by.y="case.id", sort=F)

    mir.info <- mir.info[mir.info$mir %in% mi$ID, ]
    mtch <- match(mi$ID, mir.info$mir)
    mir.info <- mir.info[mtch, ]
    ## ALL CHROMOSOMES TOGETHER ##
    d <- DGEList(counts = mi[,-1], group=group$gender,
                 genes=mir.info)
    d <- calcNormFactors(d, method="TMM")
    
    
    ## REMOVE SAMPLES WITH MISSING SEX ##
    rm.ind <- which(is.na(d$samples$group))
    if (length(rm.ind) > 0){
        d <- d[,-rm.ind]
    }
    
    ## REMOVE GENES IF > 20% OF SAMPLES HAVE 
    propGTE6 <-  apply(d$counts, 1, function(x) mean(x >= 6))
    keep <- which(propGTE6 >= .2)
    d <- d[keep,]
    d <- calcNormFactors(d, method="TMM")
    d$samples$SEX <- 1*(d$samples$group == "male")

    ## CREATE SURROGATE VARIABLES (VIA MERI'S SCRIPT ##
    svs <- generate.svs(d, covs="gender")
    
    ## hidden_factors <- hidden_factors[common,]
    mtch <- match(rownames(d$samples), rownames(svs))
    d$samples <- data.frame(d$samples, svs[mtch,])
    
    ## SPLIT MI INTO AUTOSOME AND XY ##
    xy.mirs <- get.xy.mirs()
    mi.xy <- mi[mi$ID %in% xy.mirs, -1]
    mi.auto <- mi[!mi$ID %in% xy.mirs, -1]
    
    ## d.auto IS TO BE USED FOR AUTOSOMAL GENES (NORMALIZATION WAS DONE WITHOUT
    ## SEX CHROMOSOMES ##
    ## create the DGEList object and calculate variance
    ind <- which(d$genes$chr %in% c(1:22))

    d.auto <- d[ind, keep.lib.sizes=FALSE]
    d.auto <- calcNormFactors(d.auto, method="TMM")
    d.auto <- estimateCommonDisp(d.auto,verbose=TRUE)
    d.auto <- estimateTagwiseDisp(d.auto)  

    xy.ind <- rownames(d$counts) %in% xy.mirs
    d.xy <- d[xy.ind, ]

    ## CONVERT COUNTS TO COUNTS PER MILLION AND LOG COUNT PER MILLION
    cpm.auto <- cpm(d.auto)
    lcpm.auto <- cpm(d.auto, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0
    cpm.xy <- cpm(d.xy)
    lcpm.xy <- cpm(d.xy, log=TRUE) ## NOTE: LOG=TRUE ADDS 0.25 TO EACH COUNT TO AVOID 0
    sex.col = rep("red", ncol(lcpm.auto))
    sex.col[d.auto$samples$group == "male"] = "blue"

##    pdf(paste0("plotMDS_",project,".pdf"))
##    plotMDS(lcpm.auto, labels="o", col=sex.col)
##    dev.off()
    
    ## DISTRIBUTION OF COUNT DATA ##
##    plot.file <- paste0(project, "_mir_count_dist.pdf")
##    pdf(file = plot.file)
##    library(RColorBrewer)
##    nsamples <- ncol(lcpm.auto)
##    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
##    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
##        rownames(qual_col_pals)))
##    col=sample(col_vector, nsamples, replace=TRUE)

##    plot(density(lcpm.auto[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
##         main="", xlab="")
##    title(main="Filtered data", xlab="Log-cpm")
##    abline(v=0, lty=3)
##    for (i in 2:nsamples){
##        den <- density(lcpm.auto[,i])
##        lines(den$x, den$y, col=col[i], lwd=2)
##    }
##    dev.off()

    ## ## ##
    if (0){
        ##  ## ##
        
        ## CALCULATE PEERS ##
        ## EXPECTING ROWS = SAMPLES, COLUMNS = GENES ##
        model <- PEER()
        PEER_setPhenoMean(model,as.matrix(t(lcpm.auto)))
        PEER_setCovariates(model, as.matrix(1*(group$gender=="male")))
        
        ## SET 10 CONFOUNDERS ##
        PEER_setNk(model,10)
        
        ## RUN INFERENCE ##
        PEER_update(model)
        
        factors = PEER_getX(model)
        weights = PEER_getW(model)
        precision = PEER_getAlpha(model)
        residuals = PEER_getResiduals(model)
        
        factors <- data.frame(factors[,-1]) ## remove sex since alread in d$samples
        names(factors) <- paste0("PEER",1:ncol(factors))
        ## ADD FIRST 10 PEERS TO D$SAMPLE
        d.auto$samples <- cbind(d.auto$samples, factors)
        d.xy$samples <- cbind(d.xy$samples, factors)
        
        ## TEST FOR ASSOCIATION BETWEEN PEER AND SEX ##
        dat <- data.frame(factors, sex=d.auto$sample$group)
        names(dat)[ncol(dat)] = "sex"
        lms <- lapply(grep("PEER", names(dat)),
                      function(x) as.vector(t(
                          summary(lm(dat[,x] ~ dat$sex))$coef[2,])))
        lms.mtx <- do.call(rbind, lms)
        
        ## KEEP ALL PEERS, EVEN IF THEY ARE SIGNIFICANTLY ASSOCIATED WITH SEX ##
        ## rm.peer <- which(lms.mtx[,4] < 0.01)
        ##if (length(rm.peer) >0){
        ##    d.auto$samples <- d.auto$samples[,-which(names(d.auto$samples)
        ##                                            %in% paste0("PEER",rm.peer))]
        ##    d.xy$samples <- d.xy$samples[,-which(names(d.xy$samples)
        ##                                           %in% paste0("PEER",rm.peer))]
        ##}
        
        group <- d.auto$samples$group
        group.col <- rep("blue", length(group))
        group.col[group=="female"] = "red"
        mds <- c()
        for (i in 1:5){
            out <- plotMDS(lcpm.auto, labels="o", dim.plot=(i-1)*2 + c(1,2), col=group.col, plot=FALSE)
            mds <- cbind(mds, out$x, out$y)
        }
        
        mds <- as.data.frame(mds)
        names(mds) <- paste0("MDS",c(1:10))
        d.auto$samples <- cbind(d.auto$samples, mds)
        d.xy$samples <- cbind(d.xy$samples, mds)
        
        mds.sign <- c()
        for (i in 1:10){
            out <- summary(lm(mds[,i] ~ group))
            mds.sign <- rbind(mds.sign, c(i, out$coefficients[2,4], out$r.squared))
            if (out$coefficients[2,4] <= 0.01){
                d.auto$samples$tmp <- mds[,i]
                d.xy$samples$tmp<- mds[,i]
                names(d.auto$samples)[grep("tmp",names(d.auto$samples))] <- paste0("MDS",i)
                names(d.xy$samples)[grep("tmp",names(d.xy$samples))] <- paste0("MDS",i)
            }
        }
        
    }
    
    ## SAVE DATA ##
    file <- paste0("Limma_Voom_Use/",project,"_mirna_autosome_DGEList.RDS")
    saveRDS(d.auto, file=file)
    file <- paste0("Limma_Voom_Use/",project,"_mirna_sex_DGEList.RDS")
    saveRDS(d.xy, file=file)
    file<- paste0("Limma_Voom_Use/",project, "_mirna.RDS")
    saveRDS(d, file=file)
    
    print(paste0("Wrote ", file))
    print(paste0("Wrote ", gsub("sex","autosome", file)))
    
    return(1)

}


Get_Gene_Info <- function(){
  
    genes <- get.gene.annotation()
    
    ## ADD INFORMATION ABOUT GENES THAT HAVE WEAK/STRONG INTERACTION WITH MIRNA ##
    mi.gene <- read.table("/gpfs/data/stranger-lab/askol/TCGA/hsa_MTI-4.txt",
                          as.is=T, header=T, sep="\t")
    ## REMOVE 
    mi.gene$miRNA <- gsub("-3p|-5p","",mi.gene$miRNA)
    mi.gene$miRNA <- gsub("miR", "mir", mi.gene$miRNA)

   
    ## ADD INFORMATION ABOUT LOCATION TO MI.GENE ##
    mir.info <- get.mir.info()
    mi.gene <- merge(mi.gene, mir.info, by.x = "miRNA", by.y="mir",
                     all.x=T, all.y=F)

    mir.info.mirs <- strsplit(split="-", mir.info$mir)
    mir.info.mirs <- sapply(mir.info.mirs, function(x){paste(x[1:3], collapse="-")})
    ind <- which((mir.info$mir %in% mi.gene$miRNA) == F)
    sum(mir.info$mir %in% mi.gene$miRNA)
    
    
    ## CREATE STRONG AND WEAK MIR INTERACTION INDICATORS ON A GENE/MIR BASIS
    weak.ind <- grep("Weak", mi.gene$Support.Type)
    mi.gene$mir.weak <- mi.gene$mir.strong <- 0
    mi.gene$mir.weak[weak.ind] <- 1
    mi.gene$mir.strong[-weak.ind] <- 1

    ## CREATE A LIST OF WEAK AND STRONG INTERACTION GENES
    ## NOTE: WEAK GENES IS ANY GENE THAT DOESN'T HAVE A STRONG INDICATOR
    ## 
    strong.genes <- unique(mi.gene$Target.Gene[-weak.ind])
    weak.genes <- setdiff(unique(mi.gene$Target.Gene[weak.ind]),
                          strong.genes)
    genes$intx <- NA
    genes$intx[genes$SYMBOL %in% strong.genes] <- "strong"
    genes$intx[genes$SYMBOL %in% weak.genes] <- "weak"
    genes$intx[genes$SYMBOL %in% c(strong.genes, weak.genes) == FALSE] <- "none"  

    ## CREATE STRONG AND WEAK GENE-MIR INTACTION INDICATORS ##
    ## STRONG INTERACTIN GENES HAVE AT LEAST ONE STRONG INTERACTION
    mi.gene$gene.weak <- mi.gene$gene.strong <- 0
    mi.gene$gene.weak[mi.gene$Target.Gene %in% weak.genes] <- 1
    mi.gene$gene.strong[mi.gene$Target.Gene %in% strong.genes] <- 1

    ## ADD NAMES OF STRONG AND WEAK INTERACTING MIRS FOR EACH PROTEIN MIRS,
    ## THEIR CHROMOSOMES AND START POSITIONS
    ind.weak <- which(mi.gene$mir.weak == 1)
    ind.strong <- which(mi.gene$mir.strong == 1)
    mirs.weak <- with(mi.gene[ind.weak,],
                      tapply(miRNA, Target.Gene, paste, collapse=";"))
    mirs.weak <- data.frame(SYMBOL = rownames(mirs.weak), mir.weak = as.character(unlist(mirs.weak)),
                            stringsAsFactors=FALSE)
    mirs.weak.chr <- with(mi.gene[ind.weak,],
                          tapply(chr, Target.Gene, paste, collapse=";"))
    mirs.weak.chr <- data.frame(SYMBOL = rownames(mirs.weak.chr),
                                mir.weak.chr = as.character(unlist(mirs.weak.chr)))
    mirs.weak.pos <- with(mi.gene[ind.weak,],
                          tapply(start, Target.Gene, paste, collapse=";"))
    mirs.weak.pos <- data.frame(SYMBOL = rownames(mirs.weak.pos),
                                mir.weak.pos = as.character(unlist(mirs.weak.pos)))
    
    mirs.strong <- with(mi.gene[ind.strong,],
                      tapply(miRNA, Target.Gene, paste, collapse=";"))
    mirs.strong <- data.frame(SYMBOL = rownames(mirs.strong), mir.strong = as.character(unlist(mirs.strong)),
                            stringsAsFactors=FALSE)
    mirs.strong.chr <- with(mi.gene[ind.strong,],
                          tapply(chr, Target.Gene, paste, collapse=";"))
    mirs.strong.chr <- data.frame(SYMBOL = rownames(mirs.strong.chr),
                                mir.strong.chr = as.character(unlist(mirs.strong.chr)))
    mirs.strong.pos <- with(mi.gene[ind.strong,],
                          tapply(start, Target.Gene, paste, collapse=";"))
    mirs.strong.pos <- data.frame(SYMBOL = rownames(mirs.strong.pos),
                                mir.strong.pos = as.character(unlist(mirs.strong.pos)))
    
    genes <- merge(genes, mirs.weak, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.weak.chr, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.weak.pos, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong.chr, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    genes <- merge(genes, mirs.strong.pos, by="SYMBOL", all.x = TRUE, all.y=FALSE)
    cols <- c("mir.weak.chr", "mir.weak.pos", "mir.strong.chr", "mir.strong.pos")
    for (col in cols){
        genes[,col] <- as.character(genes[,col])
    }     

    ## SEX.IND : IS GENE ON X/Y? SEX.MIR*:DOES GENE INTERACT WITH A MIRNA ON X/Y?
    genes$sex.ind <- genes$sex.mir.weak.ind <-  genes$sex.mir.strong.ind <- 0
    xy.ind <- which(genes$chr %in% c("X","Y"))
    mir.xy.weak.ind <- grep("X|Y", genes$mir.weak.chr)
    mir.xy.strong.ind <- grep("X|Y", genes$mir.strong.chr)
    
    genes$sex.ind[xy.ind] <- 1
    genes$sex.mir.weak.ind[mir.xy.weak.ind] <- 1
    genes$sex.mir.strong.ind[mir.xy.strong.ind] <- 1
    
    ## COLLAPSE GENES CHRS / POSITIONS TOGETHER FOR SAME GENE SYMBOL ##
    dupe.ind.all <- which(duplicated(genes$ENSEMBL))
    if (length(dupe.ind.all) > 0){
        dupes <- unique(genes$ENSEMBL[dupe.ind.all])
        col.names.local <- c("SYMBOL", "chr", "start","end", "intx",
                             "sex.ind", "sex.mir.weak.ind", "sex.mir.strong.ind")
        col.names.global <- c("mir.weak", "mir.weak.chr", "mir.weak.pos",
                              "mir.strong", "mir.strong.chr", "mir.strong.pos")

        for (dupe in dupes){
            dupe.ind <- which(genes$ENSEMBL == dupe)
            ## ARE ALL THE GENE SYMBOLS THE SAME? ##
            genes.differ <- length(unique(genes$SYMBOL[dupe.ind])) > 1
            for (col.name in col.names.local){
                if (col.name == "SYMBOL" & genes.differ==F) { next }
                col.ind <- which(names(genes) == col.name)            
                tmp <- genes[dupe.ind, col.ind]
                tmp <- paste(tmp, collapse=";")
                genes[dupe.ind, col.ind] <- tmp
            }
            
            if (genes.differ == F){ next }
            for (col.name in col.names.global){
                col.ind <- which(names(genes) == col.name)            
                tmp <- genes[dupe.ind, col.ind]
                tmp <- paste(as.character(tmp), collapse=")(")
                tmp <- paste0("(",tmp, ")")
                genes[dupe.ind, col.ind] <- tmp
            }        
        }
    
    
        genes <- genes[-dupe.ind.all,]
    }

    return(genes)
}

get.gene.annotation <- function(){

    gene.info <- grch38
    ens.ind <- grep("ensgene", names(gene.info))
    names(gene.info)[ens.ind] <- "ENSEMBL"

    sym.ind <- grep("symbol", names(gene.info))
    names(gene.info)[sym.ind] <- "SYMBOL"     
    
    ## REMOVE DUPLICATES DUE TO HAVING MULTIPLE ENTREZ NAMES BUT SAME ENSEMBL ID AND POSITION
    rm.ind <- which(duplicated(
        paste(gene.info$ENSEMBL, gene.info$chr, gene.info$SYMBOL,
              sep=".")))
    if (length(rm.ind) > 0 ){
        gene.info <- gene.info[-rm.ind,]
    }

    ## REMOVE GENES NOT ON 1-22 OR X,Y
    keep <- which(gene.info$chr %in% c(1:22,"X","Y"))
    gene.info <- gene.info[keep, ]
    return(gene.info)
}


make.mi.map <- function(project){

    METADIR <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_MetaData"
      
    meta.miRNA.file <- paste0(METADIR,"/",project,"_miRNA_meta.txt")
      
    mi.map <- read.table(file=meta.miRNA.file, header=T, as.is=T, sep="\t")
    names(mi.map) <- gsub("file_id", "file_id_mi", names(mi.map))
    names(mi.map) <- gsub("cases.0.samples.0.","", names(mi.map))

    if (!any(grep("gender", names(mi.map)))){
        mi.map$cases.0.demographic.gender = NA
    }
    
    mi.map <- gender.fill(mi.map)
    
    return(mi.map)
}

## formally make.m.map 
get.covs <- function(project, MetaDir){
    
    ## meta.mRNA.file <- paste0(METADIR,"/",project,"_mRNA_meta.txt")
    meta.mRNA.file <- paste0(MetaDir,"/",project,"_mRNA_Count_meta.txt")
    
    covs <- read.table(file=meta.mRNA.file, header=T, as.is=T, sep="\t")
    
    names(covs) <- gsub("file_name", "file_name_m", names(covs))
    covs$file_name_m = gsub("\\..*","",covs$file_name_m)
    
    names(covs) <- gsub("cases.0.samples.0.","", names(covs))

    birth.ind <- grep("birth", names(covs))
    names(covs)[birth.ind] <- "BirthYear"
    covs <- gender.fill(covs)
    
    return(covs)
}

gender.fill <- function(covs){

    ## LOOK FOR PROBLEMS WITH GENDER ASSIGNMENT IN MAP FILE ##
sex.tbl <- table(covs[,c("cases.0.case_id","cases.0.demographic.gender")])

    ## IF SEX IS MISSING, THEN THERE WILL BE A COLUMN WITH NO COLNAME THAT HAS
    ## A 1 WHERE NO MALE OR FEMALE VALUE WAS STATED. REMOVE THIS COLUMN SINCE THE
    ## MALE AND FEMALE COLUMNS WILL EACH HAVE 0
    ind <- which(colnames(sex.tbl) == "")
    if (length(ind) >0){
        sex.tbl = sex.tbl[, -ind]
    }
    sex.tbl <- 1*(sex.tbl > 0)
    ind <- rownames(sex.tbl)[(rowSums(sex.tbl) > 1)]
    if (length(ind) > 0){
        print("Case ids ",rownames(sex.tbl)[ind], "have more than one sex assignment")
        print("Removing these samples . . .")       
    }

    ## TAKES CARE OF CASE of MULTIPLE SAMPLES AND ONE IS MISSING
    ## A GENDER ASSIGNMENT. IT WILL BE ASSIGNED THE SAME AS THE OTHER SAMPLES WITH
    ## GENEDER ASSIGNMENTS
    sex.update <- data.frame(case.id = rownames(sex.tbl),
                             gender = NA)
    ind = sex.tbl%*%c(1,2)
    sex.update$gender[ind !=0] <- colnames(sex.tbl)[ind[ind!=0]]

    covs <- merge(covs, sex.update, by.x = "cases.0.case_id", by.y="case.id", all.x=T,
                 all.y=FALSE)

    return(covs)
}

make.map.files <- function(project){
    
    m.map <- make.mi.map.file(project)
    mi.map <- make.m.map.file(project)
    
    return(list(m.map = m.map, mi.map = mi.map))

}

if (0){

    dupe.ids <- unique(m.map$cases.0.samples.0.sample_id[duplicated( m.map$cases.0.samples.0.sample_id)])
    
    for (id in dupe.ids){
        print(m.map[m.map$cases.0.samples.0.sample_id == id,])
    }
}

get.data <- function(m.file, mi.file, m.map, mi.map, gene.info){

    ## GET MI_RNA DATA ##
    mi <- get.mi.data(mi.file, mi.map)
    
    ## GET MRNA DATA ##
    m <- get.m.data(m.file, m.map, gene.info)
    return(list(m = m, mi = mi))
}    
    
    
get.mi.data <- function(file, mi.map){

    head <- read.table(file = file, as.is=T, nrow=1, header=F)
    mi <- read.table(file = file, as.is=T, header=F, skip=1)
    names(mi) <- c(paste0("A",head))
    ind <- grep("miRNA_ID", mi[,1])
    if (length(ind)>0){
        mi <- mi[-ind,]
    }
    mi[,-1] <- sapply(mi[,-1], as.numeric)
    rownames(mi) <- mi[,1]
        
    ## CLEAN UP DATA  (REMOVE NORMAL SAMPLES; AVERAGE OR REMOVE DUPES) ## 
    mi <- straighten.up(mi, "", mi.map, data.type="mi")

    return(mi)
}

get.m.data <- function(countFile, covs, gene.info,
                       min.reads = 5, min.samples=0.02){

    head <- read.table(file = countFile, as.is=T, nrow=1)
    head <- gsub("\\.htseq","",head)
    head <- gsub("\\.counts\\.gz","",head)

    m <- read.table(file = countFile, as.is=T, header=F, skip=1)
    names(m) <- c(head[1], paste0("A",head[-1]))
    m[,1] <-  gsub("\\.[0-9]+", "", m[,1])
    rownames(m) <- m[,1]
    ## THERE ARE COLUMNS LABELS __NO_FEATURE, __AMBIGUOUS, ETC THAT MUST BE
    ## REMOVED
    m <- m[grep("ENSG",rownames(m)), ]
    m.and.cov <- straighten.up(m, covs, data.type="m", min.reads, min.samples)
     
    return(m.and.cov)
}
 
straighten.up <- function(data,  covs, data.type, min.reads, min.samples){

    ## TRIM EXTRA COLUMNS IF DATATYPE IS MIRNA ##
    covs$file_id_use = ""
    
    if (data.type == "mi"){
        names(data)[1] = "ID"
        covs$file_id_use <- covs$file_id_mi
    }else{
        names(data)[1] <- "ensmbl"
        covs$file_id_use <- covs$file_name_m
        
    }
    
    ## REMOVE NORMAL ##
    norm.ind <- grep("Normal", covs$sample_type)
    norm.id <- covs$file_id_use[norm.ind]

    if (length(norm.ind)>0){

        ## REMOVE FROM COVS FILE ##
        covs <- covs[-norm.ind,]
        norm.ind <- which(gsub("^A","",names(data)) %in% norm.id)
        data <- data[, -norm.ind]
        print(paste0("Removing ",length(norm.ind), " normal samples."))
    }

    dupe.inds <- which(duplicated(covs$sample_id))
    dupes <- unique(covs$sample_id[dupe.inds])
    if (length(dupes)>0){
        print(paste0("Found ",length(dupe.inds)," duplicates from "
                     ,length(dupes), " samples"))
    }
    
    ## TAKE CARE OF DUPLICATE SAMPLES DEPENDING ON THE AMOUNT OF MISSINGNESS ##
    for (samp.id in dupes){

        file.ids <- covs$file_id_use[covs$sample_id %in% samp.id]       
        ind <- which(gsub("^A","", names(data)) %in% file.ids)

        if (length(ind) == 0){
            print(paste0("No mir/mrna data for samp.id ", samp.id, " found"))
            next
        }

        prop.zeros <- colSums(data[,ind]>0, na.rm=TRUE)/(nrow(data))

        min.zeros = min(prop.zeros)
        dif.zeros = prop.zeros - min.zeros
        rmv.dif.zeros <- which(dif.zeros > 0.1)
        
        ## REMOVE SAMPLES IF DIFFERENCE IN THE NUMBER OF ZEROS > 10% ##
        ## REMOVES SAMPLES WITH THE MOST ZEROS ##
        if (length(rmv.dif.zeros) > 0){

            print(paste0("Remove one sample because differences in missingness ",
                  "was greater than 10%"))
            data <- data[,-ind[rmv.dif.zeros]]
            ind <- ind[-rmv.dif.zeros]
        }

        ## AVERAGE EXPRESSION IF NUMBER OF REPS TWO OR MORE ##
        if (length(ind) > 1){

            print("Averaging duplicate count data")
            data[, ind[1]] = rowMeans(data[,ind], na.rm=T)
            data <- data[,-ind[-1],]
        }       
    }

    ## UPDATE COVs FILE ##
    ids <- gsub("^A", "", names(data))
    if (data.type == "mi"){
        covs$file_id_use <- covs$file_id_mi
    }else{
        covs$file_id_use <- covs$file_name_m
    }

    ## REMOVE INFO FOR SAMPLES NOT IN EXPRESSION DATA ##
    covs <- covs[covs$file_id_use %in% ids,]
    
    pref.order <- c("Primary Tumor", "Additional - New Primary",
                    "Recurrent Tumor", "Metastatic", "Additional Metastatic",
                    "Primary Blood Derived Cancer - Peripheral Blood")
    
    ## FIND DUPLICATE CASE IDS ##
    dupe.ind <- which(duplicated(covs$cases.0.case_id))
    dupe.ids <- unique(covs$cases.0.case_id[dupe.ind])

    if (length(dupe.ids) > 0){

        print(paste0(length(dupe.ind), " duplicates by case, made up of ",
                     length(dupe.ids), " cases."))

        ## PICK PREFFED SAMPLE FOR EACH DUPLICATE CASE ##
        for (id in dupe.ids){
                
            inds <- grep(id, covs$cases.0.case_id)
            file.ids <- covs$file_id_use[inds]
            
            ind.use <- match(covs$sample_type[inds], pref.order)
            ind.use <- which(ind.use == min(ind.use))
            if (length(ind.use) > 1){

                print(paste("Randomly selecting a sample for case ",id))
                ind.use = max(ind.use)
            }

            file.id.rm <- covs$file_id_use[inds[-ind.use]]

            data <- data[,-which(gsub("^A","", names(data)) %in% file.id.rm)]
            covs <- covs[-inds[-ind.use], ]
        }
    }    
    
    ## CLEAN DATA TO EXCLUDE GENES WITH LOW COUNT INFORMATION ##
    
    if (data.type == "mi"){
        data <- clean.data(data, data.type=data.type, min.reads = "", min.samples="")
    }

    ## REPLACE FILE ID (CURRENT HEADER) WITH CASE ID ##
    ord <- match(gsub("A","",names(data)), covs$file_id_use)
    names(data)[is.na(ord)==F] <- paste0("A",covs$cases.0.case_id[ord[is.na(ord)==F]])

    ## CONSOLIDATE COVS TO HAVE ONLY A SINGLE ENTRY PER CASE.ID ##
    covs <- covs[covs$cases.0.case_id %in% gsub("A","",names(data)),]
    tbl <-  table(covs$cases.0.demographic.gender)
    print(paste("Number", paste0(names(tbl),"s : ", tbl), collapse="; "))

    ## REMOVE CASES WITH NO GENDER ##
    ind <- covs$cases.0.case_id[is.na(covs$cases.0.demographic.gender)]
    ind <- which(names(data) %in% ind)
    if (length(ind) > 0){
        print(paste0("Removing ", length(ind), " cases with no designated sex"))
        data <- data[,-ind]
    }
        
    return(list(CountData = data, covs = covs))
}
    
clean.data <- function(data, data.type = "m",
                       min.reads=5,
                       min.samples=0.02){

    num.cols <- which(sapply(data, is.numeric))
    rna.id.col <- which(names(data) %in% c("ID","ensmbl"))
    
    if (data.type == "mi"){
        cpm <- sweep(data[,-1], 2, colSums(data[,-1]), FUN="/")*10^6
        keep <- which( rowMeans(cpm <=.1, na.rm=T)  < .8)
    }

    ## GTEX CRITERIA IS EXCLUDE IF 20% OF SAMPLES HAVE 6 OR FEWER READS OR
    ## IF 20% OF SAMPLES HAVE TPM OF <= 0.10 {CHANGED TO 80% BECAUSE GOAL IS
    ## IS DIFFERENTIAL EXPRESSION}
    if (data.type == "m"){

        rpm <- data[, sapply(data, is.numeric)]
        ## GTEX V6 KEPT GENES IF THEY HAD READS/MILLION > 5 IN > 2% OF SAMPLES ##
        rpm.factor <- colSums(rpm, na.rm=T)/(10^6)
        rpm <- sweep(rpm, 2, rpm.factor, "/")

        keep <- rowSums(rpm >= min.reads, na.rm=T)/rowSums(is.na(rpm)==F)
        keep <- which(keep > min.samples)       
    }
    
    data <- data[keep , ]

    return(data)
}

concat.dupes <- function(gi){

    dupe.ind <- which(duplicated(gi$ensembl_gene_id))
    dupe.names <- unique(gi$ensembl_gene_id[dupe.ind])

    for (name in dupe.names){

        ind <- which(gi$ensembl_gene_id == name)
        gi$hgnc_symbol[ind[1]] <- paste(gi$hgnc_symbol[ind], collapse=";")
        gi <- gi[-ind[-1],]
    }

    return(gi)
}

process.mi.gene <- function(data){
    
    data$miRNA = gsub("-",".",data$miRNA)
    data$miRNA = gsub("miR","mir", data$miRNA)
    data$miRNA = gsub("\\.5p|\\.3p","", data$miRNA)
    data$miRNA.alt1 <- paste0(data$miRNA, ".1")
    data$miRNA.alt2 <- paste0(data$miRNA, ".2")
    data$miRNA.alt3 <- paste0(data$miRNA, ".3")
    data$miRNA.alt4 <- paste0(data$miRNA, ".4")

    return(data)
}

get.xy.mirs <- function(){

    ## GET MIR LOCATIONS ##
    mir.locs <- read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/hsa.gff3", skip = 13,
                           header=FALSE, sep="\t")
    mir.locs[,9] = gsub(".*Name=","",mir.locs[,9])
    mir.locs[,9] = gsub(";.*","", mir.locs[,9])
    mir.locs <- mir.locs[,c(9,1)]
    names(mir.locs) <- c("mir","chr")
    xy.mirs <- mir.locs$mir[grep("x|y|X|Y",mir.locs$chr)]
    xy.mirs <- unique(xy.mirs)
    xy.mirs <- gsub("-","\\.",xy.mirs)

    ## REMOVE 3P/5P FROM END OF MIRNA NAME ##
    xy.mirs <- gsub("miR","mir", xy.mirs)
    xy.mirs <- gsub("\\.3p$|\\.5p$","", xy.mirs)

    return(xy.mirs)
}

get.mir.info <- function(){

    ## GET MIR LOCATIONS FOR PRIMARY TRANSCRIPTS ##
    mir.locs <- read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/hsa.gff3", skip = 13,
                           header=FALSE, sep="\t")
    ## KEEP ONLY PRIMARY TRANSCRIPT ENTRIES ##
    keep.ind <- grep("primary",mir.locs[,3])
    mir.locs <- mir.locs[ keep.ind, c(9,1,4,5)]
    names(mir.locs) <- c("mir","chr","start","end")
    mir.locs$mir = gsub(".*Name=","",mir.locs$mir)
    mir.locs$mir = gsub(";.*","", mir.locs$mir)
    
    return(mir.locs)
}

get.mir.info.trim.name <- function(){

    mir.locs <- get.mir.info()

    ##    
    mirs <- strsplit(split="-", mir.locs$mir)
    mirs <- sapply(mirs, function(x){paste(x[1:3], collapse="-")})
    mir.locs$mir.long <- mir.locs$mir
    mir.locs$mir <- mirs

    chrs <- with(mir.locs, tapply(chr, mir, paste, collapse=","))
    chrs <- data.frame(mir = rownames(chrs), chrs=chrs)
    
    starts <- with(mir.locs, tapply(start, mir, paste, collapse=","))
    starts <- data.frame(mir = rownames(starts), starts=starts)

    ends <- with(mir.locs, tapply(end, mir, paste, collapse=","))
    ends <- data.frame(mir = rownames(ends), ends=ends)

    mir.longs <- with(mir.locs, tapply(mir.long, mir, paste, collapse=","))
    mir.longs <- data.frame(mir = rownames(mir.longs), ends=ends)
    
    mir.locs <- merge(mir.locs, chrs, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, starts, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, ends, by="mir", all.x=T, all.y=F)
    mir.locs <- merge(mir.locs, mir.longs, by="mir", all.x=T, all.y=F)
    
    ## REMOVE DUPES ##
    dupe.ind <- which(duplicated(mir.locs$mir))
    mir.locs <- mir.locs[-dupe.ind , -which(names(mir.locs) %in%
                                            c("chr","start","end"))]
    
    return(mir.locs)
}



get.mir.gene.combos <- function(mir.info, gene.info, dist.from=20000){

    ord <- order(mir.info$chr, mir.info$start)
    mir.info <- mir.info[ord,]

    ord <- order(gene.info$chr, gene.info$start)
    gene.info <- gene.info[ord,]

    left.idx <- 0
    right.idx <- 0

    gene.mir.mtx <- data.frame(matrix(0,nrow(mir.info), nrow(gene.info)))
    names(gene.mir.mtx) <- gene.info$ENSEMBL
    rownames(gene.mir.mtx) <- mir.info$mir
    mtx.count <- 1
    
    for (chr in 1:22){

        mir.inds <- which(mir.info$chr == chr)
        mir.poss <- mir.info$start[mir.inds]
        gene.inds <- which(gene.info$chr == chr)
        gene.poss <- gene.info$start[gene.inds]

        ind <- which(abs(mir.poss[1] - gene.poss) < dist.from)
        if (length(ind) == 0){ ind = 1 }
        left.idx <- min(ind)
        right.idx <- max(ind)
        if (length(ind) > 0){
            gene.mir.mtx[mtx.count , left.idx:right.idx] = 1
        }
        mir.count = 1
        
        while(mir.count < length(mir.poss)){

            mtx.count <- mtx.count + 1
            mir.count <- mir.count + 1
            mir.pos <- mir.poss[mir.count]
            
            ## FIND NEXT BEGINNING ##

            dist <- gene.poss[right.idx] - mir.pos
            incr = 1
            while(dist < dist.from & (right.idx + incr) <= length(gene.poss)){
                incr <- incr + 1 
                dist <- gene.poss[right.idx+incr] - mir.pos
            }
            if (incr > 1 ){ incr <- incr - 1}
            right.idx <- right.idx + incr

            dist <- gene.poss[left.idx] - mir.pos
            while(abs(dist) > dist.from & left.idx < right.idx){
                left.idx <- left.idx + 1
                dist <- gene.poss[left.idx] - mir.pos
            }           

            if (left.idx != right.idx){
                gene.mir.mtx[mtx.count, left.idx:right.idx] = 1
            }
            
        }

    }

    return(gene.mir.mtx)
}
    

invNorm <- function(vals){

    trans <- qnorm((rank(vals , na.last="keep")-0.5)/sum(!is.na(vals)))
    return(trans)
}


create.quant.norm.inv.norm <- function(data, xy.genes){

    d.auto.q <- list()
    d.xy.q <- list()

    ind.xy <- which(rownames(data$counts) %in% xy.genes)
    ind.auto <- which(!rownames(data$counts) %in% xy.genes)
    
    d.auto.q$counts <- q.norm.inv.norm(data$count[ind.auto,])
    rownames(d.auto.q$counts) <- rownames(data)
    colnames(d.auto.q$counts) <- colnames(data)
    
    d.auto.q$samples <- data$sample
    svs.auto <- get.svs(d.auto.q$count)
    
    num.svs.auto.quant <- num.sv(d.auto.quant, mod,
                                 method = "be", vfilter = NULL, B = 30, seed = NULL)
    svs.auto.quant <- sva(d.auto.quant, mod, mod0, n.sv=num.svs.auto.quant, method="irw")
    svs.quant <- as.data.frame(svs.auto.quant$sv)
    names(svs.quant) <- paste0("sv.",1:num.svs.auto.quant)
    
    all <- q.norm.inv.norm(data$count[ind.xy,])
        
    d.xy.quant <- d.all.quant[xy.ind,]

    mod <- model.matrix(~gender, data=d.auto$samples)
    mod0 <- model.matrix(~1, data=d.auto$samples) 
    num.svs.auto.quant <- num.sv(d.auto.quant, mod,
                           method = "be", vfilter = NULL, B = 30, seed = NULL)
    svs.auto.quant <- sva(d.auto.quant, mod, mod0, n.sv=num.svs.auto.quant, method="irw")
    svs.quant <- as.data.frame(svs.auto.quant$sv)
    names(svs.quant) <- paste0("sv.",1:num.svs.auto.quant)
    

    ## CALCULATE SVS IN XY ##
    num.svs.xy <- num.sv(lcpm.xy, mod,
method = "be", vfilter = NULL, B = 30, seed = NULL)
svs.xy <- sva(lcpm.xy, mod, mod0, n.sv=num.svs.auto, method="irw")
}

q.norn.inv.norm <- function(data){
    
    qn <- normalize.quantiles(as.matrix(data))
    inorm <-  t(apply(qn, 1, invNorm))
    return(inorm)
}

get.svs <- function(data, samples){

    mod <- model.matrix(~gender, data=samples)
    mod0 <- model.matrix(~1, data=samples)
    
    num.svs <- num.sv(data, mod,
                      method = "be", vfilter = NULL, B = 30, seed = NULL)
 

    svs <- sva(data, mod, mod0, n.sv=num.svs, method="irw")
    svs <- as.data.frame(svs$sv)
    names(svs) <- paste0("sv.",1:num.svs)

    return(svs)
}



## FROM MERI'S CODE IN: ##
## /gpfs/data/gtex-group/sex_biased_regulation_v8/sexDE_v8_final/meri/
## software/generate.svs.R    
generate.svs <- function(data, covs=c()) {
    
    common=rownames(covs)[rownames(covs)%in%colnames(data$samples)]
    cov.ind <- which(names(data$samples) %in% covs)
    
    ## Return edgeR normalized/rescaled CPM (counts per million)
    form=paste("~", paste(covs, collapse="+"))
    design <- model.matrix(eval(parse(text=form)), data=data$sample)
    count.log2.cpm <- voom(counts=data, design=design)
    
    ## Inverse-rank-normalize counts
    ## norm_counts <- inverse_quantile_normalization(count.log2.cpm$E)
    norm_counts <- scale(count.log2.cpm$E)
    
    ## Calculate smartSVs blocking for covs
    ## Determine the number of SVs
    cmd <- paste0('Y.r <- t(resid(lm(t(norm_counts) ',form,
                  ' , data=data$samples)))')
    eval(parse(text=cmd))
    
    ## Add one extra dimension to compensate potential loss of 1 degree of freedom
    ## in confounded scenarios (very important)
    n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
    mod <- model.matrix(eval(parse(text=form)),data=data$samples)
    
    ## Modify the default parameters: iteration numbers (B) and learning rate (alpha)
    sv.obj <- smartsva.cpp(norm_counts, mod, mod0=NULL, n.sv=n.sv, B=1000, alpha=1)
    
    allSv <- sv.obj$sv
    colnames(allSv) <- paste0("SV", 1:n.sv)
    rownames(allSv) <- colnames(data)

    return(allSv)
}

inverse_quantile_normalization <- function(gct) {
        gct = t(apply(gct, 1, rank, ties.method = "average"));
        gct = qnorm(gct / (ncol(gct)+1));
        return(gct)
    }


make.aracne.files <- function(lcpm, data, project, OutDir){
    
    file <-  paste0(OutDir, project, "_lcpm_auto_male.txt")
    male.ind <- data$sample$ID[which(data$samples$group=="male")]
    lcpm.male <- lcpm[, colnames(lcpm) %in% male.ind]
    lcpm.male <- cbind(rownames(lcpm.male), lcpm.male)
    colnames(lcpm.male)[1] <- "gene"
    write.table(file=file, lcpm.male, quote=F, row.names=F, col.names=T, sep="\t")

    file <-  paste0(OutDir, project,"_lcpm_auto_female.txt")
    lcpm.female <- lcpm[,!colnames(lcpm) %in% male.ind]
    lcpm.female <- cbind(rownames(lcpm.female), lcpm.female)
    colnames(lcpm.female)[1] <- "gene"
    write.table(file=file, lcpm.female, quote=F, row.names=F, col.names=T, sep="\t")
}

make.wgcna.files <- function(lcpm, data, project, OutDir){

    file <-  paste0(OutDir,project,"_lcpm_auto_male.txt")
    male.ind <- data$sample$ID[which(data$samples$group=="male")]
    lcpm.male <- lcpm[,male.ind]
    write.table(file=file, lcpm.male, quote=F, row.names=T, col.names=T, sep="\t")

    file <-  paste0(OutDir, project,"_lcpm_auto_female.txt")
    lcpm.female <- lcpm[,!colnames(lcpm) %in% male.ind]
    write.table(file=file, lcpm.female, quote=F, row.names=T, col.names=T, sep="\t")

}
