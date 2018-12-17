options(stringsAsFactors = FALSE)
library('pheatmap')

##############################################################################################

######### RNA feature extraction
# my.fpkm.file: filename/path to FPKM file
# my.3,my.12,my.29: vectors of column ranges for each age
# my.tissue: name of tissue
# creates a RData object and returns the gene list to be used

get_RNA_features <- function(my.fpkm.file, my.tissue, my.3=1:3, my.12=4:6, my.29=7:9) {
  
  # correct the offeset of gene name colum
  my.3 = my.3 + 1
  my.12 = my.12 + 1
  my.29 = my.29 + 1
  
  my.rna <- read.csv(my.fpkm.file,sep="\t",header=T)
  my.genes <- my.rna[,1]
  
  my.rna.med <- data.frame(FPKM_3m = apply(my.rna[,my.3],1,median),
                           FPKM_12m = apply(my.rna[,my.12],1,median),
                           FPKM_29m = apply(my.rna[,my.29],1,median))
  rownames(my.rna.med) <- my.genes
  
  my.filename <- paste(Sys.Date(),my.tissue,"RNA_feature_object.RData",sep="_")
  save(my.rna.med,file=my.filename)
  
  return(my.genes)
}



######### global aging feature extract from level feat matrix
# 
get_DE_features <- function (my.gene.list,my.RNAseq.process,my.tissue,my.fdr = 0.1) {
  
  my.genes <- intersect(rownames(my.RNAseq.process[[1]]),my.gene.list)
  
  # initialize resulting feature matrix
  my.feat.mat <- data.frame(matrix(NA,length(my.genes),4))
  colnames(my.feat.mat) <- c("GeneName", "age_FC", "FDR" )
  my.feat.mat$GeneName <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    my.gene <- my.genes[i]
    my.de.ix <- which(rownames(my.RNAseq.process[[1]]) %in% my.gene)
    
    # get significance information
    my.feat.mat[i,2] <- my.RNAseq.process[[1]][my.de.ix,]$log2FoldChange # age FC from DEseq
    
    my.feat.mat[i,3] <- my.RNAseq.process[[1]][my.de.ix,]$padj # FDR for change
    
  }
  
  # get categorical result on FC
  my.feat.mat$age_change <- 'UNCLEAR'
  my.feat.mat$age_change[intersect(which(my.feat.mat$FDR < my.fdr),which(my.feat.mat$age_FC > 0))] <- "UP"
  my.feat.mat$age_change[intersect(which(my.feat.mat$FDR < my.fdr),which(my.feat.mat$age_FC < 0))] <- "DOWN"
  
  ##### restrict the constant class to the genes within half a sd of 0 (no change)
  my.feat.mat$age_change[ intersect(which(my.feat.mat$FDR > my.fdr),which(abs(my.feat.mat$age_FC) < 0.5 * sd(my.feat.mat$age_FC))) ] <- "CONSTANT"
  
  
  my.feat.mat$age_change <- as.factor(my.feat.mat$age_change) # correct data type
  
  my.filename <- paste(Sys.Date(),my.tissue,"FOR_FC_feature_object.RData",sep="_")
  
  save(my.feat.mat,file=my.filename)
  
}



######### H3K4me3 breadth data feature extraction
# my.tissue: sample type name
# my.genes: gene list to be extracted

get_breadth_features_slopes <- function(my.tissue, my.genes, my.breadth.mat.file,
                                        my.cols = c(10,13,16,19,22,25), my.ages) {
  
  # my.cols <- c(10,13,16,19,22,25) : default for all NON NPC sample (other is below)
  # number of bp of metapeak covered
  my.breadth.mat <- read.csv(my.breadth.mat.file,header=T,sep="\t")
  my.new.data <- my.breadth.mat[,my.cols]
  
  my.y.cols <- which(my.ages %in% min(my.ages))
  
  my.K4breadth.qt <- apply(my.new.data,2,get_quantiles)
  my.BD.cat.y <- identify_BDs(apply(my.K4breadth.qt[,my.y.cols],1,mean))
  
  #### 1. extract absolute breadth data
  my.K4breadth.max <- data.frame(matrix(0,length(my.genes),4))
  colnames(my.K4breadth.max) <- c("maxbreadth_y","BD_qt_y","breadth_qt_slope","BD_stat_y")
  rownames(my.K4breadth.max) <- my.genes

  for (i in 1:length(my.genes)) {
    
    my.id <- which(my.breadth.mat$Gene.Name %in% my.genes[i])
    
    if (length(my.id) == 0) {
      # no corresponding signal: breadth of 0
      my.K4breadth.max[i,1] <- 0
      my.K4breadth.max[i,2] <- 0
      my.K4breadth.max[i,3] <- 0
      my.K4breadth.max[i,4] <- "None"
        
    } else if (length(my.id) == 1) {
      # only one domain
      # get mean of replicates
      my.K4breadth.max[i,1] <- mean(as.numeric(my.new.data[my.id,my.y.cols]))
      
      mini.mat <- data.frame( age=my.ages, 
                              prom_int = as.numeric(my.K4breadth.qt[my.id,]) )
      my.reg <- lm(prom_int ~ age , data = mini.mat)
      
      # get slope of mark with age
      my.K4breadth.max[i,2] <- mean(as.numeric(my.K4breadth.qt[my.id,my.y.cols]))
      my.K4breadth.max[i,3] <- my.reg$coefficients[2] # breadth_qt_slope
      my.K4breadth.max[i,4] <- as.character(my.BD.cat.y[my.id]) # BD_stat_y
      
    } else {
      # more than one domain
      # select max breadth to report
      my.breadth.means <- apply(as.vector(my.new.data[my.id,]),1,mean)
      my.max.id <- which(my.breadth.means == max(my.breadth.means))
      
      if (length(my.max.id)) {
        # if more than one equal to max, just pick the first one
        my.max.id <- my.max.id[1]
      }
      
      my.K4breadth.max[i,1] <- mean(as.numeric(my.new.data[my.id[my.max.id],my.y.cols]))
      
      mini.mat <- data.frame( age=my.ages, 
                              prom_int = as.numeric(my.K4breadth.qt[my.id[my.max.id],]) )
      my.reg <- lm(prom_int ~ age , data = mini.mat)
      
      # get slope of mark with age
      my.K4breadth.max[i,2] <- mean(as.numeric(my.K4breadth.qt[my.id[my.max.id],my.y.cols]))
      my.K4breadth.max[i,3] <- my.reg$coefficients[2] # breadth_qt_slope
      my.K4breadth.max[i,4] <- as.character(my.BD.cat.y[my.id[my.max.id]]) # BD_stat_y
      
    }
    
  }
  
  my.K4breadth.max$GeneName <- my.genes
  my.K4breadth.max$BD_stat_y <- factor(my.K4breadth.max$BD_stat_y)
  
  my.filename1 <- paste(Sys.Date(),my.tissue,"K4breadth_aging_with_ALL_features_object.RData",sep="_")
  save(my.K4breadth.max,file=my.filename1)  
  
}


# breadth helper functions:

get_quantiles <- function(x) {
  
  my.null <- which(x %in% 0)
  my.ecdf <- ecdf(x[-my.null]) # remove null domains (clogging up ecdf)
  return(my.ecdf(x))
}


identify_BDs <- function(x) {
  my.new.x <- rep("Sharp",length(x))
  
  my.bds <- (x > 0.95)
  my.new.x[my.bds] <- "Broad"
  
  my.null <- (x == 0)
  my.new.x[my.null] <- "None"
  
  return(factor(my.new.x))
}



######### global aging feature extract from level feat matrix
# modify to get median RNAseq FPKM across ages
# modify to get slopes for each type of data instead of values
# get values for y
get_age_slopes <- function (my.mark, my.tissue, my.file, my.genes.in, my.ages) {
    
  my.hist.data <- read.csv(my.file,header=F, sep="\t")
  
  my.genes <- intersect(unique(my.hist.data$V4), my.genes.in) # gene name
  
  my.young.age.idx <- which(my.ages %in% min(my.ages))
  
  my.feat.mat <- as.data.frame(matrix(0,length(my.genes),3))
  colnames(my.feat.mat) <- c("GeneName",paste(my.mark,"_averageIntensity_",min(my.ages),"m",sep=""),paste(my.mark,"age_slope",sep="_"))
  my.feat.mat$GeneName <- my.genes
  
  # get the mean intensity in young, and get the highest for each gene
  my.mean.ints.yg <- apply(my.hist.data[,-c(1:4)][,my.young.age.idx],1,mean)
  
  for ( i in 1:length(my.genes)) {
    
    my.gene.ix <- which (my.hist.data$V4 %in% my.genes[i])
    
    my.max.id <- which(my.mean.ints.yg[my.gene.ix] %in% max(my.mean.ints.yg[my.gene.ix]))
    
    mini.mat <- data.frame(age=my.ages, 
                           prom_int = as.numeric(my.hist.data[my.gene.ix[my.max.id[1]],-c(1:4)])[1:length(my.ages)]
    )
    my.reg <- lm(prom_int ~ age , data = mini.mat)
    
    # get average
    my.feat.mat[i,2] <- my.mean.ints.yg[my.gene.ix[my.max.id[1]]] # age_intensity
    
    # get slope of mark with age
    my.feat.mat[i,3] <- my.reg$coefficients[2] # age_slope
    
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,my.mark,"chromatin_aging_feature_object.RData",sep="_")
  save(my.feat.mat,file=my.filename)
  
  #return(my.feat.mat)
  
}


######### Changed Nucleosome feature extraction
# HOMER annotated DANPOS/DiNUP calls
# consensus nucleosomes + DiNUP called fold changes

get_max_abs_change <- function(my.vector){
  
  my.min <- min(my.vector) # most lost
  my.max <- max(my.vector) # most gained
  
  if (1/my.min > my.max) {
    return(my.min)
  } else {
    return(my.max)
  }  
}

get_Nuc_features <- function(my.tissue, my.genes, my.nuc.gain.file, my.nuc.lost.file, my.homer.dinup.file) {
  
  # read in annotated significant nucleosomes
  my.nuc.g <- read.csv(my.nuc.gain.file, header=T, sep="\t")
  my.nuc.l <- read.csv(my.nuc.lost.file, header=T, sep="\t")
  my.nuc.change <- read.csv(my.homer.dinup.file, header=T, sep="\t")
  
  # create output dataframe
  my.nuc.feats <- data.frame(GeneName=my.genes,
                             Increased_Nucs=rep(0,length(my.genes)),
                             Decreased_Nucs=rep(0,length(my.genes)),
                             Max_occupancy_log2FC_age=rep(0,length(my.genes)) # 0 is no change
                             
  )
  
  for (i in 1:length(my.genes)) {
    
    my.gene <- my.genes[i]
    
    # count number of differential nucleosomes attributed to eachh gene
    my.nuc.feats$Increased_Nucs[i] <- sum(my.nuc.g$Gene.Name %in% my.gene)
    my.nuc.feats$Decreased_Nucs[i] <- sum(my.nuc.l$Gene.Name %in% my.gene)
    my.tmp <- get_max_abs_change(my.nuc.change$Peak.Score[my.nuc.change$Gene.Name %in% my.gene])
    
    if ( (my.nuc.feats$Increased_Nucs[i] +  my.nuc.feats$Decreased_Nucs[i] > 0) && (my.tmp < 50000) ) { 
      # max in data is 11; if over 50000, infinity was returned
      # if no detected consistent change, keep a 0 log2FC
      # otehrwise, record log FC
      my.nuc.feats$Max_occupancy_log2FC_age[i] <- log2(my.tmp)
    }
    
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"Differential_Nucleosomes_3vs29m_feature_object.RData",sep="_")
  
  save(my.nuc.feats,file=my.filename)
  
}


######### SE-Lite data feature extraction (categorical status in young)
# my.tissue: sample type name
# my.genes: gene list to be extracted
# my.se.3m1.file, my.se.3m2.file, my.se.12m1.file, my.se.12m2.file, my.se.29m1.file, my.se.29m2.file: SE call files

get_SELite_stat <- function(my.tissue, my.genes, my.se.3m1.file, my.se.3m2.file) {
  
  my.se.3m1 <- read.csv(my.se.3m1.file,sep="\t",header=F) 
  my.se.3m2 <- read.csv(my.se.3m2.file,sep="\t",header=F) 
  
  # WARNING: some genes have multiple detected enhancers, some have none !!! 
  my.se.list <- list(my.se.3m1,my.se.3m2)
  
  my.unique.se.genes <- intersect(unique(my.se.3m1$V7),my.genes) # keep only genes for which we have RNA information (all the se files have smae gene list)
  
  my.SE.tmp <- data.frame(matrix('None',length(my.genes),2))
  colnames(my.SE.tmp) <- c("SE_3m","SE_3m")
  rownames(my.SE.tmp) <- my.genes
  
  for (i in 1:length(my.unique.se.genes)) {
    
    my.current.gene <- my.unique.se.genes[i]
    my.row.ix <- which(my.genes %in% my.current.gene)
    
    for (j in 1:2) {
      
      my.id <- which(my.se.list[[j]]$V7 %in% my.current.gene)
      
      if (length(my.id) != 0) {
        # if one of the enhancers is super, super
        row.idx <- grep("Super", my.se.list[[j]]$V8[my.id])
        
        if (length(row.idx) >= 1) {
          my.SE.tmp[my.row.ix,j] <- "Super"
          
        } else {
          my.SE.tmp[my.row.ix,j] <- "Typical"
          
        }
        
      }
      
    }
    
  }
  
  # now resolve bio replicates per age
  my.SE <- data.frame("Gene_Name"=my.genes,
                      "SE_3m"=rep('None',length(my.genes)))
  
  # will make it so that only genes that are reliably super are super !!!
  # if one super one typical: typical
  # if they are different, they cannot be none (because all the samples are assessed on pn the enhancers)
  
  for ( i in 1:length(my.genes) ){
    
    if (my.SE.tmp[i,1] %in% my.SE.tmp[i,2])   {
      
      my.SE[i,2] <- my.SE.tmp[i,1] # if identical Super or identical typical, just take the first
      
    } else {
      my.SE[i,2] <- "Typical" # if not identical, mark as Typical
    }
        
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"_SE_categorical_feature_object.RData",sep="_")
  save(my.SE,file=my.filename)  
}

######### SE-Lite data feature extraction
# my.tissue: sample type name
# my.genes: gene list to be extracted
# my.se.mat.file: meta SE call files

get_SE_max_Score_feature <- function (my.tissue, my.genes, my.se.mat.file, my.ages) {
  
  my.score.se.mat <- read.csv(my.se.mat.file,sep="\t",header=T) 
  
  ### compute SE scores (see SE function files)
  # all samples have the same number of replicates for K27ac
  my.lib.size <- apply(my.score.se.mat[,8:19],2,sum)
  
  # apply lib size normalization
  my.libSnorm <- my.score.se.mat[,8:19]
  for (i in 1:12) {
    my.libSnorm[,i] <- 1e6 * my.libSnorm[,i]/my.lib.size[i]
  }
  
  # normalization of libs per bp
  my.intnorm <- data.frame(matrix(0,dim(my.score.se.mat)[1],6))
  my.enh.breadth <- my.score.se.mat$End - my.score.se.mat$Start
  my.intnorm$Gene_Name <- my.score.se.mat$Gene.Name
  
  for (i in 1:6) {
    my.intnorm[,i] <- my.libSnorm[,i]/my.enh.breadth
  }
  
  # now build score matrix
  my.SE.score <- data.frame(matrix('0',length(my.genes),2))
  colnames(my.SE.score) <- c("SE_score_3m","SE_aging_slope")
  rownames(my.SE.score) <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    my.idx <- which(my.intnorm$Gene_Name %in% my.genes[i])
    
    if (length(my.idx) == 0) {
      # no corresponding signal: SE score of 0
      my.SE.score[i,] <- 0
      
    } else if (length(my.idx) == 1) {
      my.SE.score[i,1] <- mean(as.numeric(my.intnorm[my.idx,1:2]))
      
      
      mini.mat <- data.frame( age=my.ages, 
                              prom_int = as.numeric(my.intnorm[my.idx,1:6]) )
      my.reg <- lm(prom_int ~ age , data = mini.mat)
      
      # get slope of mark with age
      my.SE.score[i,2] <- my.reg$coefficients[2] # breadth_qt_slope      
      
    } else if (length(my.idx) != 0) { # more than one
      my.score.means <- as.numeric(apply(my.intnorm[my.idx,1:6],1,mean))
      my.max.id <- which(my.score.means == max(my.score.means))
      
      if (length(my.max.id)) {
        # if more than one equal to max, jsut pick the first one
        my.max.id <- my.max.id[1]
        
      }
      
      my.SE.score[i,1] <- mean(as.numeric(my.intnorm[my.idx[my.max.id],1:2]))
     
            
      mini.mat <- data.frame( age=my.ages, 
                              prom_int = as.numeric(my.intnorm[my.idx[my.max.id],1:6]) )
      my.reg <- lm(prom_int ~ age , data = mini.mat)
      
      # get slope of mark with age
      my.SE.score[i,2] <- my.reg$coefficients[2] # age SE slope
      
    }
  
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"_SE_scores_feature_object.RData",sep="_")
  save(my.SE.score,file=my.filename)  
}


##########################
get_SE_dist <- function (my.se.mat.file, my.se.score.file,my.tissue) {
  
  my.se.mat <- read.csv(my.se.mat.file,sep="\t")
  my.se.score <- read.csv(my.se.score.file,sep="\t", header=F)
  colnames(my.se.mat)[1] <- "PeakID"
  
  my.merged <- merge(my.se.mat[,c(1,10)],my.se.score, by.x = "PeakID", by.y = 'V4')
  
  my.genes <- unique(my.merged$V9)
  
  # HOMER gets distance to middle of peak
  # edge distance is:
  #   dist + 1/2length if upstream, 
  #   dist - 1/2length if downstream, 
  
  my.res <- rep(0,length(my.genes))
  
  for ( i in 1:length(my.genes)) {
    
    my.idx <- which(my.merged$V9 %in%  my.genes[i] )
    
    my.lengths <- my.merged$V3[my.idx] - my.merged$V2[my.idx]
    my.tss.dists <- my.merged$Distance.to.TSS[my.idx]
    my.max.se <- which( my.merged$V5[my.idx] == max(my.merged$V5[my.idx])) # distance to strongest SE
    my.res[i] <- my.tss.dists[my.max.se[1]]
  }
  
  my.SE.dist.feat <- data.frame('GeneName'=my.genes,'SE_TSS_Dist'=my.res)
  
  my.filename <- paste(Sys.Date(),my.tissue,"SE_TSS_distance_feature_object.RData",sep="_")
  save(my.SE.dist.feat,file=my.filename)
  
  return (my.SE.dist.feat)
}



######### Pol2 peaks feature extraction
# HOMER annotated Pol2 peak calls (From LICR)
get_Pol2_numbers <- function(my.tissue, my.genes, my.Pol2.file) {
  
  # read in annotated significant nucleosomes 3vs29m DANPOS HOMER
  my.Pol2 <- read.csv(my.Pol2.file, header=T, sep="\t")
  
  # create output dataframe
  my.Pol2.feats <- data.frame(Gene_Name=my.genes,
                             Pol2_peaks=rep(0,length(my.genes)),
                             Pol2_abs_TSS_dist=rep(0,length(my.genes)),
                             Pol2_MACS2_max_score=rep(0,length(my.genes))
  )
  
  for (i in 1:length(my.genes)) {
    
    my.gene <- my.genes[i]
    
    # count number of peaks attributed to eachh gene
    my.Pol2.feats$Pol2_peaks[i] <- sum(my.Pol2$Gene.Name %in% my.gene)
    my.Pol2.feats$Pol2_abs_TSS_dist[i] <- min(min(abs(my.Pol2$Distance.to.TSS[my.Pol2$Gene.Name %in% my.gene])),500000) # will get 500000 if in to min
    my.Pol2.feats$Pol2_MACS2_max_score[i] <- max(my.Pol2$Peak.Score[my.Pol2$Gene.Name %in% my.gene], 0)
    
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"Pol2_young_feature_object.RData",sep="_")
  
  save(my.Pol2.feats,file=my.filename)
  
}

######### promoter data feature extraction
# my.prom.density.file: 
# my.genes:

get_prom_features <- function(my.prom.density.file, my.mark, my.column.names) {
  my.prom.all <- read.csv(my.prom.density.file,sep="\t",header=F) 
  # lots of gene with alternate TSSs. Will extract the one with most signal
  
  # gene with RNA info but not TSS info: Will put NAs there!
  # columns
  
  my.genes <- unique(my.prom.all$V4)
  
  my.prom.av <- data.frame(matrix(0,length(my.genes),length(my.column.names)))
  colnames(my.prom.av) <- my.column.names
  rownames(my.prom.av) <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    my.id <- which(my.prom.all$V4 %in% my.genes[i])
    
    if (length(my.id) == 0) {
      # no corresponding signal: NA
      my.prom.av[i,] <- NA
      
    } else if (length(my.id) == 1) {
      # only one TSS
      for (j in 1:length(my.column.names)){
        my.prom.av[i,j] <- mean(as.numeric(my.prom.all[my.id,j+4]))
      }
      
    } else {
      # more than one TSS
      # select most intense TSS to report ?
      # get the one with on average the most signal
      # print(i)
      my.tss.means <- apply(as.vector(my.prom.all[my.id,-c(1:4)]),1,mean)
      my.max.id <- which(my.tss.means == max(my.tss.means))
      
      if (length(my.max.id)) {
        # if more than one equal to max, jsut pick the first one
        my.max.id <- my.max.id[1]
        
      }
      
      for (j in 1:length(my.column.names)){
        my.prom.av[i,j] <- mean(as.numeric(my.prom.all[my.id[my.max.id],j+4]))
      }
    }
    
  }
  
  my.filename <- paste(Sys.Date(),my.mark,"prom_feature_object.RData",sep="_")
  save(my.prom.av,file=my.filename)  
}


######### Traveling Ratio data feature extraction
# my.prom.density.file: 
# my.genes:
get_TR_features <- function(my.density.file, my.genes, my.tissue) {
  my.tr.all <- read.csv(my.density.file,sep="\t",header=T) 
  # lots of gene with alternate TSSs. Will extract the max TR
  
  my.tr.av <- data.frame(matrix(0,length(my.genes),3))
  colnames(my.tr.av) <- c("Promoter","GeneBody","TR")
  rownames(my.tr.av) <- my.genes
  
  for (i in 1:length(my.genes)) {
    
    my.id <- which(my.tr.all$GeneName %in% my.genes[i])
    
    if (length(my.id) == 0) {
      # no corresponding signal: NA
      my.tr.av[i,] <- NA
      
    } else if (length(my.id) == 1) {
      # only one TSS
      my.tr.av[i,1] <- as.numeric(my.tr.all$Prom_density[my.id])
      my.tr.av[i,2] <- as.numeric(my.tr.all$GeneBody_density[my.id])
      my.tr.av[i,3] <- as.numeric(my.tr.all$Traveling_ratio[my.id])
      
    } else {
      # more than one TSS
      # select most intense TSS to report ?
      # get the max signal
      # print(i)
      my.tss.means <- as.vector(my.tr.all$Traveling_ratio[my.id])
      my.max.id <- my.id[which(my.tss.means == max(my.tss.means))]
      
      if (length(my.max.id)) {
        # if more than one equal to max, just pick the first one
        my.max.id <- my.max.id[1]
        
      }
      
      my.tr.av[i,1] <- as.numeric(my.tr.all$Prom_density[my.max.id])
      my.tr.av[i,2] <- as.numeric(my.tr.all$GeneBody_density[my.max.id])
      my.tr.av[i,3] <- as.numeric(my.tr.all$Traveling_ratio[my.max.id])
    }
    
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"_TravelingRatio_feature_object.RData",sep="_")
  save(my.tr.av,file=my.filename)  
  
}

######### H3K4me3/H3K27me3 bivalency data feature extraction
get_bivalent <- function(my.tissue, my.genes, my.biv.file) {
  
  my.biv <- read.csv(my.biv.file, header=T, sep="\t")
  
  # create output dataframe
  my.biv.feats <- data.frame(Gene_Name=my.genes,
                              Bivalent_status=rep('None',length(my.genes)))
  
  my.biv.genes <- my.genes %in% my.biv$Gene.Name 
  my.biv.feats[my.biv.genes,2] <- "Bivalent"
  my.biv.feats$Bivalent_status <- factor(my.biv.feats$Bivalent_status)

  my.filename <- paste(Sys.Date(),my.tissue,"Bivalency_young_feature_object.RData",sep="_")
  
  save(my.biv.feats,file=my.filename)
  
}


######### CTCF peak distance feature extraction
# distance of closest CTCF peak
get_CTCF_feats <- function (my.tissue, my.genes.tissue, my.ctcf.file) {
  
  my.ctcf.mat <- read.csv(my.ctcf.file,sep="\t")
  
  # create output dataframe
  my.ctcf.feats <- data.frame(Gene_Name=my.genes.tissue,
                              CTCF_peaks=rep(0,length(my.genes.tissue)),
                              CTCF_abs_TSS_dist=rep(0,length(my.genes.tissue)),
                              CTCF_MACS2_max_score=rep(0,length(my.genes.tissue))
  )
  
  
  for (i in 1:length(my.genes.tissue)) {
    
    my.gene <- my.genes.tissue[i]
    
    # count number of peaks attributed to eachh gene
    my.ctcf.feats$CTCF_peaks[i] <- sum(my.ctcf.mat$Gene.Name %in% my.gene)
    my.ctcf.feats$CTCF_abs_TSS_dist[i] <- min(min(abs(my.ctcf.mat$Distance.to.TSS[my.ctcf.mat$Gene.Name %in% my.gene])),500000) # will get 500000 if in to min
    my.ctcf.feats$CTCF_MACS2_max_score[i] <- max(my.ctcf.mat$Peak.Score[my.ctcf.mat$Gene.Name %in% my.gene], 0)
    
  }
  
  my.filename <- paste(Sys.Date(),my.tissue,"CTCF_young_feature_object.RData",sep="_")
  
  save(my.ctcf.feats,file=my.filename)

}


##############
# remove NAs
rm_nas <- function(my.features){
  
  my.na <- rep(F,dim(my.features)[1])
  
  for (i in 1:dim(my.features)[1]) {
    
    if(sum(is.na(my.features[i,]) > 0)) {
      my.na[i] <- T
    }
  }
  #sum(my.na)
  return(my.features[!my.na,])
}


######### global aging feature extract from level feat matrix
# modify to get median RNAseq FPKM across ages
# modify to get slopes for each type of data instead of values
# get values for 3m
get_feat_matsv4 <- function (my.tissue,
                             my.rna.med,
                             my.exp.feat.mat,
                             my.K4me3.feat.mat,
                             my.K27ac.feat.mat,
                             my.SE,
                             my.SE.score,
                             my.SE.dists,
                             my.K4breadth.max,
                             my.tr.av,
                             my.nuc.feats,
                             my.k4me1.data,
                             my.ctcf.data,
                             my.Pol2.data,
                             my.Pol2.feats,
                             my.ctcf.feats,
                             my.exon.data.final.v2,
                             my.ucsc.cpg.all,
                             my.homer.prom.CG.v2,
                             my.ChIP.feat.matrix,
                             my.LOF.feat.matrix) {
  
  my.tr.av$GeneName <- rownames(my.tr.av)
  my.rna.med$GeneName <- rownames(my.rna.med)
  my.SE.score$GeneName <- rownames(my.SE.score)
  
  my.merge.1 <- merge(my.exp.feat.mat,my.tr.av[,-c(1:2)],by="GeneName")
  my.merge.1b <- merge(my.merge.1,my.rna.med[,c(1,4)],by="GeneName")
  my.merge.1c <- merge(my.merge.1b,my.K4me3.feat.mat,by="GeneName")
  my.merge.1d <- merge(my.merge.1c,my.K27ac.feat.mat,by="GeneName")
  my.merge.1e <- merge(my.merge.1d,my.SE,,by.x="GeneName", by.y="Gene_Name")
  my.merge.1f <- merge(my.merge.1e,my.SE.score,by="GeneName")
  my.merge.1g <- merge(my.merge.1f,my.SE.dists,by="GeneName")
  my.merge.1h <- merge(my.merge.1g,my.K4breadth.max,by="GeneName")
    
  my.tissue.dat <- my.k4me1.data[,c("GeneName",my.tissue)]
  colnames(my.tissue.dat) <- c("GeneName","H3K4me1_prom")
  my.merge.2 <- merge(my.merge.1h,my.tissue.dat,by="GeneName")
  
  my.tissue.dat <- my.ctcf.data[,c("GeneName",my.tissue)]
  colnames(my.tissue.dat) <- c("GeneName","CTCF_prom")
  my.merge.3 <- merge(my.merge.2,my.tissue.dat,by="GeneName")
  
  my.tissue.dat <- my.Pol2.data[,c("GeneName",my.tissue)]
  colnames(my.tissue.dat) <- c("GeneName","Pol2_prom")
  my.merge.4 <- merge(my.merge.3,my.tissue.dat,by="GeneName")
  
  my.merge.4a <- merge(my.merge.4,my.Pol2.feats,by.x="GeneName", by.y="Gene_Name")
  my.merge.4b <- merge(my.merge.4a,my.nuc.feats,by="GeneName")
  
  colnames(my.homer.prom.CG.v2)[1:3] <- c("GeneName","CpG_promoter_percentage","CG_promoter_percentage")
  my.merge.5 <- merge(my.merge.4b,my.homer.prom.CG.v2[,1:3],by="GeneName")
  
  my.merge.5a <- merge(my.merge.5,my.ctcf.feats,by="GeneName",by.y="Gene_Name")
  my.merge.5b <- merge(my.merge.5a,my.exon.data.final.v2,by="GeneName",by.y="Associated.Gene.Name")
  my.merge.5c <- merge(my.merge.5b,my.ucsc.cpg.all,by="GeneName",by.y="Gene.Name")

  my.merge.6 <- merge(my.merge.5c,my.ChIP.feat.matrix,by="GeneName")
  
  my.genename <- grep("GeneName",colnames(my.LOF.feat.matrix))
  colnames(my.LOF.feat.matrix)[-my.genename] <- paste("LOF_",colnames(my.LOF.feat.matrix)[-my.genename],sep="")
  my.merge.7 <- merge(my.merge.6,my.LOF.feat.matrix,by="GeneName")
  
  my.merge.7b <- my.merge.7[(my.merge.7$age_change != 'UNCLEAR'),] # remove unclear entries to have a tight constant set
  my.merge.7b$age_change <- factor(my.merge.7b$age_change)

  my.merge.7b$SE_score_3m <- as.numeric(my.merge.7b$SE_score_3m)
  my.merge.7b$SE_aging_slope <- as.numeric(my.merge.7b$SE_aging_slope)
  my.merge.7b$CpG_islands <- as.numeric(my.merge.7b$CpG_islands)
  
  return(my.merge.7b)
}


get_additional_features <- function(my.tissue,my.feat.mat,my.K27me3.data,my.access.data,my.biv.feats) {
  
  my.tissue.dat <- my.K27me3.data[,c("GeneName",my.tissue)]
  colnames(my.tissue.dat) <- c("GeneName","H3K27me3_prom")
  
  my.merge.1 <- merge(my.feat.mat,my.tissue.dat,by="GeneName")
  
  my.tissue.dat <- my.access.data[,c("GeneName",my.tissue)]
  colnames(my.tissue.dat) <- c("GeneName","DNAseI_prom")
  
  my.merge.2 <- merge(my.merge.1,my.tissue.dat,by="GeneName")
  
  my.merge.3 <- merge(my.merge.2,my.biv.feats,by.x="GeneName",by.y="Gene_Name")
  
  return(my.merge.3)
  
}
