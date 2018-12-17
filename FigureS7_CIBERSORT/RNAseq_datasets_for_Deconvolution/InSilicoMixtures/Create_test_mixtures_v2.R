setwd('/Volumes/MyBook_5 1//RNAseq_datasets_for_Deconvolution/Process_for_deconvolution/InSilicoMixtures/')
options(stringsAsFactors=F)

library('preprocessCore')

# 2016-10-27
# try to create textbook mixtures to determine sensitivity of analysis
# talk with Ktaja
#load('2017-01-20_aggregated_counts_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.RData')

# 2017-01-20
# create textbook mixtures to determine sensitivity of analysis
# try #2

load('../2017-01-18/2017-01-18_aggregated_counts_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.RData')


################################################
############   1. PreProcess Data   ############
################################################
# remove sparse datasets
my.sparse <- which(apply(my.count.matrix.2 == 0, 2, sum) > 17000) # 9 samples

my.count.matrix.process <- data.frame(my.count.matrix.2[,-c(1:2,my.sparse)])
rownames(my.count.matrix.process) <- my.count.matrix.2$GeneName

my.meta.data.process <- my.meta.data.2[-my.sparse,]


###### Try quantile normalization before making mixtures
# do simple library size normalization to help scaling before adding
my.count.matrix.process.NORM2 <- data.frame(normalize.quantiles(as.matrix(my.count.matrix.process),copy=TRUE))
rownames(my.count.matrix.process.NORM2) <- rownames(my.count.matrix.process)
colnames(my.count.matrix.process.NORM2) <- colnames(my.count.matrix.process)

# visualize data spread
pdf("2017-01-20_boxplot_RNASeq_counts_PSEUDONORM_Quantile.pdf", width=25, height=6)
boxplot(my.count.matrix.process.NORM2+0.01, outline=F,
        log = 'y', las = 2, cex.axis = 0.5, col="tomato",
        ylab="Raw counts quantile normalization")
dev.off()
###### 
# quantile looks good


################################################
############   2. Create Mixtures   ############
################################################


### get column numbers for major cells types and macrophages
my.astro <- which(as.character(my.meta.data.process$Cell_type) %in% 'Astrocytes')
my.neurons <- which(as.character(my.meta.data.process$Cell_type) %in% 'Neurons')
my.mph <- which(as.character(my.meta.data.process$Cell_type) %in% 'Macrophages')
my.mgl <- which(as.character(my.meta.data.process$Cell_type) %in% 'Microglia')
my.hepa <- which(as.character(my.meta.data.process$Cell_type) %in% 'Hepatocytes')
my.Fib <- which(as.character(my.meta.data.process$Cell_type) %in% 'Dermal_Fibroblasts')
my.cardC <- which(as.character(my.meta.data.process$Cell_type) %in% 'Cardiomyocytes')
my.cardFib <- which(as.character(my.meta.data.process$Cell_type) %in% 'Cardiac_Fibroblasts')

# get a representative expression per cell type
my.astro.exp    <- apply(my.count.matrix.process.NORM2[,my.astro   ],1,mean)
my.neurons.exp  <- apply(my.count.matrix.process.NORM2[,my.neurons ],1,mean)
my.mph.exp      <- apply(my.count.matrix.process.NORM2[,my.mph     ],1,mean)
my.mgl.exp      <- apply(my.count.matrix.process.NORM2[,my.mgl     ],1,mean)
my.hepa.exp     <- apply(my.count.matrix.process.NORM2[,my.hepa    ],1,mean)
my.Fib.exp      <- apply(my.count.matrix.process.NORM2[,my.Fib    ],1,mean)
my.cardC.exp    <- apply(my.count.matrix.process.NORM2[,my.cardC   ],1,mean)
my.cardFib.exp  <- apply(my.count.matrix.process.NORM2[,my.cardFib ],1,mean)


################################################
get_wt_mean <- function(exp.matrix,my.weights){
  
  my.result <- rep(0,dim(exp.matrix)[1])
  
  for (i in 1: dim(exp.matrix)[1]) {
    
    my.result[i] <- weighted.mean(exp.matrix[i,], my.weights)
  }
  
  return(my.result)
}
################################################


#### simulate mixtures
# Brain
my.astro.neurons.50_50 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp), c(0.5,0.5))
my.astro.neurons.30_70 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp), c(0.3,0.7))
my.astro.neurons.70_30 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp), c(0.7,0.3))

my.astro.neurons.mgl.0.1 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.499,0.499,0.001))
my.astro.neurons.mgl.1 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.495,0.495,0.01))
my.astro.neurons.mgl.5 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.475,0.475,0.05))
my.astro.neurons.mgl.10 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.45,0.45,0.1))
my.astro.neurons.mgl.15 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.425,0.425,0.15))
my.astro.neurons.mgl.20 <- get_wt_mean(cbind(my.astro.exp,my.neurons.exp,my.mgl.exp), c(0.4,0.4,0.2))


# Heart
my.cardio_fibro.90_10 <- get_wt_mean(cbind(my.cardC.exp,my.cardFib.exp), c(0.9,0.1))
my.cardio_fibro.95_5 <- get_wt_mean(cbind(my.cardC.exp,my.cardFib.exp), c(0.95,0.05))
my.cardio_fibro.99_1 <- get_wt_mean(cbind(my.cardC.exp,my.cardFib.exp), c(0.99,0.01))

my.cardio.mph.0.1 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.999,0.001))
my.cardio.mph.1 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.99,0.01))
my.cardio.mph.5 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.95,0.05))
my.cardio.mph.10 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.9,0.1))
my.cardio.mph.15 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.85,0.15))
my.cardio.mph.20 <- get_wt_mean(cbind(my.cardC.exp,my.mph.exp), c(0.8,0.2))


# Liver
my.liver_fibro.90_10 <- get_wt_mean(cbind(my.hepa.exp,my.Fib.exp), c(0.9,0.1))
my.liver_fibro.95_5 <- get_wt_mean(cbind(my.hepa.exp,my.Fib.exp), c(0.95,0.05))
my.liver_fibro.99_1 <- get_wt_mean(cbind(my.hepa.exp,my.Fib.exp), c(0.99,0.01))

my.liver.mph.0.1 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.999,0.001))
my.liver.mph.1 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.99,0.01))
my.liver.mph.5 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.95,0.05))
my.liver.mph.10 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.9,0.1))
my.liver.mph.15 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.85,0.15))
my.liver.mph.20 <- get_wt_mean(cbind(my.hepa.exp,my.mph.exp), c(0.8,0.2))



### get a mixture matrix
my.count.mixtures <- data.frame('Neuron_Astrocytes.50_50'= my.astro.neurons.50_50, ## name is swapped, astrocytes was first in mix
                                'Neuron_Astrocytes.30_70'= my.astro.neurons.30_70, ## name is swapped, astrocytes was first in mix
                                'Neuron_Astrocytes.70_30'= my.astro.neurons.70_30, ## name is swapped, astrocytes was first in mix
                                'Neuron_Astrocytes.Microglia_0.1'= my.astro.neurons.mgl.0.1,
                                'Neuron_Astrocytes.Microglia_1'= my.astro.neurons.mgl.1,
                                'Neuron_Astrocytes.Microglia_5'= my.astro.neurons.mgl.5,
                                'Neuron_Astrocytes.Microglia_10'= my.astro.neurons.mgl.10,
                                'Neuron_Astrocytes.Microglia_15'= my.astro.neurons.mgl.15,
                                'Neuron_Astrocytes.Microglia_20'= my.astro.neurons.mgl.20,
                                'Cardiomyocytes_CardiacFibro.90_10'= my.cardio_fibro.90_10,
                                'Cardiomyocytes_CardiacFibro.95_5'= my.cardio_fibro.95_5,
                                'Cardiomyocytes_CardiacFibro.99_1'= my.cardio_fibro.99_1,
                                'Cardiomyocytes_Macrophages_0.1'= my.cardio.mph.0.1,
                                'Cardiomyocytes_Macrophages_1'= my.cardio.mph.1,
                                'Cardiomyocytes_Macrophages_5'= my.cardio.mph.5,
                                'Cardiomyocytes_Macrophages_10'= my.cardio.mph.10,
                                'Cardiomyocytes_Macrophages_15'= my.cardio.mph.15,
                                'Cardiomyocytes_Macrophages_20'= my.cardio.mph.20,
                                'Hepatocytes_Fibro.90_10'= my.liver_fibro.90_10,
                                'Hepatocytes_Fibro.95_5'= my.liver_fibro.95_5,
                                'Hepatocytes_Fibro.99_1'= my.liver_fibro.99_1,
                                'Hepatocytes_Macrophages_0.1'= my.liver.mph.0.1,
                                'Hepatocytes_Macrophages_1'= my.liver.mph.1,
                                'Hepatocytes_Macrophages_5'= my.liver.mph.5,
                                'Hepatocytes_Macrophages_10'= my.liver.mph.10,
                                'Hepatocytes_Macrophages_15'= my.liver.mph.15,
                                'Hepatocytes_Macrophages_20'= my.liver.mph.20
                                )

rownames(my.count.mixtures) <- my.count.matrix.2$GeneName

### see summary to figure out if any multiplication will be necessary before rounding
# because of the mutliplication (simulates higher depth)
summary(my.count.mixtures)

my.count.mixtures.clean <- round(10*my.count.mixtures)

### run VST normalization for CIBERSORT
library('DESeq2')

my.big.type <- c(rep("Brain",9), rep("Heart",9), rep("Liver",9))

# design matrix
dataDesign = data.frame( row.names = colnames( my.count.mixtures.clean ), my.big.type = my.big.type)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.count.mixtures.clean,
                              colData = dataDesign,
                              design = ~ my.big.type)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

pdf("2017-01-20_Dsipersion_plot_for_Mixtures.pdf", width=25, height=6)
plotDispEsts(dds)
dev.off()

vsd <- data.frame(getVarianceStabilizedData(dds))

vsd$GeneName <- my.count.matrix.2$GeneName

vsd <- vsd[,c(28,1:27)]

save(dds,vsd,file=paste(Sys.Date(),"DEseq_processed_object_inSilicoMixtures.RData", sep=""))
write.table(vsd,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution_withPseudoBulkPool_VST_normalized_inSilicoMixtures.txt", sep="_"),
            quote=F,row.names=F, sep="\t")
##########################################################################################################################################################
