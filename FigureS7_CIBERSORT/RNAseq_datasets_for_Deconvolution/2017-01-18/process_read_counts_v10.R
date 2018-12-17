setwd('/Volumes/MyBook_5 1//RNAseq_datasets_for_Deconvolution/Process_for_deconvolution/2017-01-18//')
options(stringsAsFactors=F)

### All feature counts files have been provided for code checking

# file with dataset name and accession and notes
my.datasets <- read.csv('Samples_2017-01-03.txt',header=T, sep="\t") # 376 datasets

# files with reads counted over mm9
my.files <- list.files("../Feature_count_output/",pattern = "\\.txt$") # 547 datasets

# exclude non overlapping
my.real.files <- my.datasets$read.counts[my.datasets$read.counts %in% my.files] # 1563 datasets

# contruct ouput dataframe
my.first.file <- read.csv(paste("../Feature_count_output",my.real.files[1],sep="/"),skip=1,header=T,sep="\t",stringsAsFactors=F)
my.count.matrix <- data.frame(GeneName=my.first.file$Geneid,
                              Length= my.first.file$Length)
my.meta.data <- data.frame(Sample=rep('',length(my.real.files)),
                           Cell_type=rep('',length(my.real.files)),
                           NOTES= rep('',length(my.real.files)),
                           Type= rep('',length(my.real.files))
                           )
rownames(my.count.matrix) <- my.first.file$Geneid


# loop over files and fill matrix
for (i in 1:length(my.real.files)) {
  
  my.current.file <- read.csv(paste("../Feature_count_output",my.real.files[i],sep="/"),skip=1,header=T,sep="\t",stringsAsFactors=F)
    
  my.split.name <- strsplit(colnames(my.current.file)[7], ".", fixed = TRUE)
  
  if (my.split.name[[1]][(length(my.split.name[[1]])-1)] %in% 'srt') {
    
    my.comp.name <- paste(my.split.name[[1]][(length(my.split.name[[1]])-3)], # length-3 give the third to last element (file name)
                          my.split.name[[1]][(length(my.split.name[[1]])-2)], 
                          my.split.name[[1]][(length(my.split.name[[1]])-1)],my.split.name[[1]][(length(my.split.name[[1]]))],
                          sep=".")
    my.root.name <- my.split.name[[1]][(length(my.split.name[[1]])-3)]
    
  } else {
    my.comp.name <- paste(my.split.name[[1]][(length(my.split.name[[1]])-2)], # length-2 give the third to last element (file name)
                          my.split.name[[1]][(length(my.split.name[[1]])-1)],my.split.name[[1]][(length(my.split.name[[1]]))],
                          sep=".")
    my.root.name <- my.split.name[[1]][(length(my.split.name[[1]])-2)]
  }
  

  my.line.num <- which(my.datasets$STAR.mapped.file %in% my.comp.name)
  
  my.cur.dataset <- my.datasets[my.line.num,]
  
  # add dataframe column
  my.count.matrix[,my.root.name] <- my.current.file[,7]
  
  my.meta.data[i,1] <- my.root.name
  my.meta.data[i,2] <- my.cur.dataset$Cell.type
  my.meta.data[i,3] <- my.cur.dataset$Notes
  my.meta.data[i,4] <- my.cur.dataset$Bulk_Single
    
}

save(my.count.matrix,my.meta.data,
     file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution.RData", sep="_"))

write.table(my.count.matrix,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution.txt", sep="_"),
           quote=F,row.names=F )
write.table(my.meta.data,file=paste(Sys.Date(),"metadata_matrix_for_Deconvolution.txt", sep="_"),
            quote=F,row.names=F )


# aggregate single cell data ???
my.singles <- which(my.meta.data$Type %in% 'Single')

my.count.matrix.sc <- my.count.matrix[,my.singles+2]

my.cell.types <- unique(as.character(my.meta.data$Cell_type[my.singles]))

# [1] "aNSCs"            "Astrocytes"       "Cardiomyocytes"   "Ependymal_cells"  "Macrophages"      "Microglia"        "Neuroblasts"      "Neurons"          "Oligodendrocytes"
# [10] "OPCs"    

##### create pools for pseudo bulk
my.cardio <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Cardiomyocytes')
my.count.matrix.cardio <- data.frame('Cardiomyocyte_pseudo_Bulk'= apply(my.count.matrix.sc[,my.cardio[1:19]],1,sum) )

my.meta.data.cardio <- data.frame(cbind(colnames(my.count.matrix.cardio),
                                        rep("Cardiomyocytes",length(colnames(my.count.matrix.cardio))),
                                        rep("Cardiomyocytes_PseudoBulk",length(colnames(my.count.matrix.cardio))),
                                        rep("PseudoBulk",length(colnames(my.count.matrix.cardio)))) )
colnames(my.meta.data.cardio) <- colnames(my.meta.data)


#
my.neurons <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Neurons')
#my.meta.data[my.singles[my.neurons],]; seggregated by neuron subtype; # 837 from steve quake
my.count.matrix.neuron <- data.frame('Neurons_5HT_pseudoBulk'= apply(my.count.matrix.sc[,my.neurons[1:22]],1,sum),
                                     'Neuron_olfactory_sensory_pseudoBulk'= apply(my.count.matrix.sc[,my.neurons[23:26]],1,sum),
                                     'Neurons_pyramidal_cortex_pseudoBulk'= apply(my.count.matrix.sc[,my.neurons[27:45]],1,sum),
                                     'Neurons_pyramidal_hippo_pseudoBulk'= apply(my.count.matrix.sc[,my.neurons[46:63]],1,sum),
                                     'Neurons_Striatum_pseudoBulk_1'= apply(my.count.matrix.sc[,my.neurons[64:213]],1,sum), # 150 cell bins
                                     'Neurons_Striatum_pseudoBulk_2'= apply(my.count.matrix.sc[,my.neurons[214:364]],1,sum),
                                     'Neurons_Striatum_pseudoBulk_3'= apply(my.count.matrix.sc[,my.neurons[365:515]],1,sum),
                                     'Neurons_Striatum_pseudoBulk_4'= apply(my.count.matrix.sc[,my.neurons[516:666]],1,sum),
                                     'Neurons_Striatum_pseudoBulk_5'= apply(my.count.matrix.sc[,my.neurons[667:817]],1,sum),
                                     'Neurons_Striatum_pseudoBulk_6'= apply(my.count.matrix.sc[,my.neurons[818:900]],1,sum)
                                      )


my.meta.data.neuron <- data.frame(cbind(colnames(my.count.matrix.neuron),
                                     rep("Neurons",length(colnames(my.count.matrix.neuron))),
                                     rep("Neurons_PseudoBulk",length(colnames(my.count.matrix.neuron))),
                                     rep("PseudoBulk",length(colnames(my.count.matrix.neuron)))) )
colnames(my.meta.data.neuron) <- colnames(my.meta.data)

# other from Steve quake
my.aNSCs <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'aNSCs')
#my.meta.data[my.singles[my.aNSCs],] # 7 from steve quake
my.count.matrix.aNSCs <- data.frame('aNSCs_Striatum_pseudoBulk_1'= apply(my.count.matrix.sc[,my.aNSCs[1:4]],1,sum),
                                     'aNSCs_Striatum_pseudoBulk_2'= apply(my.count.matrix.sc[,my.aNSCs[4:7]],1,sum)
)

my.meta.data.aNSCs <- data.frame(cbind(colnames(my.count.matrix.aNSCs),
                                        rep("aNSCs",length(colnames(my.count.matrix.aNSCs))),
                                        rep("aNSCs_PseudoBulk",length(colnames(my.count.matrix.aNSCs))),
                                        rep("PseudoBulk",length(colnames(my.count.matrix.aNSCs)))) )
colnames(my.meta.data.aNSCs) <- colnames(my.meta.data)

apply(my.count.matrix.aNSCs == 0, 2, sum)

#
my.astro <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Astrocytes')
#my.meta.data[my.singles[my.astro],] # 109 from steve quake and other signle cell dataset
my.count.matrix.astro <- data.frame('Astrocytes_ABISOLID_pseudoBulk'= apply(my.count.matrix.sc[,my.astro[1:2]],1,sum),
                                    'Astrocytes_Striatum_pseudoBulk_1'= apply(my.count.matrix.sc[,my.astro[3:55]],1,sum),
                                    'Astrocytes_Striatum_pseudoBulk_2'= apply(my.count.matrix.sc[,my.astro[56:109]],1,sum)
)

my.meta.data.astro <- data.frame(cbind(colnames(my.count.matrix.astro),
                                       rep("Astrocytes",length(colnames(my.count.matrix.astro))),
                                       rep("Astrocytes_PseudoBulk",length(colnames(my.count.matrix.astro))),
                                       rep("PseudoBulk",length(colnames(my.count.matrix.astro)))) )
colnames(my.meta.data.astro) <- colnames(my.meta.data)

#
my.epend <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Ependymal_cells')
#my.meta.data[my.singles[my.epend],] # 47 from steve quake and other signle cell dataset
my.count.matrix.epend <- data.frame('Ependymal_ABISOLID_pseudoBulk_1'= apply(my.count.matrix.sc[,my.epend[1:4]],1,sum),
                                    'Ependymal_ABISOLID_pseudoBulk_2'= apply(my.count.matrix.sc[,my.epend[5:8]],1,sum),
                                    'Ependymal_Striatum_pseudoBulk_1'= apply(my.count.matrix.sc[,my.epend[9:26]],1,sum),
                                    'Ependymal_Striatum_pseudoBulk_2'= apply(my.count.matrix.sc[,my.epend[27:47]],1,sum)
)

my.meta.data.epend <- data.frame(cbind(colnames(my.count.matrix.epend),
                                       rep("Ependymal",length(colnames(my.count.matrix.epend))),
                                       rep("Ependymal_PseudoBulk",length(colnames(my.count.matrix.epend))),
                                       rep("PseudoBulk",length(colnames(my.count.matrix.epend)))) )
colnames(my.meta.data.epend) <- colnames(my.meta.data)


###my.mphi <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Macrophages') # Will not use these 

#
my.microglia <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Microglia')
#my.meta.data[my.singles[my.microglia],] # 48 from steve quake 
my.count.matrix.microglia <- data.frame('microglia_pseudoBulk_1'= apply(my.count.matrix.sc[,my.microglia[1:12]],1,sum),
                                    'microglia_pseudoBulk_2'= apply(my.count.matrix.sc[,my.microglia[13:24]],1,sum),
                                    'microglia_pseudoBulk_3'= apply(my.count.matrix.sc[,my.microglia[25:36]],1,sum),
                                    'microglia_pseudoBulk_4'= apply(my.count.matrix.sc[,my.microglia[37:48]],1,sum)
)

my.meta.data.microglia <- data.frame(cbind(colnames(my.count.matrix.microglia),
                                       rep("Microglia",length(colnames(my.count.matrix.microglia))),
                                       rep("Microglia_PseudoBulk",length(colnames(my.count.matrix.microglia))),
                                       rep("PseudoBulk",length(colnames(my.count.matrix.microglia)))) )
colnames(my.meta.data.microglia) <- colnames(my.meta.data)

#
##my.neuroblasts <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Neuroblasts')

#
my.oligos <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'Oligodendrocytes')
#my.meta.data[my.singles[my.oligos],] # 43 from steve quake 
my.count.matrix.oligos <- data.frame('Oligodendrocytes_pseudoBulk_1'= apply(my.count.matrix.sc[,my.oligos[1:11]],1,sum),
                                        'Oligodendrocytes_pseudoBulk_2'= apply(my.count.matrix.sc[,my.oligos[12:22]],1,sum),
                                        'Oligodendrocytes_pseudoBulk_3'= apply(my.count.matrix.sc[,my.oligos[23:33]],1,sum),
                                        'Oligodendrocytes_pseudoBulk_4'= apply(my.count.matrix.sc[,my.oligos[34:43]],1,sum)
)

my.meta.data.oligos <- data.frame(cbind(colnames(my.count.matrix.oligos),
                                           rep("Oligodendrocytes",length(colnames(my.count.matrix.oligos))),
                                           rep("Oligodendrocytes_PseudoBulk",length(colnames(my.count.matrix.oligos))),
                                           rep("PseudoBulk",length(colnames(my.count.matrix.oligos)))) )
colnames(my.meta.data.oligos) <- colnames(my.meta.data)

#my.opcs <- which(as.character(my.meta.data$Cell_type[my.singles]) %in% 'OPCs')



########## also sum up low coverage Hepatocytes bulk samples 
# (otherwise, not enough covered genes for CIBERSORT)
my.hepa <- which(as.character(my.meta.data$Cell_type) %in% 'Hepatocytes')
my.meta.data[my.hepa,]

my.count.matrix.hepa <- data.frame('Hepatocyte_pooled_Bulk1'= apply(my.count.matrix[,my.hepa[1:2]+2],1,sum),
                                  'Hepatocyte_pooled_Bulk2'= apply(my.count.matrix[,my.hepa[3:5]+2],1,sum),
                                  'SRR1286091Aligned' = my.count.matrix[,my.hepa[6]+2],
                                  'SRR2040945Aligned' = my.count.matrix[,my.hepa[7]+2],
                                  'SRR2040946Aligned' = my.count.matrix[,my.hepa[8]+2],                                  
                                  'SRR2040947Aligned' = my.count.matrix[,my.hepa[9]+2]
                                  )

## Add together the separate runs for the cardiomyocytes samples
# (otherwise, not enough covered genes for CIBERSORT)
my.cardio.bulk <- setdiff(which(as.character(my.meta.data$Cell_type) %in% 'Cardiomyocytes'),my.singles)
my.meta.data[my.cardio.bulk,]

my.count.matrix.cardiobulk <- data.frame('SRR1390711Aligned'= my.count.matrix[,my.cardio.bulk[1]+2],
                                   'SRR1390712Aligned'= my.count.matrix[,my.cardio.bulk[2]+2],
                                   'SRR1390713Aligned' = my.count.matrix[,my.cardio.bulk[3]+2],
                                   'Cardiomyocytes_added_runs_1' = apply(my.count.matrix[,my.cardio.bulk[4:5]+2],1,sum),
                                   'Cardiomyocytes_added_runs_2' = apply(my.count.matrix[,my.cardio.bulk[6:7]+2],1,sum),                        
                                   'Cardiomyocytes_added_runs_3' = apply(my.count.matrix[,my.cardio.bulk[8:9]+2],1,sum),
                                   'SRR964799Aligned' = my.count.matrix[,my.cardio.bulk[10]+2],
                                   'SRR964800Aligned' = my.count.matrix[,my.cardio.bulk[11]+2]
                                   )


### add together low coverage replicates B_cells
# (otherwise, not enough covered genes for CIBERSORT)
my.bcell.bulk <- which(as.character(my.meta.data$Cell_type) %in% 'B_cells')
my.meta.data[my.bcell.bulk,]



my.count.matrix.bcellbulk <- data.frame('SRR1536411Aligned'= my.count.matrix[,my.bcell.bulk[1]+2],
                                         'SRR1536412Aligned'= my.count.matrix[,my.bcell.bulk[2]+2],
                                         'SRR1976588Aligned' = my.count.matrix[,my.bcell.bulk[3]+2],
                                        'SRR1976593Aligned' = my.count.matrix[,my.bcell.bulk[4]+2],
                                         'B_cells_pseudoBulk_1' = apply(my.count.matrix[,my.bcell.bulk[5:6]+2],1,sum),
                                         'B_cells_pseudoBulk_2' = apply(my.count.matrix[,my.bcell.bulk[7:9]+2],1,sum),
                                        'ERR674979Aligned' = my.count.matrix[,my.bcell.bulk[10]+2]
                                          )

######
### add together fibroblasts runs
my.fibro.bulk <- grep("Fibroblasts",as.character(my.meta.data$Cell_type) )

my.meta.data[my.fibro.bulk,]


my.count.matrix.fibro.bulk <- data.frame('SRR1390714Aligned'= my.count.matrix[,my.fibro.bulk[10]+2],
                                        'SRR1390715Aligned'= my.count.matrix[,my.fibro.bulk[11]+2],
                                        'SRR1390716Aligned' = my.count.matrix[,my.fibro.bulk[12]+2],
                                        'SRR1015749Aligned' = my.count.matrix[,my.fibro.bulk[13]+2],
                                        'SRR1015750Aligned'= my.count.matrix[,my.fibro.bulk[14]+2],
                                        'SRR1015751Aligned'= my.count.matrix[,my.fibro.bulk[15]+2],
                                        'SRR964795Aligned' = my.count.matrix[,my.fibro.bulk[16]+2],
                                        'SRR964796Aligned' = my.count.matrix[,my.fibro.bulk[17]+2],

                                        '3m1_cohort1_dermal_fibroAligned' = my.count.matrix[,my.fibro.bulk[18]+2],
                                        '3m2_cohort1_dermal_fibroAligned'= my.count.matrix[,my.fibro.bulk[19]+2],
                                        '3m3_cohort4_dermal_fibroAligned'= my.count.matrix[,my.fibro.bulk[20]+2],
                                        '3m4_cohort1_dermal_fibroAligned' = my.count.matrix[,my.fibro.bulk[21]+2],
                                        '3m6_cohort4_dermal_fibroAligned' = my.count.matrix[,my.fibro.bulk[22]+2],
                                        '3m7_cohort4_dermal_fibroAligned' = my.count.matrix[,my.fibro.bulk[23]+2],
                                        
                                        'Bone_fibroblasts_pseudo_bulk' = apply(my.count.matrix[,my.fibro.bulk[1:9]+2],1,sum),
                                        'Dermal_fibroblasts_pseudo_bulk' = apply(my.count.matrix[,my.fibro.bulk[24:36]+2],1,sum),
                                        'Thymus_fibroblasts_pseudo_bulk' = apply(my.count.matrix[,my.fibro.bulk[37:44]+2],1,sum)
                                        )


# summed the experiment with several types of fibroblasts because of low coverage

##################
my.meta.data.hepa <- data.frame(cbind(colnames(my.count.matrix.hepa),
                                     rep("Hepatocytes",length(colnames(my.count.matrix.hepa))),
                                     c( rep("Hepatocyte_pooled_Bulk",2) ,  my.meta.data$NOTES[my.hepa[6:9]] ),
                                     c( rep("pooled_Bulk",2),  my.meta.data$Type[my.hepa[6:9]] ) ) )
colnames(my.meta.data.hepa) <- colnames(my.meta.data)


my.meta.data.cardiobulk <- data.frame(cbind(colnames(my.count.matrix.cardiobulk),
                                      rep("Cardiomyocytes",length(colnames(my.count.matrix.cardiobulk))),
                                      c( my.meta.data$NOTES[my.cardio.bulk[1:3]] , my.meta.data$NOTES[my.cardio.bulk[c(4,6,8)]],my.meta.data$NOTES[my.cardio.bulk[10:11]] ),
                                      c( my.meta.data$Type[my.cardio.bulk[1:3]] , my.meta.data$Type[my.cardio.bulk[c(4,6,8)]],my.meta.data$Type[my.cardio.bulk[10:11]] ) ) )
colnames(my.meta.data.cardiobulk) <- colnames(my.meta.data)


my.meta.data.bcell <- data.frame(cbind(colnames(my.count.matrix.bcellbulk),
                                            rep("B_cells",length(colnames(my.count.matrix.bcellbulk))),
                                            c( my.meta.data$NOTES[my.bcell.bulk[1:4]] , my.meta.data$NOTES[my.bcell.bulk[c(5,7)]],my.meta.data$NOTES[my.bcell.bulk[9]] ),
                                            c( my.meta.data$Type[my.bcell.bulk[1:4]] , my.meta.data$Type[my.bcell.bulk[c(5,7)]],my.meta.data$Type[my.bcell.bulk[9]] ) ) )
colnames(my.meta.data.bcell) <- colnames(my.meta.data)


my.meta.data.fibro <- data.frame(cbind(colnames(my.count.matrix.fibro.bulk),
                                      c(rep("Cardiac_Fibroblasts",8),rep("Dermal_Fibroblasts",6),"Bone_Fibroblasts","Dermal_Fibroblasts","Thymus_Fibroblasts"),
                                       c( my.meta.data$NOTES[my.fibro.bulk[10:23]] , my.meta.data$NOTES[my.fibro.bulk[c(1,24,37)]]),
                                       c( my.meta.data$Type[my.fibro.bulk[10:23]] , my.meta.data$Type[my.fibro.bulk[c(1,24,37)]] ) ) )
colnames(my.meta.data.fibro) <- colnames(my.meta.data)



my.count.matrix.2 <- cbind(my.count.matrix[,-c(my.singles+2,my.hepa+2, my.cardio.bulk+2,my.bcell.bulk+2,my.fibro.bulk+2)],
                           my.count.matrix.cardio,
                           my.count.matrix.cardiobulk,
                           my.count.matrix.neuron,
                           my.count.matrix.hepa,
                           my.count.matrix.bcellbulk,
                           my.count.matrix.fibro.bulk,
                           my.count.matrix.aNSCs,
                           my.count.matrix.astro,
                           my.count.matrix.epend,
                           my.count.matrix.microglia,
                           my.count.matrix.oligos)


my.meta.data.2 <- rbind(my.meta.data[-c(my.singles,my.hepa, my.cardio.bulk,my.bcell.bulk,my.fibro.bulk),],
                        my.meta.data.cardio,
                        my.meta.data.cardiobulk,
                        my.meta.data.neuron,
                        my.meta.data.hepa,
                        my.meta.data.bcell,
                        my.meta.data.fibro,
                        my.meta.data.aNSCs,
                        my.meta.data.astro,
                        my.meta.data.epend,
                        my.meta.data.microglia,
                        my.meta.data.oligos)

##### write to output file
save(my.count.matrix.2,my.meta.data.2,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.RData", sep="_"))
write.table(my.count.matrix.2,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.txt", sep="_"),
            quote=F,row.names=F )
write.table(my.meta.data.2,file=paste(Sys.Date(),"metadata_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.txt", sep="_"),
            quote=F,row.names=F )


##################################################################
######################   VST normalization  ######################
##################################################################
#load('2017-01-18_aggregated_counts_matrix_for_Deconvolution_withPseudoBulk_pooledBulk.RData')

library('DESeq2')

###### visualize DATA in terms of spread
pdf("2017-01-18_boxplot_RNASeq_counts.pdf", width=25, height=6)
boxplot(my.count.matrix.2[,-c(1:2)]+0.01, outline=F,
        log = 'y', las = 2, cex.axis = 0.5, col="tomato",
        ylab="Raw counts")
dev.off()


#### clean DATA to apply VST transformation
my.count.matrix.process <- data.frame(my.count.matrix.2[,-c(1:2)])
rownames(my.count.matrix.process) <- my.count.matrix.2$GeneName

# too many rows with zeros are creating issues for DESeq2
# remove samples with too many zeros?
my.sparse <- which(apply(my.count.matrix.process == 0, 2, sum) > 17000) # 7 samples
my.meta.data.2[my.sparse,]

# 13                SRR1873552Aligned  Adipocytes                               gonadal white adipose tissue, male        Bulk
# 343               SRR1634678Aligned Macrophages C57BL6/J, RNA-Seq, peritoneal_thioglycolate-elicited_macrophages        Bulk
# 241         Hepatocyte_pooled_Bulk1 Hepatocytes                                           Hepatocyte_pooled_Bulk pooled_Bulk
# 252            B_cells_pseudoBulk_2     B_cells                        C57BL/6, Ex vivo, Germinal_center_B_cells        Bulk
# 273  Astrocytes_ABISOLID_pseudoBulk  Astrocytes                                            Astrocytes_PseudoBulk  PseudoBulk
# 276 Ependymal_ABISOLID_pseudoBulk_1   Ependymal                                             Ependymal_PseudoBulk  PseudoBulk
# 277 Ependymal_ABISOLID_pseudoBulk_2   Ependymal                                             Ependymal_PseudoBulk  PseudoBulk


my.count.matrix.process.filter <- my.count.matrix.process[,-my.sparse]

pdf("2017-01-18_boxplot_RNASeq_counts_samples_filterSparse.pdf", width=25, height=6)
boxplot(my.count.matrix.process.filter+0.01, outline=F,
        log = 'y', las = 2, cex.axis = 0.5, col="tomato",
        ylab="Counts")
dev.off()

# remove null/low everywhere genes (see deseq2 vignette)
my.null <- which(apply(my.count.matrix.process.filter, 1, sum) <= 50) # 1605 genes

# Now pull out the spike in genes
spikes.idx <- grep("ERCC-", rownames(my.count.matrix.process.filter))
my.exclude <- union(my.null,spikes.idx) # 1691
 
my.filtered.matrix <- my.count.matrix.process.filter[-my.exclude,]

my.filtered.metadata <- my.meta.data.2[-my.sparse,]

# cell type extract
# remove the removed data
cell_type <- as.character( my.meta.data.2$Cell_type[-my.sparse] ) # length: 280

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), cell_type = cell_type)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData = dataDesign,
                              design = ~ cell_type)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

pdf("2017-01-18_Dsipersion_plot.pdf", width=25, height=6)
plotDispEsts(dds)
dev.off()

vsd <- getVarianceStabilizedData(dds)

save(dds,file=paste(Sys.Date(),"DEseq_processed_object.RData", sep=""))
save(vsd,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution_withPseudoBulkPool_VST_normalized.RData", sep=""))
write.table(vsd,file=paste(Sys.Date(),"aggregated_counts_matrix_for_Deconvolution_withPseudoBulkPool_VST_normalized.txt", sep="_"),
            quote=F,row.names=F )

### visualize normalized data
pdf("2017-01-18_boxplot_RNASeq_counts_VST_normalized.pdf", width=25, height=6)
boxplot(vsd+0.01, outline=F,
        log = 'y', las = 2, cex.axis = 0.5, col="tomato",
        ylab="VST normalized Counts")
dev.off()


#### MDS analysis
mds.result <- cmdscale(1-cor(vsd,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf("2017-01-18_MDS_RNAseq_DESeq_norm_together.pdf", width=25, height=25)
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2,col=NULL)
text(x, y, cell_type,cex=0.5)
dev.off()

write.table(my.filtered.metadata,file=paste(Sys.Date(),"Filtered_metadata_file.txt", sep="_"),
            quote=F,row.names=F )


library('pvclust')

bla <- vsd
colnames(bla) <- paste(cell_type,1:length(cell_type),sep="_")

my.pv <- pvclust(bla, method.hclust="average",
                 method.dist="correlation", use.cor="pairwise.complete.obs",
                 nboot=5)

pdf("2017-01-18_PVCLUST_vsd.pdf", width=25, height=10)
plot(my.pv)
dev.off()


#################################################
# extract samples for signature matrix validation
my.unique.ct <- unique(cell_type)

# initialize output structures
my.extracted.ix <- c()
my.unique <- c()
my.validation.matrix <- data.frame(GeneName=rownames(vsd))

for( i in 1:length(my.unique.ct)) {
  
  # get corresponding cells
  my.samples <- which(cell_type %in% my.unique.ct[i])
  
  # if not at least 3, skip
  # otherwise, extract one sample
  if (length(my.samples) >= 3) {
    my.valid <- sample(my.samples, 1)
    my.extracted.ix <- c(my.extracted.ix,my.valid)
    my.validation.matrix[,my.unique.ct[i]] <- 2^vsd[,my.valid]
  } else if (length(my.samples) < 2) {
    my.unique <- c(my.unique,my.samples)
    
  }
  
  
}

dim(my.validation.matrix)
#21010    21

write.table(my.validation.matrix,file=paste(Sys.Date(),"Test_MixtureFile_for_CIBERSORT.txt", sep="_"),
            quote=F,row.names=F, sep="\t")

my.input.signature <- data.frame(cbind(GeneName=rownames(vsd[,-c(my.extracted.ix,my.unique)]),2^vsd[,-c(my.extracted.ix,my.unique)]))
dim(my.input.signature)
# 21010   259

write.table(my.input.signature,file=paste(Sys.Date(),"Reference_Sample_File_for_CIBERSORT.txt", sep="_"),
            quote=F,row.names=F, sep="\t")


# Make the phenotype matrix
# Prepare your phenotype classes file. Prepare your phenotype classes file according to the CIBERSORT formatting requirements, 
# or download the example phenotype classes file used in this tutorial.

# From online manual: https://cibersort.stanford.edu/manual.php
# A value of "1" indicates membership of the reference sample to the class as defined in that row, 
# a value of "2" indicates the class that the sample will be compared against, 
# and a value of "0" indicates that the comparison will be ignored.

my.remaining.cell.types <- cell_type[-c(my.extracted.ix,my.unique)]

# initialize phenomatrix with the different value (2)
my.pheno.matrix <- matrix(2,length(my.unique.ct),length(my.remaining.cell.types) )

for( i in 1:length(my.unique.ct)) {
  my.cur.cell_type <- my.unique.ct[i]
  
  my.pheno.matrix[i,which(my.remaining.cell.types %in% my.cur.cell_type)] <- 1
  
}
  
my.pheno.matrix.2 <- data.frame(cbind(my.unique.ct, my.pheno.matrix) )

write.table(my.pheno.matrix.2,file=paste(Sys.Date(),"PhenotypeClassesFile_for_CIBERSORT.txt", sep="_"),
            quote=F,row.names=F , col.names =F , sep="\t")
