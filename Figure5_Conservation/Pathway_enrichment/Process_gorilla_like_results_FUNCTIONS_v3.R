library('pheatmap')
library(ggplot2) 
library(scales) 


#########################################################################################################
# my.data.name <- "MSIgDB_Hallmark_Datasets"
# my.mouse.sigs = my.hallmark.mouse
# my.colnames = c("Heart","Liver","Cereb", "Brain")

get_enrich_balloons_all_species <- function(my.data.name, my.mouse.sigs, my.colnames = c("Heart","Liver","Cereb", "Brain") ) {
    
  # get files from dataset
  my.enrich.sets.mouse <- list.files("Mouse_Data/FDR5percent", pattern = my.data.name)
  
  my.enrich.sets.others <- list.files("FDR5percent",pattern = my.data.name)
  my.enrich.sets.rat <- my.enrich.sets.others[grep("RAT",my.enrich.sets.others)]
  my.enrich.sets.human <- my.enrich.sets.others[grep("Human",my.enrich.sets.others)]
  my.enrich.sets.killi <- my.enrich.sets.others[grep("Killifish",my.enrich.sets.others)]
  
  # get file names and path
  my.files.sets.mouse <- paste("Mouse_Data/FDR5percent",my.enrich.sets.mouse, sep="/")
  my.files.sets.rat   <- paste("FDR5percent",my.enrich.sets.rat ,sep="/")
  my.files.sets.human <- paste("FDR5percent",my.enrich.sets.human ,sep="/")
  my.files.sets.killi <- paste("FDR5percent",my.enrich.sets.killi ,sep="/")
  
  
  # reorder files based on colnames
  my.columns.mouse <- c() 
  for (i in 1:length(my.colnames)) {
    my.columns.mouse <- c(my.columns.mouse,grep(paste("_",my.colnames[i],sep=""),my.enrich.sets.mouse))
  }
  
  my.columns.rat <- c() 
  for (i in 1:length(my.colnames)) {
    my.columns.rat <- c(my.columns.rat,grep(paste("_",my.colnames[i],sep=""),my.enrich.sets.rat))
  }
  
  my.columns.human <- c() 
  for (i in 1:length(my.colnames)) {
    my.columns.human <- c(my.columns.human,grep(paste("_",my.colnames[i],sep=""),my.enrich.sets.human))
  }
  
  my.columns.killi <- c() 
  for (i in 1:length(my.colnames)) {
    my.columns.killi <- c(my.columns.killi,grep(paste("_",my.colnames[i],sep=""),my.enrich.sets.killi))
  }
  
  
  my.enrichment.files <- c(my.files.sets.mouse[my.columns.mouse],
                           my.files.sets.rat[my.columns.rat]  ,
                           my.files.sets.human[my.columns.human],
                           my.files.sets.killi[my.columns.killi])
  
  my.samples <- c("Mouse_heart","Mouse_liver","Mouse_cerebellum",
                  "Rat_heart","Rat_liver","Rat_whole_brain",
                  "Human_heart","Human_liver","Human_cerebellum",
                  "Killifish_liver","Killifish_whole_brain"
  )
  
  # get data from significant FDR 0.05
  my.tissues.kegg <- vector(length=length(my.enrichment.files), mode="list")
  names(my.tissues.kegg) <- my.samples
  
  for ( i in 1:length(my.samples)) {
    my.file <- my.enrichment.files[i]
    my.tissues.kegg[[i]]  <- read.csv(my.file,sep="\t", header=T)
    
  }
  
  
  my.pathways <- rownames(my.mouse.sigs)
  
  ####
  # prepapre output data
  # p-val matrix
  my.matrix <- matrix(0,length(my.pathways),length(my.samples)) # default: -log10(1) pval == 0 no enrichment
  
  # Enrichment matrix
  my.matrix2 <- matrix(0,length(my.pathways),length(my.samples)) # initialize with Enrichment = 0 if no enrich
  
  # matrix with record of significance
  my.matrix3 <- matrix(0,length(my.pathways),length(my.samples)) # to get sigificant pathways
  
  colnames(my.matrix) <- my.samples
  colnames(my.matrix2) <- my.samples
  colnames(my.matrix3) <- my.samples
  
  rownames(my.matrix) <- my.pathways
  rownames(my.matrix2) <- my.pathways
  rownames(my.matrix3) <- my.pathways
  
  # collect data from files
  for (i in 1:length(my.pathways)) {
    #print(my.pathways[i])
    
    for (j in 1:length(my.samples)) { # tissues 
      
      my.id <- which(my.tissues.kegg[[j]]$Gene_Set %in% my.pathways[i])
      
      if(length(my.id) == 1) { # if was significant in this tissue (and not on both tail ends, which would be 2)
        
        my.matrix[i,j] <- -log10(my.tissues.kegg[[j]]$p.val[my.id]) # log(0) is undefined
        
        if (my.tissues.kegg[[j]]$Direction[my.id] == 'UP') {
          my.matrix2[i,j] <- my.tissues.kegg[[j]]$Enrichment[my.id]
        } else if (my.tissues.kegg[[j]]$Direction[my.id] == 'DOWN'){
          my.matrix2[i,j] <- - my.tissues.kegg[[j]]$Enrichment[my.id]
        }
        
        my.matrix3[i,j] <- 1
        
      }
      
    }
  }
  
  
  # get into data frame (all mouse significant are plotted)
  my.res.enrich <- data.frame(my.matrix2)
  my.pval.enrich <- data.frame(my.matrix)
  
  # sort by average change in original mouse analysis (transcriptome figure), stored in my.mouse.sigs
  my.average <- apply(my.mouse.sigs,1,mean)
  my.sorted <- sort(my.average,index.return=T,decreasing=T)
  my.res.enrich2 <- my.res.enrich[my.sorted$ix,]
  
  my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])
  
  my.txtname <- paste('./Stats_tables/',
                      paste(Sys.Date(),"Enrichment_table_All_species",my.data.name,"pathways_significant_in_Mouse_data.txt", sep="_"),
                      sep="")
  
  write.table(my.res.enrich2,file=my.txtname,sep="\t",quote=F)
  
  my.res.enrich2$Pathnames <- rownames(my.res.enrich2)
  
  # format for ggplot
  my.res.enrich3 <- cbind(my.res.enrich2[,c('Pathnames',my.samples[1])],rep(my.samples[1],dim(my.res.enrich2)[1]),my.pval.enrich2[,my.samples[1]])
  colnames(my.res.enrich3) <- c('Pathnames','aging_signed_enricment','condition','minusLog10Pval')
  for ( h in 2:length(my.samples)) {
    my.new <- cbind(my.res.enrich2[,c('Pathnames',my.samples[h])],rep(my.samples[h],dim(my.res.enrich2)[1]),my.pval.enrich2[,my.samples[h]])
    colnames(my.new) <- colnames(my.res.enrich3)
    my.res.enrich3 <- rbind(my.res.enrich3, 
                            my.new)
    
  }
  
  
  my.max <- max(my.res.enrich3$aging_signed_enricment)
  my.min <- min(my.res.enrich3$aging_signed_enricment)
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  # to preserve the wanted order
  my.res.enrich3$condition <- factor(my.res.enrich3$condition, levels = unique(my.res.enrich3$condition))
  my.res.enrich3$Pathnames <- factor(my.res.enrich3$Pathnames, levels = rev(unique(my.res.enrich3$Pathnames)))
  
  my.pdfname <- paste('./MYDATA/',
                      paste(Sys.Date(),"Enrichment_BALLOON_plot_All_species",my.data.name,"pathways_significant_in_Mouse_data.pdf", sep="_"),
                      sep="")
  
  pdf(my.pdfname, onefile=F, height = max(5, length(my.pathways)/3), width=15)
  my.plot <- ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=aging_signed_enricment,size=minusLog10Pval))+ theme_bw()+ geom_point(shape = 16) 
  my.plot <- my.plot + ggtitle("Aging dysregylated pathways") + labs(x = "Tissue/condition", y = "Gene Set")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
  print(my.plot)
  dev.off()  
  
}


