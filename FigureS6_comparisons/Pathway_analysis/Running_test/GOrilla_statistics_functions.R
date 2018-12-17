library('mHG')
options(stringsAsFactors=F)


#install.packages('mHG', type = "source")



# Description Runs a minimum-hypergeometric (mHG) test as described in: Eden, E. (2007). Discovering Motifs in Ranked Lists of DNA Sequences.

# # package functions
#     mHG.test
# # Arguments
# # mHG.test(lambdas, n_max = length(lambdas))
# # lambdas  \{0,1\}^N, sorted from top to bottom.
# # n_max  the algorithm will only consider the first n_max partitions.
# 
# # Value
# # statistic  The mHG statistic.
# # p.value	The p-value for the test.
# # parameters	N - total number of white and black balls; B - number of black balls ;n_max - Max partition considered by the algorithm.
# # n	The index for which the mHG was obtained 
# # b	Sum over 1 <= i <= n of lambdas[i]


####################################################################################
get_stat_sorted <- function(my.data.process, my.tissue, my.sorting.statistic) {
  my.data.process <- data.frame(my.data.process)
  my.data.process$GeneName <- rownames(my.data.process) # needs to be added for the RNA DESeq objects
  
  # keep only one instance for each gene
  my.data.max <- aggregate(my.data.process, by = list(my.data.process$GeneName), max) # for decreasing list
  my.data.min <- aggregate(my.data.process, by = list(my.data.process$GeneName), min) # for increasing list
  
  # sort based on my.sorting.statistic, decreasing to get upregulated processes
  my.data.max.sorted <- sort(my.data.max[,my.sorting.statistic], index.return = T, decreasing = T)
  
  # sort based on my.sorting.statistic, increasing to get doanregulated processes (negative)
  my.data.min.sorted <- sort(my.data.min[,my.sorting.statistic], index.return = T, decreasing = F)
  
  return(list(my.data.max$GeneName[my.data.max.sorted$ix],my.data.min$GeneName[my.data.min.sorted$ix]) )
}

# my.run <- 3
# my.tissue <- "Bochkis_Liver"
# deseq.res <- my.bochkis.RNAseq.process[[1]]

####################################################################################
run_pathway_enrich <- function(my.tissue, deseq.res, 
                               my.sorting.statistic = "stat", 
                               my.gmt.files = my.gmt.sets[my.run],
                               my.set.names = my.gmt.set.names[my.run]) {
  
  test.results <- as.data.frame(deseq.res)
  rownames(test.results) <- toupper(rownames(test.results))
  
  # get sorted lists
  my.sorted.lists   <- get_stat_sorted(test.results, my.tissue, my.sorting.statistic)
  my.sorted.upreg   <- my.sorted.lists[[1]]
  my.sorted.downreg <- my.sorted.lists[[2]]
  
  for (i in 1:length(my.gmt.files)) {
    kegg.results <- mHG.pathway.enrichment(gene.sets.file = my.gmt.files[i],
                                           test.results = test.results,
                                           my.sorted.upreg,
                                           my.sorted.downreg,
                                           statistic = my.sorting.statistic)
    
    my.sig.num <- sum(kegg.results$adj.p.val < 0.05)
    
    my.outfile <- paste(Sys.Date(),my.tissue,my.set.names[i],my.sig.num,"significant_pathway_enrichment_FDR0.05.txt",sep="_")
    my.outfile.2 <- paste(Sys.Date(),my.tissue,my.set.names[i],my.sig.num,"significant_pathway_enrichment_ALL.txt",sep="_")
    my.outfile.3 <- paste(Sys.Date(),my.tissue,my.set.names[i],my.sig.num,"object.RData",sep="_")
    
    ## Save the significant pathways in a table :
    write.table(kegg.results[kegg.results$adj.p.val < 0.05,],
                col.names = T,
                row.names = F,
                quote = FALSE,
                sep="\t",
                file = paste("./FDR5percent/",my.outfile,sep=""))
    
    write.table(kegg.results,
                col.names = T,
                row.names = F,
                quote = FALSE,
                sep="\t",
                file =  paste("./ALL_PATHWAYS/",my.outfile.2,sep=""))
    
    save(kegg.results, file = paste("./RData/",my.outfile.3,sep="") )
    
  }
  
}


####################################################################################
parse.mHG.test <- function(my.test){
  my.N <- my.test$parameters[1]
  my.B <- my.test$parameters[2]
  my.n <- my.test$n
  my.b <- my.test$b
  
  enrich <- (my.b/my.n) / (my.B/my.N)
  my.res <- c(enrich,my.n,my.b,my.N,my.b,my.test$p.value)
  names(my.res) <- c("Enrichment","n","b","N","B","p-value")
  return(my.res)
}


##########
# gene.sets.file = my.gmt.files[i]
# test.results = test.results
# my.sorted.upreg
# my.sorted.downreg
# statistic = my.sorting.statistic



mHG.pathway.enrichment <- function(gene.sets.file, test.results, my.sorted.upreg, my.sorted.downreg, statistic = "stat") {

  ############## Import gmt file and create list of gene sets ##############
  gene.sets <- read.table(gene.sets.file,
                          sep = "\n",
                          quote = "",
                          stringsAsFactors = FALSE)
  gene.sets.l <- apply(gene.sets, 1,
                       function(x){
                         y <- strsplit(x, "\t")
                         return(as.character(unlist(y)))
                       })
  names(gene.sets.l) <- sapply(gene.sets.l, function(x) x[1])
  ## Remove empty strings:
  gene.sets.l <- lapply(gene.sets.l,
                        function(x) {
                          ind.empty <- x == ""
                          return(x[!ind.empty])
                        })
  
  ## Remove genes in gene lists that were not tested:
  # genes in upreg and downreg are the same, just different order
  gene.sets.l.red <- lapply(gene.sets.l, function(x){
    intersect(x, my.sorted.upreg)
  })
  ## Remove gene lists with less than 10 genes (not biologically meaningful)
  ind.empty <- sapply(gene.sets.l.red, length) < 10
  gene.sets.l.red <- gene.sets.l.red[!ind.empty]
  ######################################################################
  
  ## Calculate average of squared statistics per gene set and return all genes that are more extreme than the average:
  my.results.df <- data.frame("Gene_Set"='',
                              "Direction"='',
                              "Enrichment"='',
                              "n"='',
                              "b"='',
                              "N"='',
                              "B"='',
                              "p-value"='',
                              "Driver_Genes"='')
  
  for ( i in 1:length(gene.sets.l.red)) {
    # Enrichment = (b/n) / (B/N) (cf GOrilla)
    
    # upregs testing
    my.lambdas.upreg <- my.sorted.upreg %in% gene.sets.l.red[[i]]
    my.upreg.test <- mHG.test(my.lambdas.upreg, n_max = round(length(my.lambdas.upreg)/2) ) # only look at top half of list, will significantly decrease runtime
    my.upreg.res <- parse.mHG.test(my.upreg.test)
    driver.genes.up <- paste(my.sorted.upreg[my.lambdas.upreg],collapse=",")
    
    # downregs testing
    my.lambdas.downreg <- my.sorted.downreg %in% gene.sets.l.red[[i]]
    my.dwnreg.test <- mHG.test(my.lambdas.downreg, n_max = round(length(my.lambdas.upreg)/2)) # only look at top half of list, will significantly decrease runtime
    my.dwnreg.res <- parse.mHG.test(my.dwnreg.test)
    driver.genes.dwn <- paste(my.sorted.downreg[my.lambdas.downreg],collapse=",")
    ### make enrichment negative for downregulated so that we can plot
    
    my.results.df <- rbind(my.results.df,
                           c(names(gene.sets.l.red)[i],"UP",my.upreg.res,driver.genes.up),
                           c(names(gene.sets.l.red)[i],"DOWN",my.dwnreg.res,driver.genes.dwn) )
    
  }
  
  my.results.df <- my.results.df[-1,] # remove empty first line
  my.results.df$adj.p.val <- p.adjust(my.results.df$p.value, method = "BH")
  
  ind.ord <- order(my.results.df$p.val)
  my.results.df <- my.results.df[ind.ord, ]
  return(my.results.df)
}
