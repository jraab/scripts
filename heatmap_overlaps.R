#script is designed to take a list of genes and get pair-wise overlap frequency for a desired cancer. 
#2013-09-14 jr
#this script takes an arbitrary list of genes, and a cancer type and spits out a heatmap of the overlapping mutation
#frequency. It should work for relatively large matrices ok. 
#this is a good visualizaiton if you have a way to cull through potentially interesting gnes
# plots will pop up at the end, or you need to add a call to ggsave to save it
#mutations are defined as GISTIC score of +/- 2 or somatic mutation called 
###############LIBRARIES########################
library(cgdsr)
library(plyr)
library(ggplot2) 
library(reshape2)
library(biomaRt)
setwd('~/Dropbox/github/labwork/mutations/')
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
studies <- getCancerStudies(mycgds)[,1:2] #if you need to know which cancers are possible
############IMPORTANT VARIABLES######################
#only alter things in this block - no changes below
genesDF = read.table('~/Dropbox/gene_by_complex.txt', header=T, sep='\t', stringsAsFactors = F)
genes <- list(genesDF$gene)
#genes <- c('ARID1A', 'ARID1B', 'ARID2', 'SMARCA2', 'SMARCA4' )
#genes <- c('ARID1A', 'ARID2', 'ARID1B', 'APC', 'SMARCA2', 'KDM6A', 'SMARCA4', 
#        'INO80', 'PTEN', 'KRAS', 'PIK3CA', 'EZH2', 'TP53')
## genes <- c('CDKN2A', 'ARID1A', 'B2M', 'STK11', 'CASP8', 
#            'MBD3', 'TP53', 'CDKN1B',  'BAP1', 
#            'PBRM1', 'ARHGAP35', 'MAP2K4', 'ARID2', 'RASA1',
#           'MAP3K1', 'CMTR2', 'ZNRF3', 'SOX9', 'ACVR1B', 
#           'SMAD2', 'KDM6A', 'CTCF', 'RB1', 'LARP4B', 
#           'PHLPP1', 'RBM10', 'APC', 'RNF43', 'FAT1', 'ARID1B', 
#           'SETD2', 'NF1', 'BMPR2', 'NF2', 'AMER1', 'VHL', 
#           'AMOT', 'CDH1', 'KMT2D', 'PIK3R1', 'PTEN', 'SMARCA4', 
#           'CUX1', 'ASXL2', 'NOTCH1', 'MGA', 'FBXW7', 'EP300', 
#           'SMAD4', 'CREBBP', 'NCOR1', 'SACS', 'ATM', 'PIK3CA')
#genes <- c('INO80', 'ACTL6A', 'SMARCA4','TP53','BRCA1','BRCA2','PTEN','RB1',
#           'APC','CDKN2A','VHL','NF1','NF2','MEN1','RUNX1','RUNX3','DCC','KRAS',
#           'MLH1','CASP8','ARID1A','ARID1B','ARID2','ARID4A','ARID3A','PIK3CA','ATR','ATM',
#           'RAD54B','TOPBP1','NBN', 'PARP1','HAT1','RAD51','FANCD2', 'SOX2')
#chk1 chk2 tfbr2 not showing up

# controls genes queried
cancer_type <- 'ucec_tcga'
#values <- c( '#99CCFF','#666666','#CCCCCC') #controls colors of final plots

########
#SHOULD NOT NEED TO CHANGE THINGS BELOW THIS LINE
#filter the gene list - running the rest of the scirpt with invalid genes will break. 
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genes = unlist(list(getBM(attributes='hgnc_symbol',filters = 'hgnc_symbol', values=genes, mart = ensembl)$hgnc_symbol))
genes <- genes[order(as.character(genes), decreasing=T)]
######################FUNCTIONS####################
################################################
###function defs

#remove duplicate studies from a list of studies ( keeping those that are not published, i.e. provisional)
find_dups <- function(list_of_ids) { 
  index <- list()
  for (i in list_of_ids) { 
    p1 <- strsplit(i, '_')[[1]][1] 
    if (length(grep(p1, list_of_ids)) >1 ){ 
      patt <- 'pub'
      id <- as.character(grep(patt,list_of_ids,value=TRUE, invert=TRUE))
      index <- c(index,id)
    }
    else {
      i <- as.character(i)
      index <- c(index,i)
    }
  }
  return(unique(unlist(index)))
}

# returns a data frame of mutations and copy numbers for each tumor for a given gene
getInfo <- function(study_id, genes) { 
  study_id <- study_id
  #caseList = paste(study_id, '_cnaseq', sep='') #cases with gistic and mutation data
  study_info <- getCaseLists(mycgds, study_id)[,1]
  caseList <- grep('cna', study_info, value=TRUE, fixed=TRUE)
  if (length(caseList ==1)){
    if (length(grep('cnaseq', caseList) !=0)) {c
      gistic <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_gistic', sep=''), caseList=caseList)
      mutation <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_mutations', sep=''), caseList=caseList)
      m <- merge(gistic, mutation, by = 'row.names')
      return(study_id=m)
    }
  }
  if (length(grep('cna_seq', caseList)!= 0)) {
    gistic <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_cna', sep=''), caseList=caseList)
    mutation <- getProfileData(mycgds, genes=genes, geneticProfiles = paste(study_id, '_mutations', sep=''), caseList=caseList)
    m <- merge(gistic, mutation, by = 'row.names')
    return(study_id=m)
  }
  else{ 
    string <- paste('Failed-NoValid _cnaseq set', study_id, sep=' ' )
    print(string)
    return(NA)
  }
}
makeMutFrame <- function(genes, mp, gp, caseList) {
  #genes is a list of genes
  mutations <- getProfileData(mycgds, caseList=caseList, geneticProfile=mp, genes=genes )
  gistic <- getProfileData(mycgds, caseList = caseList, geneticProfile=gp, genes=genes)
  mutations$tumor <- rownames(mutations)
  gistic$tumor <- rownames(gistic)
  m_gistic <- melt(gistic, id.vars='tumor')
  m_mutations <- melt(mutations, id.vars='tumor')
  colnames(m_gistic) <- c('tumor', 'gene', 'gistic')
  colnames(m_mutations) <- c('tumor','gene','mutation')
  m <- merge(m_gistic, m_mutations, by=c('tumor', 'gene'))
  m <-as.data.frame(apply(m, c(1,2), replaceNaN))
  m$gistic <- as.numeric(as.character(m$gistic))
  print(head(m))
  tumor<- ddply(m, .(tumor, gene), summarize, mutated=isMut(c(mutation,gistic)))
  return(tumor)
  #return is: 
  #tumor #gene #mutated              
}

#for a given tumor classifies it as mutant or not
#row[2] is the gistic data ( mut is homozygous loss or high amplification)
#row[3] is the mutation data (no mutation is labelled as NA)
isMut <- function(row) {
  mut <- ifelse(as.numeric(row[2]) < -1 | as.numeric(row[2]) > 1 | !is.na(as.character(row[1])), 1, 0 )
  return(mut)
}

#classify if the mutation overlaps or not
mutType <- function(x, genes) {
  if (x[1] == 1 & x[2] ==0){ a <- eval(genes[1])}
  if (x[1] == 0 & x[2] ==1) {a <- eval(genes[2])}
  if (x[1] ==1 & x[2] ==1) { a <- 'BOTH'}
  if (x[1] ==0 & x[2] ==0 ) {a <- 'NONE'}
  return(a)
}

categorize = function(row){
  if (row[,1]== 1 & row[,2]==1) { 
    return('BOTH')
  }
  if (row[,1] == 1 & row[,2]==0) { 
    return('GENE_A')
  }
  if (row[,1] == 0 & row[,2] == 1) { 
    return('GENE_B')
  } 
  if (row[,1]== 0 &row[,2] == 0) {
    return('NEITHER')
  }
}

#converts the study_ids to the actual english names - prettier looking
mapNames <- function(study_id) { 
  #maps study_id to printable names
  require(cgdsr)
  mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
  studyMapping <- getCancerStudies(mycgds)[,c(1:2)]
  name <- studyMapping[studyMapping$cancer_study_id %in% study_id, ]$name
  #name <- unlist(strsplit(name, "\\(", ))[[1]]
  return(name)
}

replaceNaN <- function(cell){ 
  if (is.na(cell) ) {return(NA)}
  if (cell == 'NaN' ) {cell=NA}
  return(cell)
}
#####################WORKworkWORK####################################
#do the work 
#for (cancer_type in cancer_study_list){ 


#get tumor information for each gene on my list int he cancer I care about. 
study_info <- getCaseLists(mycgds, cancer_type)[,1]
caseList <- grep('cnaseq', study_info, value=TRUE, fixed=TRUE)
mp <- getGeneticProfiles(x=mycgds, cancerStudy=cancer_type )[,1]
mp <- grep('mutations', mp, value=TRUE)
gp <- getGeneticProfiles(x=mycgds, cancerStudy=cancer_type)[1,1]

df <- makeMutFrame(genes, mp, gp, caseList)
df <- dcast(df, tumor~gene, value.var='mutated')
df <- df[, !apply(is.na(df), 2, all)]
overlap_mat <- data.frame(matrix(NA, nrow=length(names(df)), ncol=length(names(df))))
colnames(overlap_mat) <- names(df)
rownames(overlap_mat) <- names(df)

#######################################################
#compare the genes pairwise to get overlap percent( that is all I care about in this case)
genes1 <- names(df)
genes2 <- names(df)


num_df <- data.frame(matrix(NA, ncol=4))
colnames(num_df) <- c('i', 'j', 'Var1', 'Freq')
numtests = 0
for (i in genes1[2:length(genes1)]){
  for (j in genes2[2:length(genes1)]) { 
    profiles <-data.frame(cbind(df[,i], df[,j]) )
    profiles$tumor = factor(df$tumor)
    colnames(profiles) = c(i, j, 'tumor')
    p = cbind(profiles[,1], profiles[,2])
    calcs <- as.data.frame(table(factor(rowSums(p), levels=0:2)))
    num <- ddply(profiles, .(tumor),  categorize)
    colnames(num)[2] <- 'class'
    overlaps <- calcs[calcs$Var1==2, 'Freq']
    percent_overlap <-overlaps/nrow(profiles)
    overlap_mat[i,j] <- percent_overlap
    num$gene_a = i
    num$gene_b = j 
    numbers = as.data.frame(table(factor(num$class, levels=c('BOTH', 'GENE_A', 'GENE_B', 'NEITHER'))))
    d = cbind(i,j,numbers)
    m <- matrix(c(numbers$Freq[4], numbers$Freq[3], numbers$Freq[2], numbers$Freq[1]), ncol=2)
    pv = fisher.test(m)$p.value
    numtests =+ 1
    if (pv < 0.05) { 
      print(paste(i,':', j, sep=''))
      print(fisher.test(m))
    }
    num_df = rbind(d, num_df)  
    #add odds ration and pvalue to print 
  }
}
#correct pvals
num_df = num_df[complete.cases(num_df), ]  
fname = c('~/Dropbox/', cancer_type, '_numbers.tsv', sep='\t')
pvalname = c('~/Dropbox/', cancer_type, '_pval.tsv', sep='\t')
#write.table(num_df, fname,row.names=F,  col.names=T, sep='\t')  
#write.table(pvals, pvalname , row.names=F, col.names=T, sep='\t')
######################################
#do some plotting
tab <- overlap_mat
tab$gene <-factor(row.names(tab))
tm <- melt(tab)
#a little magic to remove one triangle 
tm.lower <- tm[lower.tri(overlap_mat, diag=TRUE),]
#tm.lower <- subset(tm[lower.tri(overlap_mat),], gene !=variable) #this version makes the diagonal go away
#tm.upper <- subset(tm[upper.tri(overlap_mat),], gene != variable) 


g <-ggplot(tm.lower, aes(x=gene, y=variable, fill=value))+ geom_tile() 
g <- g + geom_text(label=ifelse(tm.lower$value > 0.05, round(tm.lower$value, 2), ''), sieze=2)
g <- g + scale_fill_continuous(low='white', high='#3399FF', limits=c(0,1)) + xlab('') + ylab('')
g <- g + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=8))
g<- g +  coord_fixed() + labs(fill='') + theme(axis.text.x=element_text(angle=90))
g
  