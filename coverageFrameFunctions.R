#calculate coverage frames based on a .rda file
#this should live in my scripts repo and be sourced to a .Rmd file

calcCvgFrame <- function(frameOfRegions, coverageAll, revStrand, left, right){
  # this function expects frameOfRegions to be a proper 6 column bed
  # coverageAll is a rle vector 
  i<-1
  resultFrame <-matrix(data=rep(0),nrow=nrow(frameOfRegions),ncol=right-left + 1)
  #geneNames <- list()
  if (length(colnames(frameOfRegions)) != 6 ){
    stop('frameOfRegions is not a 6 colunn bed')
  }
  colnames(frameOfRegions) <- c('chr','start','end','name','score','strand')
  chromosomes <- unique(frameOfRegions$chr)
  for (chr in chromosomes){ 
    coverage <- coverageAll[names(coverageAll)==chr]
    coverage <- coverage[[1]]
    centersOnChr <- frameOfRegions[frameOfRegions$chr==chr,]
    winCenters <- centersOnChr$start
    for (j in seq(along=winCenters) ) {
      if (winCenters[j] + left <1 | winCenters[j] + right > length(coverage)) {next}
      df <- data.frame(chr=chr, start=winCenters[j] + left, end=winCenters[j]+ right)
      v <- as.vector( seqselect(coverage, winCenters[j]+left, winCenters[j] + right))
      if (revStrand[j] ) {
        v <- rev(v)
      }
      resultFrame[i,] <- v
      i <- i +1
    }
  }
  return(resultFrame)
}

rowBins<-function(row,binSize=200, log=TRUE){
  #this function normalizes a coverageframe by row
  if (log==TRUE){
    row <- log2(row+1)
  }
  numBins<-round(length(row)/binSize)
  outRow<-array(rep(0,numBins))
  for (i in 1:numBins){
    if (i==1)
    { window<-c(1,binSize)
    }
    else{
      window <- c((binSize*(i-1)),(binSize*i))
    }
    binSum <- sum(row[window[1]:window[2]])
    outRow[i]<-as.vector(binSum/binSize)
  }
  return(outRow)
}

