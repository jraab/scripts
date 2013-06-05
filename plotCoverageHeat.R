#this file saves my favorite plotting methods. 

plotHeatmap <- function(mat, orderingMatrix,label, legend=TRUE, scaleMin,scaleMax){ 
  orderOfSums <- order(rowSums(orderingMatrix), decreasing=T)
  orderedMat <- as.data.frame(mat[orderOfSums,])
  orderedMat$order <- factor(orderOfSums, levels=unique(as.character(orderOfSums)))
  meltMat <- melt(orderedMat)
  colors <- brewer.pal(5, 'Blues')
  pal <- colorRampPalette(colors) 
  g <- ggplot(meltMat, aes(x=variable, y=order, fill=value)) + geom_tile() + ggtitle(label) 
  g  <- g + scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=pal(5))  
  g <- g + theme(axis.ticks=element_blank(), axis.text=element_blank())
  if (legend==FALSE) { 
    g <- g + theme(legend.position='none')
  } 
  return(g)            
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
