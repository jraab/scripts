#calculate the numbers needed to make a 3way venn diagram. 
#each function returns a venn diagram object that can be plotted with
#grid.newpage()
#grid.draw(venndiagram)

calc2way <- function(list1, list2, ...) { 
  require(VennDiagram)
  require(RColorBrewer)
  area1 <- length(list1)
  area2 <- length(list2)
  cross.area <- table(factor(list1 %in% list2, levels=c('TRUE', 'FALSE')))[1]
  vd <- draw.pairwise.venn(area1, area2, cross.area, euler.d = TRUE,
                     scaled = TRUE, inverted = FALSE, ext.text = TRUE, ext.percent = rep(0.05, 3),
                     lwd = rep(2, 2), lty = rep("solid", 2), col = rep("black", 2), fill =c("#E69F00", "#56B4E9"),
                     alpha = rep(0.5, 2), label.col = rep("black", 3), cex = rep(1, 3),
                     fontface = rep("plain", 3), fontfamily = rep("serif", 3), cat.pos = c(-50, 50),
                     cat.dist = rep(0.025, 2), cat.cex = rep(1, 2), cat.col = rep("black", 2),
                     cat.fontface = rep("plain", 2), cat.fontfamily = rep("serif", 2),
                     cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos = "outer",
                     cat.prompts = FALSE,
                     ext.pos = rep(0, 2), ext.dist = rep(0, 2), ext.line.lty = "solid",
                     ext.length = rep(0.95, 2), ext.line.lwd = 1, rotation.degree = 0,
                     rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0, ...)
  return(vd)
}

calc3way <- function(list1, list2, list3, ...) { 
  require(VennDiagram)                                        
  #each list should contain the ids of genes that change
  n12 <- table(factor(list1 %in% list2, levels=c('TRUE', 'FALSE')))[1] 
  n13 <- table(factor(list1 %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n23 <- table(factor(list2 %in% list3,levels=c('TRUE', 'FALSE')))[1]
  l1_l2_names <- list1[list1%in% list2] 
  n123 <- table(factor(l1_l2_names %in% list3, levels=c('TRUE', 'FALSE')))[1]
  area1 <- length(list1) 
  area2 <- length(list2) 
  area3 <- length(list3) 
  colors <- brewer.pal(3, 'Dark2')
  vd <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                   rotation = 1, reverse = FALSE, euler.d = TRUE,
                   scaled = TRUE, lwd = rep(2, 3), lty = rep("solid", 3),
                   col = rep("black", 3), fill = colors, alpha = rep(0.5, 3),
                   label.col = rep("black", 7), cex = rep(1, 7), fontface = rep("plain", 7),
                   fontfamily = rep("serif", 7), cat.pos = c(-40, 40, 180),
                   cat.dist = c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
                   cat.cex = rep(1, 3), cat.fontface = rep("plain", 3),
                   cat.fontfamily = rep("serif", 3),
                   cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos = "outer",
                   cat.prompts = FALSE, rotation.degree = 0, rotation.centre = c(0.5, 0.5),
                   ind = TRUE, sep.dist = 0.05, offset = 0, ...)
  return(vd)
}

calc4way <- function(list1, list2, list3, list4, ...) {
  require(VennDiagram)                                        
  #each list should contain the ids of genes that change
  n12 <- table(factor(list1 %in% list2, levels=c('TRUE', 'FALSE')))[1]
  n13 <- table(factor(list1 %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n14 <- table(factor(list1 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n23 <- table(factor(list2 %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n24 <- table(factor(list2 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n34 <- table(factor(list3 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  l1_l2_names <- list1[list1%in% list2] 
  l1_l3_names <- list1[list1%in% list3] 
  l1_l4_names <- list1[list1%in% list4]
  l2_l3_names <- list2[list2 %in% list3] 
  l2_l4_names <- list2[list2 %in% list4 ] 
  n123 <- table(factor(l1_l2_names %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n124 <- table(factor(l1_l2_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n134 <- table(factor(l1_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n234 <- table(factor(l2_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  l1_l2_l3_names <- list1[l1_l2_names %in% list3]
  n1234 <- table(factor(l1_l2_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  area1 <- length(list1) 
  area2 <- length(list2) 
  area3 <- length(list3) 
  area4 <- length(list4) 
  vd <- draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124,
                n134, n234, n1234, lwd = rep(2, 4), lty = rep("solid", 4),
                col = rep("black", 4), fill = brewer.pal(4, 'Dark2'), alpha = rep(0.5, 4),
                label.col = rep("black", 15), cex = rep(1, 15), fontface = rep("plain", 15),
                fontfamily = rep("serif", 15), cat.pos = c(-15, 15, 0, 0),
                cat.dist = c(0.22, 0.22, 0.11, 0.11), cat.col = rep("black", 4),
                cat.cex = rep(1, 4), cat.fontface = rep("plain", 4),
                cat.fontfamily = rep("serif", 4), cat.just = rep(list(c(0.5, 0.5)), 4),
                rotation.degree = 0, rotation.centre = c(0.5, 0.5), ind = TRUE, ...)
  return(vd)
}

calc5way <- function(list1, list2, list3, list4, list5, ...) {
  require(VennDiagram)                                        
  #each list should contain the ids of genes that change
  n12 <- table(factor(list1 %in% list2, levels=c('TRUE', 'FALSE')))[1]
  n13 <- table(factor(list1 %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n14 <- table(factor(list1 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n15 <- table(factor(list1 %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n23 <- table(factor(list2 %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n24 <- table(factor(list2 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n25 <- table(factor(list2 %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n34 <- table(factor(list3 %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n35 <- table(factor(list3 %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n45 <- table(factor(list4 %in% list5, levels=c('TRUE', 'FALSE')))[1]
  l1_l2_names <- list1[list1%in% list2] 
  l1_l3_names <- list1[list1%in% list3] 
  l1_l4_names <- list1[list1%in% list4]
  l2_l3_names <- list2[list2 %in% list3] 
  l2_l4_names <- list2[list2 %in% list4 ] 
  l3_l4_names <- list3[list3%in% list4] 
  n123 <- table(factor(l1_l2_names %in% list3, levels=c('TRUE', 'FALSE')))[1]
  n124 <- table(factor(l1_l2_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n125 <- table(factor(l1_l2_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n134 <- table(factor(l1_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n135 <- table(factor(l1_l3_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n145 <- table(factor(l1_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n234 <- table(factor(l2_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n235 <- table(factor(l2_l3_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n245 <- table(factor(l2_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n345 <- table(factor(l3_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  l1_l2_l3_names <- l1_l2_names[l1_l2_names %in% list3]
  l1_l2_l4_names <- l1_l2_names[l1_l2_names %in% list4]
  l1_l3_l4_names <- l1_l3_names[l1_l3_names %in% list4] 
  l2_l3_l4_names  <- l2_l3_names[l2_l3_names %in%list4]
  n1234 <- table(factor(l1_l2_l3_names %in% list4, levels=c('TRUE', 'FALSE')))[1]
  n1235 <- table(factor(l1_l2_l3_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n1245 <- table(factor(l1_l2_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n1345 <- table(factor(l1_l3_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  n2345 <- table(factor(l2_l3_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  l1_l2_l3_l4_names <- l1_l2_l3_names[l1_l2_l3_names %in% list4] 
  n12345 <- table(factor(l1_l2_l3_l4_names %in% list5, levels=c('TRUE', 'FALSE')))[1]
  area1 <- length(list1) 
  area2 <- length(list2) 
  area3 <- length(list3) 
  area4 <- length(list4) 
  area5 <- length(list5)
  colors = brewer.pal(5, 'Dark2')
  vd <- draw.quintuple.venn(area1, area2, area3, area4, area5, n12, n13, n14, n15, n23, n24, n25,
                      n34, n35, n45, n123, n124, n125, n134, n135, n145, n234, n235, n245, n345, n1234,
                      n1235, n1245, n1345, n2345, n12345,lwd = rep(2, 5),
                      lty = rep("solid", 5), col = rep("black", 5), fill = colors, alpha = rep(0.5, 5),
                      label.col = rep("black", 31), cex = rep(1, 31), fontface = rep("plain", 31),
                      fontfamily = rep("serif", 31), cat.pos = c(0, 287.5, 215, 145, 70),
                      cat.dist = rep(0.2, 5), cat.col = rep("black", 5), cat.cex = rep(1, 5),
                      cat.fontface = rep("plain", 5), cat.fontfamily = rep("serif", 5),
                      cat.just = rep(list(c(0.5, 0.5)), 5), rotation.degree = 0,
                      rotation.centre = c(0.5, 0.5), ind = TRUE, ...)
  return(vd)
}

