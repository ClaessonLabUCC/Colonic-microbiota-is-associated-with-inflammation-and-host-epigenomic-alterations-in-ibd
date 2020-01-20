heatplot_editd=function (dataset, dist_in, dend = c("both", "row", "column", "none"), 
                         cols.default = TRUE, lowcol = "green", highcol = "red", scale = "none", 
                         classvec = NULL, classvecCol = NULL, classvec2 = NULL, distfun = NULL, 
                         returnSampleTree = FALSE, method = "ave", dualScale = TRUE, 
                         zlim = c(-3, 3), scaleKey = TRUE, ...) 
{
  library(gplots)
  data <- array2ade4(dataset)
  data <- as.matrix(data)
  if (dualScale) {
    print(paste("Data (original) range: ", round(range(data), 
                                                 2)[1], round(range(data), 2)[2]), sep = "")
    data <- t(scale(t(data)))
    print(paste("Data (scale) range: ", round(range(data), 
                                              2)[1], round(range(data), 2)[2]), sep = "")
    data <- pmin(pmax(data, zlim[1]), zlim[2])
    print(paste("Data scaled to range: ", round(range(data), 
                                                2)[1], round(range(data), 2)[2]), sep = "")
  }
  distEisen <- function(x, use = "pairwise.complete.obs") {
    co.x <- cor(x, use = use,method= "spearman")
    dist.co.x <- 1 - co.x
    return(as.dist(dist.co.x))
  }
  cols <- function(low = lowcol, high = highcol, ncolors = 123) {
    low <- col2rgb(low)/255
    if (is.character(high)) 
      high <- col2rgb(high)/255
    col <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                        high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
    return(col)
  }
  cols.gentleman <- function() {
    library(RColorBrewer)
    hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    return(rev(hmcol))
  }
  if (cols.default) 
    plotcols = cols.gentleman()
  else plotcols = cols()
  if (is.null(distfun)) 
    distf = distEisen
  else distf = function(x) dist(t(x), method = distfun)
  Colv <- FALSE
  Rowv <- FALSE
  dend <- match.arg(dend)
  dend = tolower(dend)
  if (dend %in% c("row", "r")) {
    print(dend)
    Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
  }
  if (dend %in% c("column", "col", "c")) {
    print(dend)
    Colv = as.dendrogram(hclust(dist_in, method = method))
    dend = "column"
  }
  if (dend %in% c("both", "TRUE")) {
    Colv = as.dendrogram(hclust(dist_in, method = method))
    Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
    dend = "both"
  }
  RSideColors = CSideColors = NULL
  if (any(!is.null(classvec), !is.null(classvec2))) {
    proc.classvec <- function(classvec) {
      classvec = as.factor(classvec)
      if (is.null(classvecCol)) 
        classvecCol = getcol(length(levels(classvec)))
      SideCols = factor(classvec, labels = classvecCol)
      print(cbind(Class = levels(classvec), Color = levels(SideCols)))
      SideCols = as.character(SideCols)
      nSC = length(SideCols)
      return(list(nSC, SideCols))
    }
    if (!is.null(classvec)) {
      out = proc.classvec(classvec)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data)) 
        print("Error: classvec length not equal to nrow or ncol in data")
      if (nSC == nrow(data)) 
        RSideColors = SideCols
      if (nSC == ncol(data)) 
        CSideColors = SideCols
    }
    if (!is.null(classvec2)) {
      out = proc.classvec(classvec2)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data)) 
        print("Error: classvec2 length not equal to nrow or ncol in data")
      if (nSC == nrow(data)) 
        RSideColors = SideCols
      if (nSC == ncol(data)) 
        CSideColors = SideCols
    }
  }
  if (all(is.null(RSideColors), is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              zlim = zlim, dendrogram = dend, ...)
  if (all(!is.null(RSideColors), is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              RowSideColors = RSideColors, dendrogram = dend, zlim = zlim, 
              ...)
  if (all(is.null(RSideColors), !is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              ColSideColors = CSideColors, dendrogram = dend, zlim = zlim, 
              ...)
  if (all(!is.null(RSideColors), !is.null(CSideColors))) 
    heatmap.2(data, Colv = Colv, Rowv = Rowv, col = plotcols, 
              scale = scale, trace = "none", density.info = "none", 
              RowSideColors = RSideColors, ColSideColors = CSideColors, 
              zlim = zlim, dendrogram = dend, ...)
  if (returnSampleTree) 
    return(Colv)
}