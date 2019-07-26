library(pheatmap)

library(RColorBrewer)
myheatmap <- function(x,ncol=49,...) {
    mx <- max(abs(x))
    cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(ncol)
    cols[(ncol+1)/2] <- "#cccccc"
    pheatmap(x,
             color = cols,
             breaks=seq(-mx,mx,length.out=ncol+1),
             ...)
}

