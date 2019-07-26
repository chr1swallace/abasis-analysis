##' MAKE FIGURE 2
##' clusters of significant projected traits

## basis p values
## install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(ggplot2)
library(cupcake)
library(parallel)
library(pheatmap)
library(cluster)
library(dendextend)
library(factoextra)
library(cowplot)
library(gridGraphics) # recordPlot

renamer <- function(x) {
    x <- sub("_19","",x)
    x[x=="jia_case"] <- "jia_combined"
    x[x=="mpo"] <- "vasculitis"
    w <- grep("myogen",x)
    x[w]  %<>%  sub("_myogen","",.)  %>% paste0("myositis_",.)
    x[x=="myositis_myositis"] <- "myositis_combined"
    x[x=="egpa"] <- "EGPA_combined"
    x[x %in% c("anca_Neg","mpo_Pos")]  %<>%  paste0("EGPA_",.)
    x %<>% sub("_Pos|Pos|pos$","+",.)
    x %<>% sub("_Neg|Neg|neg$","-",.)
    x[grepl("copd",x)] <- "COPD"
    x
}


patt=NULL      

dplotter <- function(b,z,p,fixed=NULL,patt=NULL,what=c("signfdr","b"),
                     dist.method="euclidean",hclust.method="ward.D2",k=4,
                     pal=redblue,show.legend=TRUE) {
    what <- match.arg(what)
    if(!is.null(patt))
        fixed <- c(fixed,grep(paste(patt,collapse="|"),rownames(b),value=TRUE))
    if(!is.null(fixed)) {
        fixed <- intersect(fixed,rownames(b))
        b <- b[fixed,]
        z <- z[fixed,]
        p <- p[fixed,]
    }
    db <- get_dist(b,method=dist.method)
    d <- hclust(db,method=hclust.method)
    if(what=="signfdr") {
        zp <- cut(sign(z) * (1-p),breaks=c(-2,-0.99,-0.95,-0.9,-0.8,0,0.8,0.9,0.95,0.99,2),include.lowest=TRUE)  %>%
          as.numeric()  
        cols <- matrix(pal(10)[zp], 
                       nrow(z),ncol(z),dimnames=dimnames(z))
    } else {
        tb <- t(b)
        tb <- (tb/apply(abs(tb),1,max))  %>% t()
        zp <- cut(tb ,breaks=seq(-1,1,by=0.1),include.lowest=TRUE)  %>%  t()  %>% as.numeric()  
        cols <- matrix(pal(21)[zp], nrow(b),ncol(b),dimnames=dimnames(b))
    }
    
    par(mfrow=c(1,1))
    
    dd <- as.dendrogram(d)
    cuts <- cutree(dd,k=k)
    leafcols <- tol7qualitative[cuts[labels(dd)]]
    labels_colors(dd) <-  trait2col[ labels(dd) ]
    labels(dd)  %<>% renamer()
    rownames(cols)  %<>% renamer()
    dd %<>%
      set("leaves_pch", 19)  %>%
      ## set("leaves_col", trait2col[ labels(dd) ])  %>%
      set("leaves_col", leafcols) %>%
      set("branches_k_color", unique(leafcols), k = k)
    plot(dd,axes=FALSE)

    ## add coloured bars
    M <- ncol(b)
    ord <- structure(d$order,names=d$order.lab)
    colored_boxes(colors = cols[labels(dd),M:1],dend=dd,sort_by_labels_order=FALSE)
    
    ## tapply(p[fixed,],zp,summary)
    if(show.legend) {
        if(what=="signfdr") {
        legend("topleft",fill=pal(10)[c(1:4,10:7)],
               ncol=2,bty="n",
               text.width=0,
               legend=c(rep("",4),c("<0.01","<0.05","<0.1","<0.2")),
               title="FDR component association")
    } else {
        legend("topright",fill=pal(21)[c(1,11,21)],
               bty="n",
               ## legend=c("1","0.5","0",""),
               legend=c("-","0","+"),
               title="component value")
    }
        }
    invisible(d)
} 
## dplotter(b,z,p,fixed=tr)

################################################################################

J <- 1
source("~/Projects/auto-basis/scripts/read-results.R")
W <- reader("weight")
NW <- reader("noweight")

data <- W # so colours can use this
source("~/Projects/auto-basis/scripts/colours.R")
## BB vs basis on weight vs noweight
cats.keep <- c("bb_disease","bb_medications")
cats.drop <- setdiff(W$proj$category,c("bb_disease","bb_medications"))
THR <- 1 
L <- list(weight=W,noweight=NW)
for(nm in names(L)) {
    data <- L[[nm]]
    val <- "delta"
    b=getter(data,val,pthr=THR,
         cats.keep=cats.keep,
         cats.drop=cats.drop
         )
    rows.keep <- c(rbind(BB_LU,names(BB_LU)))  %>% unlist()  %>% intersect(.,rownames(b))
    ## dplotter(b[rows.keep,],z[rows.keep,],p[rows.keep,],what="b",hclust.method="ward.D2",
    m <- match(names(bb.renames),rows.keep)
    rows.keep[m] <- bb.renames
    m <- match(names(bb.renames),rownames(b))
    rownames(b)[m] <- bb.renames
    pdf(paste0("~/basis-hclust-",nm,".pdf"), height=6,width=6)
    par(mar = c(15,2,1,1))
    dplotter(b[rows.keep,],NULL,NULL,what="b",hclust.method="ward.D2",
             pal=grnvi,
             k=1,show.legend=nm=="weight")
    dev.off()
}
