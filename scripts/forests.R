library(cowplot)
forest_plot <- function(proj.dat,basis.dat=basis.DT,pc,focal=NULL,title="",fdr_thresh=0.01){
    if(is.list(focal))
        focal  %<>%  unlist(.)
    dat <- proj.dat[variable==pc & ((p.adj<fdr_thresh ) |
                                    trait %in% unique(c(btraits,focal))),]
    ## dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
    ## for(i in seq_along(tsubs))
    ##     dat[trait==names(tsubs)[[i]],trait:=tsubs[i]]
    ## for(i in seq_along(csubs))
    ##     dat[category==names(csubs)[[i]],category:=csubs[i]]
    ocat <- c("Basis","UKBB/GWAS")
    ocat <- c(setdiff(unique(dat$category),ocat),ocat)
    dat[,category:=factor(category,levels=ocat)]
    dat[,newcat:=category]
    dat[category=="Basis",newcat:="UKBB/GWAS"]
  dat[,trait:=factor(trait,levels=unique(dat[order(newcat,delta,decreasing=TRUE),]$trait))]
    ## if(is.na(theme)){
  ##     theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),
  ##                    legend.position="none")
  ## }
        shapes <- ifelse(names(cols)=="Basis",17,19)
    names(shapes) <- names(cols)
ggplot(dat,aes(x=trait,y=delta,colour=category,pch=category,linetype=p.adj>fdr_thresh)) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=lower,ymax=upper)) +
    coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) +
      ggtitle(paste("Component",sub("PC","",pc))) +
      theme(plot.title=element_text(hjust = 0),
            legend.position=c(0,0),legend.justification=c(0,0),legend.box.background=element_rect(fill="lightyellow"),legend.title=element_blank()) +
      background_grid() +
      scale_colour_manual(values=cols) +
scale_linetype_discrete(guide="none") +
    xlab("Trait") + ylab("Basis score - control")
}

forest_everything <- function(proj.dat,basis.dat=basis.DT,pc,focal=NULL,title="",fdr_thresh=0.01){
    if(is.list(focal))
        focal  %<>%  unlist(.)
    dat <- proj.dat[variable==pc,]
                    ## & trait %in% unique(c(btraits,focal)),]
    ## dat[trait %in% focal,category:='aa_focal_diseases']
    dat[,ci:=1.96 * sqrt(variance)]
    dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
    dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
    ## for(i in seq_along(tsubs))
    ##     dat[trait==names(tsubs)[[i]],trait:=tsubs[i]]
    ## for(i in seq_along(csubs))
    ##     dat[category==names(csubs)[[i]],category:=csubs[i]]
    ## ocat <- c("Basis","UK Biobank disease","UK Biobank medications","UKBB/GWAS")
    ocat <- c("Basis","UKBB/GWAS")
    ocat <- c(setdiff(unique(dat$category),ocat),ocat)
    dat[,category:=factor(category,levels=ocat)]
    dat[,newcat:=category]
    dat[category=="Basis",newcat:="UKBB/GWAS"]
    dat[,trait:=factor(trait,levels=unique(dat[order(newcat,delta,decreasing=TRUE),]$trait))]
    ## if(is.na(theme)){
    ##     theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),
    ##                    legend.position="none")
    ## }
    shapes <- ifelse(names(cols)=="Basis",17,19)
    names(shapes) <- names(cols)
    ggplot(dat,aes(x=trait,y=delta,colour=category,pch=category,linetype=p.adj>fdr_thresh)) +
      geom_point(size=2) +
      geom_errorbar(aes(ymin=lower,ymax=upper)) +
      coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) +
      ggtitle(paste("Component",sub("PC","",pc))) +
      theme(plot.title=element_text(hjust = 0),
            legend.position=c(0,0),legend.justification=c(0,0),legend.box.background=element_rect(fill="lightyellow"),legend.title=element_blank()) +
      background_grid() +
      scale_shape_manual(values=shapes) +
      scale_colour_manual(values=cols) +
      scale_linetype_discrete(guide="none") +
      xlab("Trait") + ylab("Basis score - control")
}

