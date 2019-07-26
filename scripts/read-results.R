##' read results of basis construction and projection
##' 
##' loads res.DT, basis.DT and getter() to extract significant /
##' interesting subset of traits & make projection matrix

library(data.table)
library(magrittr)
library(cupcake)
library(parallel)
################################################################################

## RESULTS.FILE <- '~/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
## name of each element below
ti <- c("+vit (latest)","ten","+psa","vit+t2d","ten noshrink")
RESULTS.FILE <- c("~/share/as_basis/GWAS/RESULTS/17_07_19_0619_summary_results.RDS",
                  "~/share/as_basis/GWAS/RESULTS/03_07_19_0619_summary_results.RDS",
                  "~/share/as_basis/GWAS/RESULTS/24_06_19_0619_summary_results.RDS",
                  "~/share/as_basis/GWAS/RESULTS/09_04_19_summary_results.RDS",
                  "~/share/as_basis/GWAS/RESULTS/08_04_19_PSA_summary_results.RDS",
                "~/share/as_basis/GWAS/RESULTS/18_06_19__vit_t2d_summary_results.RDS")
BASIS_FILE <- c('~/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS', # latest
                '~/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS', # latest
                '~/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
                '~/share/as_basis/GWAS/support/psa_ss_basis_gwas.RDS',
                '~/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS',
                '~/share/as_basis/GWAS/support/ss_basis_noshrink.RDS') # noshrink
if(!exists("J"))
    J <- 1# use vit, not psa or t2d
message("loading basis ",ti[J])


bb.renames <- c("SRD.asthma"="UKBB_asthma",
             "SRD.crohns.disease"="UKBB_CD",
             "SRD.malabsorption.coeliac.disease"="UKBB_CEL",
             "SRD.multiple.sclerosis"="UKBB_MS", 
             "SRD.rheumatoid.arthritis"="UKBB_RA",
             "SRD.systemic.lupus.erythematosis.sle"="UKBB_SLE", 
             "SRD.type.1.diabetes"="UKBB_T1D",
             "SRD.ulcerative.colitis"="UKBB_UC",
             "SRD.vitiligo"="UKBB_VIT")

## UNWEIGHTED BASIS
proj.noweight <- "~/share/as_basis/GWAS/bb_projections/noweight_0619/" # - bb projections w/o weights
basis.noweight <- "~/share/as_basis/GWAS/support/noweight_basis_gwas_0619.RDS" # - unweighted basis


## renaming function
bbnames <- function(x) {
    as.character(x)  %>% sub("bb_","",.)  %>% make.names()
}

dt2mat <- function(dt,...) {
    tmp <- dcast(dt,...)
    rn <- tmp[[1]]
    m <- as.matrix(tmp[,-1])
    rownames(m) <- rn
    m
}
BB_LU <- c(
  CD = 'SRD.crohns.disease',
  CEL = 'SRD.malabsorption.coeliac.disease',
  MS = 'SRD.multiple.sclerosis',
  RA = 'SRD.rheumatoid.arthritis',
  SLE = 'SRD.systemic.lupus.erythematosis.sle',
  T1D = 'SRD.type.1.diabetes',
  UC = 'SRD.ulcerative.colitis',
  PBC = 'PBC',#SRM.ursodeoxycholic.acid',
  PSC = 'PSC',
VIT = 'SRD.vitiligo',
  asthma = 'SRD.asthma'
)
BB_LU <- BB_LU[order(names(BB_LU))]

centre.basis <- function(f) {
        pc.emp <- readRDS(f)
        t(pc.emp$x)[,"control"]
}
read.basis <- function(f) { # calculates delta
        pc.emp <- readRDS(f)
        t( t(pc.emp$x) - t(pc.emp$x)[,"control"] )[,1:(ncol(pc.emp$x)-1)]
}
mangle.proj <- function(tmp,cats.keep) {
        tmp[,trait:=bbnames(trait)]
        ## remove two duplicates
        tmp[,.dup:=.N,by=c("category","trait","variable")]
        tmp <- tmp[.dup==1,][,.dup:=NULL]
        ## deal with two same named traits, different cats
        tmp[trait %in% c("shingles","pneumonia"), trait:=paste(sub("_.*","",category),trait)]
        ## tmp[,star:=ifelse(p.adj<0.01, "*","")]
        copy(tmp)
}

reader <- function(what=c("weight","noweight"),
                   J=1,
                   cats.keep=c("ferreira_asthma",
                               "tian_infectious_disease",
                               "bowes_jia_2019", "bb_cancer",
                               "bb_disease", "bb_medications",
                               "brown_as",
                               "astle_blood",
                               "rhodes_pah",
                               "taylor_mtx",
                               "kiryluk_iga_neph",
                               "psyc_consortium",
                               "kuiper_bs",
                               "li_as",
                               "cousminer_lada",
                               "geneatlas_srd",
                               "myogen",
                               "estrada_NMO",
                               "lyons_egpa",
                               "bowes_psa")) {
    what <- match.arg(what)
    if(what=="weight") {
        tmp <- readRDS(RESULTS.FILE[[J]])
        basis <-  read.basis(BASIS_FILE[[J]])
        message("categories found:")
        print(unique(tmp$category))
        message("categories kept:")
        print(intersect(cats.keep,unique(tmp$category)))
        message("use arg cats.keep to change")
        return(list(basis=basis,proj=mangle.proj(tmp)[category %in% cats.keep,]))
    } else {
        nobasis <- read.basis(basis.noweight)
        ctl <- centre.basis(basis.noweight)
        tmp <- list.files(proj.noweight,full=TRUE)  %>%
          lapply(.,readRDS)  %>%
          rbindlist()  %>%
          melt(.,c("trait"))
        tmp[,trait:=make.names(trait)]
        tmp[,category:="bb_disease"]
        tmp[,control.loading:=ctl[variable]]
        tmp[,delta:=value - control.loading]
        return(list(basis=nobasis,proj=mangle.proj(tmp)))
    }
}
    

getter <- function(data,vv="delta",limit=TRUE,
                   cats.keep=c("bowes_jia_2019", "bowes_psa", "estrada_NMO", "myogen", "lyons_egpa","taylor_mtx","kiryluk_iga_neph", "taylor_mtx"),
                   cats.drop=c("astle_blood",
                               "bb_medications","wong_aav","Ig.titre","lee_CD_prognosis"),
                   pthr=0.01) {
    if(length(intersect(cats.keep,cats.drop)))
        stop("cats must be in either keep or drop, not both: ",intersect(cats.keep,cats.drop))
    if("p.adj" %in% names(data$proj))
        sig <-  data$proj[p.adj<pthr,]$trait  %>% unique()
    else
        sig <- ""
    tmp <- copy(data$proj)
    if(limit)
        tmp <- tmp[(trait %in% sig | category %in% cats.keep) &
                   !(category %in% cats.drop) ,]
    proj <- dt2mat(tmp,
                   trait ~ variable,value.var=vv)[,colnames(data$basis)]
    if(vv %in% c("z","p.adj")) {
        nullv <- if(vv=="z") { sign(basis)*0.1 } else { 0.5 }
        basis <- matrix(nullv,nrow(data$basis),ncol(data$basis),
                        dimnames=dimnames(data$basis))
    } else {
        basis <- data$basis
    }
    res <-  rbind(proj,basis)
    if("any_asthma_eczema_rhinitis" %in% data$proj$trait && vv %in% c("Z","delta")) {
        sw <- dt2mat(data$proj[trait=="any_asthma_eczema_rhinitis",],
                     trait ~ variable,value.var=vv)[,colnames(basis)]  %>% c()
        sw[4] <- -sw[4]
        res <- res * matrix(sign(sw),nrow(proj) + nrow(basis),
                            ncol(proj),byrow=TRUE)
    }
    res
}


hm <- function(data,i) {
    use <- with(data$proj[[i]],
                trait %in% sig[[i]] &
                !(category %in% c("astle_blood","bb_medications")))
    proj <- dt2mat(data$proj[[i]][use,],
                   trait ~ variable,value.var="delta")[,colnames(basis[[i]])]
    P <-  dt2mat(data$proj[[i]][use,],
                 trait ~ variable,value.var="star")[,colnames(basis[[i]])]
    proj <- rbind(proj,basis[[i]])
    P <- rbind(P,matrix("",nrow(basis[[i]]),ncol(basis[[i]]),dimnames=dimnames(basis[[i]])))
    ann <- as.data.frame(unique(data$proj[[i]][use,.(trait,category)]))
    rownames(ann) <- ann$trait
    ann <- ann[rownames(P),"category",drop=FALSE]
    pheatmap(proj,cluster_cols=FALSE,display_numbers=P,annotation_row=ann,method="ward")
}
