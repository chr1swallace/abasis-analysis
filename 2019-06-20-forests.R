library(data.table)
library(magrittr)
library(ggplot2)
source("~/Projects/auto-basis/scripts/read-results.R")
res <- reader(cats.keep=c(
"ferreira_asthma"  ,       "tian_infectious_disease",
"bowes_jia_2019"    ,      "bb_cancer"              ,
"bb_disease"         ,     #"bb_medications"         ,
#"astle_blood"         ,    "Ig.titre"               ,
"lee_CD_prognosis"     ,   "psyc_consortium"        ,
"myogen"                ,  "estrada_NMO"            ,
"lyons_vasculitis","lyons_egpa"    ,         # "ad-pid"                 ,
"bowes_psa"      ,         #"rhodes_pah"             ,
"taylor_mtx"      ,        "kiryluk_iga_neph"       ,
"brown_as"         ,       "kuiper_bs"              ,
"cousminer_lada"     ,     "li_as"             
))
res.DT <- res$proj
## unique(res.DT$trait)


## add in basis traits for comparison
basis.DT <- data.table(trait=rownames(res$basis),res$basis) %>%
  melt(.,id.vars='trait',value.name="delta")
basis.DT[,category:="basis"]
## rename to match projected data
nm <- c(CD="crohns disease",
        UC="ulcerative colitis",
        PSC="primary sclerosing cholangitis",
        CEL="coeliac disease",
        MS="multiple sclerosis",
        RA="rheumatoid arthritis",
        PBC="primary billiary cholangitis",
        T1D="type 1 diabetes",
        VIT="vitiligo",
        SLE="systemic lupus erythematosis")
for(i in seq_along(nm))
    basis.DT[trait==names(nm)[i], trait:=nm[i]]
## res.DT <- rbind(res.DT,basis.DT,fill=TRUE)


## mangle trait names
tsubs <- c("jdm myogen"="Juvenile dermatomyositis",
           "pm myogen"="Polymyositis",
           "dm myogen"="Dermatomyositis",
           "myositis myogen"="Myositis",
           "eo"="Eosiniphil count",
           "eo baso sum"="Eosiniphil + basophil",
           "lymph"="Lymphocyte count",
           "wbc"="White blood cell count",
           mpo="Vasculitis: MPO+",
           "anca Neg"="EGPA: ANCA-",
           "egpa"="EGPA",
           "mpo Pos"="EGPA: MPO+")
for(i in seq_along(tsubs))
    res.DT[trait==names(tsubs)[i], trait:=tsubs[i]]
res.DT[,trait:=as.character(trait)]
res.DT[,trait:=sub("bb_","",trait)]
res.DT[,trait:=sub("SRD.|SRM.","",trait)]
res.DT[,trait:=gsub("."," ",trait,fixed=TRUE)]
res.DT[,trait:=gsub("_"," ",trait,fixed=TRUE)]
res.DT[,trait:=sub("^malabsorption ","",trait)]
res.DT[trait=="systemic lupus erythematosis sle",trait:="systemic lupus erythematosis"]

## group ank spond
res.DT[category=="brown_as",trait:="AS (Brown)"]
res.DT[category=="li_as",trait:="AS (Li)"]
res.DT[category %in% c("brown_as","li_as"),category:="ank.spond"]


unique(res.DT$category)
rares <- c( "lee_CD_prognosis" ,"bowes_jia_2019", "estrada_NMO", "lyons_egpa", "lyons_vasculitis","myogen", "kiryluk_iga_neph", "cousminer_lada" , "ank.spond"  )
## everything else will be merged
res.DT[ !(category %in% c(rares,"basis")), category:="UKBB/GWAS"]
sort(unique(res.DT$category))
csubs <- c(basis="Basis",
           "UKBB/GWAS",
           bb_medications="UK Biobank medications",
           bowes_psa="PSA",
           psyc_consortium="UKBB/GWAS",
           "bowes_jia_2019"="JIA",
           estrada_NMO="Neuromyelitis optica",
           myogen="MYOGEN",
           "tian_infectious_disease"="Infectious disease",
           lyons_vasculitis="Vasculitis",
           lyons_egpa="EGPA",
           astle_blood="Blood cells")
for(i in seq_along(csubs)) {
    res.DT[category==names(csubs)[i], category:=csubs[i]]
    basis.DT[category==names(csubs)[i], category:=csubs[i]]
}
## cols["Vasculitis"] <- cols["UK Biobank cancer"]
ctl.DT <- unique(res.DT[,.(variable,control.loading)])
all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## plot

## set one colour per category
library(RColorBrewer)
## these will have their own colour
cols <- c(brewer.pal(9,"Set1")[-6])
cols <- c(cols,"#333333","#ff0000") # put Basis=black,GWAS=grey
if(length(unique(res.DT$category)) > length(cols)) {
    print(sort(unique(res.DT$category)))
    print(length(cols))
    stop("colour mismatch")
}
names(cols) <- c(sort(unique(res.DT$category)),"Basis")

btraits <- unique(basis.DT$trait)
source("~/Projects/auto-basis/scripts/forests.R")
theme_set(theme_cowplot(16))
## bigr <- function(p,size.rel = 1) {
##   p +
##     theme(
##       plot.title    = element_text(size = rel(2.0 * size.rel), hjust = 0),
##       plot.subtitle = element_text(size = rel(1.5 * size.rel), hjust = 0),
##       legend.title  = element_text(size = rel(1.8 * size.rel)),
##       legend.text   = element_text(size = rel(1.3 * size.rel)),
##       axis.title    = element_text(size = rel(1.5 * size.rel)),
##       axis.text     = element_text(size = rel(1.2 * size.rel)),
##       strip.text.x  = element_text(size = rel(1.8 * size.rel)),
##       strip.text.y  = element_text(size = rel(1.8 * size.rel)),
##       legend.key.height = unit(1.3 * size.rel, "line"),
##       legend.key.width  = unit(1.3 * size.rel, "line")
##     )
## }

## all basic 
traits <- unique(res.DT[(p.adj<0.01 | category=="Basis" | category!="UKBB/GWAS")]$trait)# | category %in% rares,]$trait)


pdf("~/everything-forest.pdf",height=15,width=12)
for(i in 1:ncol(basis.DT)) {
    ipc=paste0("pc",i)
iPC=paste0("PC",i)
    p <- forest_everything(res.DT[trait %in% c(traits),],#,traits.i),],
                     basis.DT,pc=iPC,
                     focal=traits)
    print(p)
}
dev.off()
    

################################################################################

## individual customized
forest_plot(proj.dat=res.DT[category %in% c("UKBB/GWAS"),],
            basis.dat=basis.DT,
            pc="PC1",
            title="")  #%>% bigr()

## MYGEN, PC1
forest_plot(proj.dat=res.DT[category %in% c("UKBB/GWAS","MYOGEN")],
                  basis.DT,pc="PC1",focal=all.traits[c('MYOGEN')])
## AS, PC1
forest_plot(proj.dat=res.DT[category %in% c("UKBB/GWAS","ank.spond")],
                  basis.DT,pc="PC1",focal=all.traits[c('ank.spond')])
## JIA, PC3
forest_plot(proj.dat=res.DT[category %in% c("UKBB/GWAS","JIA")],
            basis.DT,pc="PC1",focal=all.traits[c('JIA')])

