#!/usr/bin/env Rscript
#
# subtract_featureCount.R -a process/featureCounts/exon_genecounts.txt -b process/featureCounts/all_genecounts.txt
library(optparse)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

# debug start.
#fa <- fread("/scratch-cbe/users/jaritz/aid/2023Mar29/featureCounts/exon_genecounts.txt",sep="\t")
#fas <- fread("/scratch-cbe/users/jaritz/aid/2023Mar29/featureCounts/exon_genecounts.txt.summary",sep="\t")
#fb <- fread("/scratch-cbe/users/jaritz/aid/2023Mar29/featureCounts/all_genecounts.txt",sep="\t")
#fbs <- fread("/scratch-cbe/users/jaritz/aid/2023Mar29/featureCounts/all_genecounts.txt.summary",sep="\t")
#fc <- fread("/groups/busslinger/Kimon/genome_mm10/refGene.mm10.2018_1206.withIgtcr.final.annotation.IMP_IGTCR_only.with_introns.gtf")
#outname <- "test.txt"
#setwd("/scratch-cbe/users/jaritz/aid/2023Mar29/tmp")
# debug end.

option_list = list (
  make_option(
    c("-a","--fname_a"),
    type="character",
    default=NULL,
    help="first FeatureCount result table",
    metavar="character"
  ),
  make_option(
    c("-b","--fname_b"),
    type="character",
    default=NULL,
    help="second FeatureCount result table",
    metavar="character"
  ),
  make_option(
    c("-c","--fname_c"),
    type="character",
    default=NULL,
    help="GTF annotation file (exons)",
    metavar="character"
  ),
  make_option(
    c("-m","--metric"),
    type="character",
    default="subtract",
    help="output metric, subtract or concat",
    metavar="character"
  ),
  make_option(
    c("-o","--fname_o"),
    type="character",
    default=NULL,
    help="Output filename (without erxtension)",
    metavar="character"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$fname_a) ||
    is.null(opt$fname_b) ||
    is.null(opt$fname_c)) {
  print_help(opt_parser)
  stop("At least three files must be supplied ...\n",call.=FALSE)
}

#
# read the count files and associated summary files
#
# Geneid (equal) 
# Chr, Start, End, Strand, Length, sample names ...
fa <- fread(opt$fname_a,sep="\t")
fas <- fread(paste0(opt$fname_a,".summary"),sep="\t")
fb <- fread(opt$fname_b,sep="\t")
fbs <- fread(paste0(opt$fname_b,".summary"),sep="\t")
fc <- fread(opt$fname_c,sep="\t")

#
# Checks
#

# number of samples
if(dim(fa)[2]!=dim(fb)[2]) {
  stop("count files do not have equal number of samples! ...\n",call.=FALSE)
}

#
# we need common gene identifiers,
# so clean up diverging gene names (due to the gene merging step):
#

# clean up genes names unique in second (all gene featurecount file)
# 

# these are common, so no problem for merging both files ...
fb.common <- intersect(fb$Geneid,fa$Geneid)
print("common gene names:")
length(fb.common)

# these are unique in the second (genes), need to clean up to be able to merge
fb.unique <- setdiff(fb$Geneid,fa$Geneid)
print("unique gene names:")
length(fb.unique)

# handle _NM_, _NR_ : remove transcript id
print("handle _NM_, _NR_ in unique names ...")
nm.enames <- str_split(fb.unique[str_detect(fb.unique,"_NM_")],"_NM_",n=2,simplify = T)[,1]
nm.gnames <- fb.unique[str_detect(fb.unique,"_NM_")]
nr.enames <- str_split(fb.unique[str_detect(fb.unique,"_NR_")],"_NR_",n=2,simplify = T)[,1]
nr.gnames <- fb.unique[str_detect(fb.unique,"_NR_")]

# ig-like molecules: replace with names obtained from exons annotation
print("handle Ig-like (Gm|Tr|Ig|Tc|C92|B23) in unique names ..." )
xx <- fb.unique[str_detect(fb.unique,"Gm|Tr|Ig|Tc|C92|B23")]

# 1. original gtf: select annotation lines from exons
s <- (tibble(fc) %>% filter(fc$V3=="exon") %>% dplyr::select(V9))$V9
#length(s)

# 2. original gtf: extract genes and names
genes <- str_remove_all(str_remove_all(str_split(s,";",simplify = T)[,1],"gene_id "),"\"")
names <- str_remove_all(str_remove_all(str_split(s,"gene_name ",simplify = T)[,2],"\""),";")

# 3. original gtf: concat and make unique
genes.names <- unique(paste(genes,names))
#length(genes.names)
#head(genes.names)

# 4. original gtf: search
print("... searching for gene names in original annotation table... (~10 sec)")
m <- sapply(xx, grepl, genes.names, ignore.case=F)
#dim(m)
xx.genes <- str_split(genes.names[rowSums(m)>0]," ",n=3,simplify=T)[,2]
xx.names <- str_split(genes.names[rowSums(m)>0]," ",n=2,simplify=T)[,1]

# join all obtained lists
print("compiling and adding final gene names ...")
df <- unique(data.frame(gname=c(nm.gnames,nr.gnames,xx.genes,fb.common),
           ename=c(nm.enames,nr.enames,xx.names,fb.common)))
dim(df)

# add to gene featurecount table
fb.ext <- merge(fb,df,by.x="Geneid",by.y="gname",all.x = T,all.y=T)
print("exon dim (a):")
dim(fa)
print("genes dim(b):")
dim(fb)
print("extended genes dim: (incl. exon names):")
dim(fb.ext)

# these cases have not been handled (correctly?) in exon table,
# we set the ename again
print("exon transcripts mapped to multiple genes (replace by simple gene name:")
table(is.na(fb.ext$ename))  # none
#fb.ext[is.na(fb.ext$ename),c(1:6)]

# replace
#geneid.enames <- str_remove(fb.ext[is.na(fb.ext$ename),"Geneid"][[1]],"__.*")
#fb.ext[is.na(fb.ext$ename),"ename"] <- geneid.enames

# now should be 0
setdiff(fb.ext$ename,fa$Geneid)  # NA in fb.ext

# check (should be 0)
table(is.na(fa$Geneid))
setdiff(fa$Geneid,fb.ext$ename)  # Pira6 in fa, occurs on plus and minus on same strand

# unique gene ids (all 1)
tibble(fa) %>% dplyr::select(Geneid) %>% group_by(Geneid) %>% dplyr::count() %>% arrange(-n)

# unique gene ids (all 1)
tibble(fb) %>% dplyr::select(Geneid) %>% group_by(Geneid) %>% dplyr::count() %>% arrange(-n)

# some gene ids are now repeated up to three times,
# we need to sum up the counts via group by after merging ...
print("most frequent gene ids:")
geneids.mult <- tibble(fb.ext) %>% dplyr::select(Geneid) %>% group_by(Geneid) %>% dplyr::count() %>% filter(n>1) %>% arrange(-n)
geneids.mult

print("examples:")
fb.ext %>% filter(Geneid=="Ighe") %>% dplyr::select(c(1:6,dim(fb.ext)[2]))

# some exon ids are up to 24 times,
# so single exons are now spread out to independent genes,
# the exon counts will be spread out 
print("most frequent exon ids:")
ge <- tibble(fb.ext) %>% dplyr::select(ename) %>% group_by(ename) %>% dplyr::count()  %>% filter(n>1) %>% arrange(-n)
ge
fb.ext[ename %in% ge$ename,c("Geneid","ename")]

#
# Finally merge the two tables
# - sum up counts if multiple exons are combined to one gene
# - take the exon counta as is, if one exon os spread to multiple genes
colnames(fa) <- paste0(colnames(fa),".exon")
colnames(fb.ext) <- paste0(colnames(fb.ext),".all")

##
## concat
##
print("merging exon and extended genes table ...")
fd <- merge(fb.ext,fa,by.x="ename.all",by.y="Geneid.exon",all = T) %>%
  relocate(Length.exon,.after=Length.all)  %>%
  relocate(Strand.exon,.after=Length.all)  %>%
  relocate(End.exon,.after=Length.all)  %>%
  relocate(Start.exon,.after=Length.all)  %>%
  relocate(Chr.exon,.after=Length.all)
print("final table:")

dim(fd)
# annotation:
fd[1:3,1:12]
fd[1:3,13:dim(fd)[2]]

# handle duplicated gene names: sum up and make uniq (remove exon info)
fd[Geneid.all %in% geneids.mult$Geneid,c(1,2,13)]
fd.sum <- fd %>% group_by(Geneid.all) %>% mutate(across(13:dim(fd)[2]-1,sum))
fd.sum
fd.sum.newgeneids <- fd.sum %>% filter(Geneid.all %in% geneids.mult$Geneid) %>% mutate(Geneid.all=paste0(Geneid.all,"__",ename.all)) %>%
  select(Geneid.all)

# replace
fd[Geneid.all %in% geneids.mult$Geneid,"Geneid.all"] <- fd.sum.newgeneids$Geneid.all

# check
fd$Geneid.all[duplicated(fd$Geneid.all)]
dim(fd)
colnames(fd)

if(opt$metric=="concat") {
  #fwrite(fd,file=paste0("concat.txt"),quote=F,sep="\t")
  con <- file(opt$fname_o,open="wt")
  writeLines(paste0("# Program:subtract_featueCount.R; Command:\"subtract_featureCount.R -a ",opt$fname_a,
    " -b ",opt$fname_b,
    " -c ",opt$fname_c,
    " -m ",opt$metric,
    " -o ",opt$fname_o),
    con=con)
  write.table(fd,file=con,append=T,quote=F,sep="\t",col.names=T,row.names=F)
  close(con)
}

##
## subtract
##
print("subtracting exons from genes ...")

# preparing exon count table
t <- dim(fb)[2]-6                    # exon table dimensions, without annot
m.exon <- fd[,(12+t+1):dim(fd)[2]]    # second half of matrix
dim(m.exon)
colnames(m.exon)                    

# prepare genes table
t <- dim(fb)[2]-6                    # all table dimensions
m.all <- fd[,(12+1):(t+12)]
dim(m.all)
colnames(m.all)

# subtract
m.rest <- m.all - m.exon

# replace negative values that occure due to 
# cases of (exon) genes in (large) introns, for example
print("replacing negative values that occure due to exons of gene A  in introns of gene B:")
m <- as.matrix(m.rest)
print("% of effected values:")
100* (sum(as.numeric(m<0)) / ( sum(as.numeric(m>=0)) + sum(as.numeric(m<0))))

tibble(x=m[m<0]) %>% ggplot(aes(log10(abs(x)))) + geom_histogram(bins=50) + ggtitle("negative counts")

x <- tibble(x=fd$Geneid.all[apply(m,1,min)<=-500]) %>% 
  filter(!grepl("Rik",x)) %>%
  filter(!grepl("^Gm",x)) %>%
  filter(!grepl("^Mir",x)) %>%
  filter(!grepl("^Igh",x))  %>%
  filter(!grepl("_NR_",x))
x

m[m<0] <- 0

# clean colnames
colnames(m) <- str_replace_all(colnames(m),".all",".rest")
colnames(m)

# add annotation (featrueCount format)
fd[1:3,c(2:7)]
m.df <- cbind(fd[,c(2:7)],m)
colnames(m.df) <- str_remove_all(colnames(m.df),".all")
colnames(m.df) <- str_remove_all(colnames(m.df),".rest")
#colnames(m.df)[1] <- "Geneid"
colnames(m.df)
m.df[1:3,c(1:8)]

if(opt$metric=="subtract") {
  #fwrite(fd,file=paste0("subtract.txt"),quote=F,sep="\t")
  con <- file(opt$fname_o,open="wt")
  writeLines(paste0("# Program:subtract_featueCount.R; Command:\"subtract_featureCount.R -a ",opt$fname_a,
                " -b ",opt$fname_b,
                " -c ",opt$fname_c,
                " -m ",opt$metric,
                " -o ",opt$fname_o),
             con=con)
  write.table(m.df,file=con,append=T,quote=F,sep="\t",col.names=T,row.names=F)
  close(con)
}

# we assume that there is also a summary file
# fas
# fbs
if(opt$metric=="concat") {
  fas.df <- fas[,-1]
  colnames(fas.df) <- paste0(colnames(fas.df),".ex")
  fbs.df <- fbs[,-1]
  colnames(fbs.df) <- paste0(colnames(fbs.df),".all")
  df <- data.frame(fas[,1],fbs.df,fas.df)
  fwrite(df,file=paste0(opt$fname_o,".summary"),quote=F,sep="\t")
} else if (opt$metric=="subtract") {
  df <- data.frame(fas[,1],fbs[,-1] - fas[,-1])
  colnames(df) <- paste0(colnames(df),".rest")
  fwrite(df,file=paste0(opt$fname_o,".summary"),quote=F,sep="\t")
}
