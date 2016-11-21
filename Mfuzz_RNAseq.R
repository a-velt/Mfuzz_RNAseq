#Author : Amandine Velt (amandine.velt@inra.fr)
#Date : 18/11/2016

#################################################################################################################################################################
# A the beginning, Mfuzz was developped to perform clustering of gene expression from microarray data.
# To adapt this method on RNA-seq data, the author suggest to do some additional preprocessing. For instance, starting from FPKMs (normalization by gene length)
# and exclude genes which do not show expression (i.e. with FPKM equals zero).
#
# This script takes as input a directory path containing all (and only) the RNA-seq raw count data tables (eg a directory with one htseq-count file per sample)
# and performs the DESeq normalization method (normalization by library size) and then divides these normalized counts by gene length (in kb) to obtain
# a RPKN (reads per kilobase number). After gene length normalization, this script performs the clustering of gene expression time-series RNA-seq data with Mfuzz.
#
# Usage in, e.g, Rstudio :
# ./Mfuzz_RNAseq -f count_files_folder -a annotation -b gene_name_attribute -t time -n nb_clusters -m membership_cutoff -o output_directory
#
# Arguments description :
# count_files_folder -> directory containing all the raw count data tables (one per sample)
# annotation -> a gtf or gff file allowing to calculate the genes length
# gene_name_attribute -> the name of the attribute in the gtf referring to the gene information
# time -> give the time value of each file by respecting the same order in the vector than the files in the folder.
#   if several files correspond to a same time (replicates), give the same time value and then the script performs the mean on the normalized counts of all the
#   samples of a same time
# nb_clusters -> number of clusters to generate with Mfuzz (empirical choice)
# membership_cutoff -> the membership cut-off to use with Mfuzz -> see the Mfuzz paper : http://www.bioinformation.net/002/000200022007.pdf
# output -> directory where store the results
#
# Examples of arguments to give :
# count_files_folder="/home/user/count_files_folder"
# annotation="HS.Genes.v2.gff"
# gene_name_attribute="gene"
# time="time1,time1,time1,time2,time2,time2,time3"
# nb_clusters = 4
# membership_cutoff = 0.7
# output = /home/user/cluster_output
#################################################################################################################################################################

# libraries dependencies
suppressMessages(library("optparse"))
suppressMessages(library("tools"))
suppressMessages(library("Mfuzz"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("DESeq"))
suppressMessages(library("edgeR"))

# options of the script
option_list = list(
  make_option(c("-f", "--folder"), type="character", default=NULL, 
    help="Directory containing all the raw count data tables (one per sample)", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
    help="A gtf or gff file allowing to calculate the genes length", metavar="character"),
  make_option(c("-b", "--gene_attribute"), type="character", default="gene", 
    help="The name of the attribute in the gtf referring to the gene information [default= %default]", metavar="character"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
    help="Give the time value of each file by respecting the same order in the vector than the files in the folder. Give a list of type 'time1,time1,time1,time2,time2,time2,time3'. If several files correspond to a same time (replicates), give the same time value and then the script performs the mean on the normalized counts of all the samples of a same time", metavar="character"),
  make_option(c("-n", "--nb_clusters"), type="integer", default=4, 
    help="Number of clusters to generate with Mfuzz (empirical choice) [default= %default]", metavar="integer"),
  make_option(c("-m", "--membership_cutoff"), type="integer", default=0.7, 
    help="The membership cut-off to use with Mfuzz -> see the Mfuzz paper : http://www.bioinformation.net/002/000200022007.pdf [default= %default]", metavar="integer"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
    help="The directory where store the results. By default, the current directory. [default= %default]", metavar="character")
); 
 
# parsing of the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test the three essential arguments and exit if one of them in not given
if (is.null(opt$folder)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input folder).n", call.=FALSE)
}
if (is.null(opt$annotation)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (GFF/GTF file).n", call.=FALSE)
}
if (is.null(opt$time)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (time values vector).n", call.=FALSE)
}

count_files_folder=opt$folder
annotation=opt$annotation
gene_name_attribute=opt$gene_attribute
time=as.vector(strsplit(opt$time,","))[[1]]
nb_clusters=opt$nb_clusters
membership_cutoff=opt$membership_cutoff
output=opt$output

if (is.null(output)){
  output=getwd()
}

annotation_ext=file_ext(annotation)
# create the object with all count files
files=list.files(count_files_folder)
raw=readDGE(files, path=count_files_folder, group=c(1:length(files)), columns=c(1,2))
# create transcripts database from gtf or gff
txdb=makeTxDbFromGFF(annotation,format=annotation_ext)
# then collect the exons per gene id
exons.list.per.gene=exonsBy(txdb,by=gene_name_attribute)
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes=as.data.frame(sum(width(reduce(exons.list.per.gene))))
colnames(exonic.gene.sizes)="gene_length_bp"
# our raw data table containing all the samples
datafile=raw$counts
# remove all the genes with a 0 count in all samples
data = datafile[apply(datafile,1,sum)!=0,]
# determine number of studied samples
nblib= dim(data)[2]
# create a factice vector for DESeq normalization
conds = factor(1:nblib)
# normalize raw read counts by library size with DESeq method
cds = newCountDataSet(data, conds)
cds = estimateSizeFactors(cds)
datanorm = t(t(data)/sizeFactors(cds))
colnames(datanorm)=paste(colnames(data),"normalized_by_DESeq", sep="_")
# merge the raw read counts and the read counts normalized by DESeq
alldata = merge(data, datanorm, by="row.names", all=T)
alldata2=alldata[,-1]
rownames(alldata2)=alldata[,1]
alldata=alldata2
# merge the raw and normalized read counts with the genes length (in bp)
alldata_ok = merge(alldata, exonic.gene.sizes, by="row.names", all.x=T)
alldata=alldata_ok[,-1]
rownames(alldata)=alldata_ok[,1]
# recover of the start of normalized read count columns and the end
start=length(files)+1
end=length(files)*2
# calcul of the RPKN -> normalized read counts / (gene length / 1000)
# gene length is divided by 1000 to transform bp on kbp
rpkn = alldata[,start:end] / (as.vector(dim(alldata)[2]/1000 ))
colnames(rpkn)=paste(colnames(data),"normalized_by_DESeq_and_divided_by_gene_length", sep="_")
# merge of rpkn with the table containing raw and normalized read counts and gene length
alldata_ok = merge(alldata, rpkn, by="row.names", all.x=T)
alldata=alldata_ok[,-1]
rownames(alldata)=alldata_ok[,1]
# write of this table containing the raw read counts and the different normalization
write.table(alldata, paste(output,"normalized_counts_all_genes.txt",sep="/"), sep="\t", quote=F, row.names=T, dec=".")
first_rpkn_column=dim(alldata)[2]-length(files)+1
# here we create a matrix containing all the RPKN columns, used by Mfuzz for the clustering
exprs=as.matrix(alldata[,first_rpkn_column:dim(alldata)[2]])
# and for each time value containing replicates, we calculate the RPKN count means
# if there are no replicates, we keep the initial RPKN
count=1
for ( i in unique(time) ){
  if ( dim(as.data.frame(exprs[,which(time==i)]))[2] == 1 ){
    mean_rpkn=data.frame(exprs[,which(time==i)])
  } else {
    mean_rpkn=data.frame(rowMeans(exprs[,which(time==i)]))
  }
  colnames(mean_rpkn)=i
  if (count == 1){
    mean_rpkn_ok=mean_rpkn
  } else {
    mean_rpkn_ok=merge(mean_rpkn_ok,mean_rpkn,by="row.names")
    rownames(mean_rpkn_ok)=mean_rpkn_ok[,1]
    mean_rpkn_ok=mean_rpkn_ok[,-1]
  }
  count=count+1
}
# here we have a RPKN matrix containing one column per time value (and not one column per sample)
exprs_with_time=as.matrix(mean_rpkn_ok, header=TRUE, sep="\t",row.names=1,as.is=TRUE)
# we create the Mfuzz object (ExpressionSet)
exprSet=ExpressionSet(assayData=exprs_with_time)

#--------------------------------------------------------------------------------------------------------------------------
# As a first step,  we exclude genes with more than 25% of the measurements missing 
# -> genes with 0 RPKN in 25% of the conditions  
exprSet.r=filter.NA(exprSet, thres=0.25)
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Fuzzy  c-means  like  many  other  cluster  algorithms,  does  not  allow  for  missing  values.
# Thus,  we  timelace  remaining  missing  values  by  the  average  values  expression  value  of  the
# corresponding gene.
# Methods for replacement of missing values. Missing values should be indicated by NA in the expression matrix.
# Mode method for replacement of missing values:
  # mean- missing values will be replaced by the mean expression value of the gene,
  # median- missing values will be replaced by the median expression value of the gene,
  # knn- missing values will be replaced by the averging over the corresponding expression values of the k-nearest neighbours,
  # knnw-same replacement method as knn, but the expression values averaged are weighted by the distance to the corresponding neighbour
exprSet.f=fill.NA(exprSet.r,mode="mean")
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Most cluster analyses published include a filtering step to remove genes which are expressed
# at  low  levels  or  show  only  small  changes  in  expression.   Different  filtering  procedures  have
# been proposed.  A popular procedure is the setting of a minimum threshold for variation
# This function can be used to exclude genes with low standard deviation.
# Thus, the value of a filtering threshold remains arbitrary.  As no stringent filtering proce-
# dure currently exists, we avoided any prior filtering of gene data.  This prevents the loss of
# genes that may be biologically important.
tmp=filter.std(exprSet.f,min.std=0, visu=FALSE)
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Since  the  clustering  is  performed  in  Euclidian  space,  the  expression  values  of  genes  were
# standardised to have a mean value of zero and a standard deviation of one.  This step ensures
# that vectors of genes with similar changes in expression are close in Euclidean space
# Importantly, Mfuzz assumes that the given expression data are fully preprocessed including  any  data  normalisation.
# The  function standardise does  not  replace  the  normalisation step (eg RPKN normalization).
exprSet.s=standardise(tmp)
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# clustering
m1=mestimate(exprSet.s)
cl=mfuzz(exprSet.s,c=nb_clusters,m=m1)
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# membership cut-off part and plot clusters
pdf(paste(output,paste(paste("clusters_Mfuzz_membership_equals_",membership_cutoff,sep=""),".pdf",sep=""), sep="/"))
mfuzz.plot2(exprSet.s,cl=cl,time.labels=unique(time),min.mem=membership_cutoff, colo="fancy", x11=FALSE)
dev.off()
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# generates one genes list per cluster
acore.list=acore(exprSet.s,cl=cl,min.acore=membership_cutoff)
for (cluster in 1:nb_clusters){
  print(paste(paste("Number of genes in cluster", cluster, sep=" "),dim(acore.list[[cluster]])[1], sep=" : "))
  cluster_table=merge(alldata,acore.list[[cluster]][2], by="row.names", all.y=TRUE)
  write.table(cluster_table,paste(output,paste(paste("list_of_genes_in_cluster",cluster,sep="_"),".txt"),sep="/"))
}
#--------------------------------------------------------------------------------------------------------------------------
