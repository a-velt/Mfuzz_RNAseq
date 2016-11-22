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
# Usage :
# Complete command : /usr/bin/Rscript Mfuzz_RNAseq.R -f count_files_folder -a annotation -b gene_name_attribute -t time -n nb_clusters -m membership_cutoff -s 0 -e 0.25 -r "mean" -o output_directory
# Minimal command : /usr/bin/Rscript Mfuzz_RNAseq.R -f count_files_folder -a annotation -t time
#
# Arguments description :
# count_files_folder -> directory containing all the raw count data tables (one per sample)
# annotation -> a gtf or gff file with transcripts/genes information allowing to calculate the genes length (sum of the exons length, overlap of exons is take into account)
# gene_name_attribute -> the name of the attribute in the gtf referring to the gene information
# time -> give the time value of each file by respecting the same order in the vector than the files in the folder.
#   if several files correspond to a same time (replicates), give the same time value and then the script performs the mean on the normalized counts of all the
#   samples of a same time
# nb_clusters -> number of clusters to generate with Mfuzz (empirical choice)
# membership_cutoff -> the membership cut-off to use with Mfuzz -> see the Mfuzz paper : http://www.bioinformation.net/002/000200022007.pdf
#                       you can give one value, eg "0.7" or several values separated by ",", eg '0.5,0.7'
# output -> directory where store the results
#
# Examples of arguments to give :
# count_files_folder="/home/user/count_files_folder"
# annotation="HS.Genes.v2.gff"
# gene_name_attribute="gene"
# time="time1,time1,time1,time2,time2,time2,time3"
# nb_clusters = 4
# membership_cutoff = '0.5,0.7'
# min_std_threshold= 0
# exclude_thres= 0.25
# replacement_mode="mean"
# output="/home/user/cluster_output"
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
    help="[REQUIRED] Directory containing all the raw count data tables (one per sample)", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
    help="[REQUIRED] A gtf or gff file with transcripts/genes information allowing to calculate the genes length (sum of the exons length, overlap of exons is take into account)", metavar="character"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
    help="[REQUIRED] Give the time value of each file by respecting the same order in the vector than the files in the folder. 
            Give a list of type 'time1,time1,time1,time2,time2,time2,time3'. 
            If several files correspond to a same time (replicates), give the same time value and then the script performs the mean on the normalized counts of all the samples of a same time", metavar="character"),
  make_option(c("-b", "--gene_attribute"), type="character", default="gene", 
    help="The name of the attribute in the gtf referring to the gene information [default= %default]", metavar="character"),
  make_option(c("-n", "--nb_clusters"), type="integer", default=as.numeric(4), 
    help="Number of clusters to generate with Mfuzz (empirical choice) [default= %default]", metavar="integer"),
  make_option(c("-m", "--membership_cutoff"), type="character", default=as.numeric(0.7), 
    help="The membership cut-off to use to generate gene lists for each cluster with Mfuzz.
          By default, genes having a membership value of 0.7 for the cluster are recovered in the list for this cluster.
          You can give one value, eg '0.7' or several values separated by ",", eg '0.5,0.7'.
          See the Mfuzz paper : http://www.bioinformation.net/002/000200022007.pdf [default= %default]", metavar="character"),
  make_option(c("-s", "--min_std"), type="double", default=as.numeric(0), 
    help="Threshold for minimum standard deviation, use by Mfuzz. If the standard deviation of a gene's expression is smaller than min.std the corresponding gene will be excluded.
          Default : no filtering. [default= %default]", metavar="double"),
  make_option(c("-e", "--exclude_thres"), type="double", default=as.numeric(0.25), 
    help="Exclude genes with more than n% of the measurements missing [default= %default] -> by default, genes with 25% of the measurements missing are excluded.", metavar="double"),
  make_option(c("-r", "--replacement_mode"), type="character", default="mean", 
    help="Mode method for replacement of missing values. Fuzzy  c-means  like  many  other  cluster  algorithms,  does  not  allow  for  missing  values.
            Thus, by default, we  timelace  remaining  missing  values  by  the  average  values  expression  value  of  the corresponding gene. [default= %default]
            Other available methods : median, knn, knnw", metavar="character"),
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

# variable assignment
count_files_folder=opt$folder
annotation=opt$annotation
gene_name_attribute=opt$gene_attribute
time=as.vector(strsplit(opt$time,","))[[1]]
nb_clusters=opt$nb_clusters
membership_cutoff=as.vector(strsplit(opt$membership_cutoff,","))[[1]]
min_std_threshold=opt$min_std
exclude_thres=opt$exclude_thres
replacement_mode=opt$replacement_mode
output=opt$output

# test of the output, if empty, give the current directory
if (is.null(output)){
  output=getwd()
}
# create output if doesn't exists
dir.create(output, showWarnings = FALSE)

###########################################################################################################################
# Normalization part
###########################################################################################################################

# determine the extension of the annotation file (may be gtf or gff)
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
alldata_tmp=alldata[,-1]
rownames(alldata_tmp)=alldata[,1]
alldata=alldata_tmp
# merge the raw and normalized read counts with the genes length (in bp)
alldata_tmp = merge(alldata, exonic.gene.sizes, by="row.names", all.x=T)
alldata=alldata_tmp[,-1]
rownames(alldata)=alldata_tmp[,1]
# recover of the start of normalized read count columns and the end
start=length(files)+1
end=length(files)*2
# calcul of the RPKN -> normalized read counts / (gene length / 1000)
# gene length is divided by 1000 to transform bp on kbp
rpkn = alldata[,start:end] / (as.vector(dim(alldata)[2]/1000 ))
colnames(rpkn)=paste(colnames(data),"normalized_by_DESeq_and_divided_by_gene_length", sep="_")
# merge of rpkn with the table containing raw and normalized read counts and gene length
alldata_tmp = merge(alldata, rpkn, by="row.names", all.x=T)
alldata=alldata_tmp[,-1]
rownames(alldata)=alldata_tmp[,1]
# write of this table containing the raw read counts and the different normalization
write.table(alldata, paste(output,"normalized_counts_all_genes.txt",sep="/"), sep="\t", quote=F, row.names=T, dec=".")

###########################################################################################################################
# Mfuzz part
###########################################################################################################################

# determine the first RPKN column in the alldata object -> RPKN are used by Mfuzz for genes clustering
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
exprSet.r=filter.NA(exprSet, thres=exclude_thres)
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
exprSet.f=fill.NA(exprSet.r,mode=replacement_mode)
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# As soft clustering is noise robust, pre-filtering can usually be avoided. 
# However, if the number of genes with small expression changes is large, such pre-filtering may be necessary to reduce noise. 
# This function can be used to exclude genes with low standard deviation.
# min.std : threshold for minimum standard deviation. 
# If the standard deviation of a gene's expression is smaller than min.std the corresponding gene will be excluded.
tmp=filter.std(exprSet.f,min.std=min_std_threshold, visu=FALSE)
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

for (membership in membership_cutoff){
  membership=as.numeric(membership)
  # create one output folder per membership
  dir=paste(output,paste("cluster_with_membership",membership, sep="_"),sep="/")
  dir.create(dir, showWarnings = FALSE)
  #--------------------------------------------------------------------------------------------------------------------------
  # membership cut-off part and plot clusters
  pdf(paste(dir,paste(paste("clusters_Mfuzz_membership_equals_",membership,sep=""),".pdf",sep=""), sep="/"))
  mfuzz.plot2(exprSet.s,cl=cl,time.labels=unique(time),min.mem=membership, colo="fancy", x11=FALSE)
  dev.off()
  #--------------------------------------------------------------------------------------------------------------------------

  #--------------------------------------------------------------------------------------------------------------------------
  # generates one genes list per cluster
  acore.list=acore(exprSet.s,cl=cl,min.acore=membership)
  for (cluster in 1:nb_clusters){
    print(paste(paste("Number of genes in cluster", cluster, sep=" "),dim(acore.list[[cluster]])[1], sep=" : "))
    cluster_table=merge(alldata,acore.list[[cluster]][2], by="row.names", all.y=TRUE)
    write.table(cluster_table,paste(dir,paste(paste("list_of_genes_in_cluster",cluster,sep="_"),".txt"),sep="/"), sep="\t",row.names=F, dec=".")
  }
  #--------------------------------------------------------------------------------------------------------------------------
}
