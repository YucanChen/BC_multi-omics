
# DNA methylation (EPIC) data

The DNA methylation data were preprocessed mainly using minfi and ChAMP R packages.

## import the data

```{r}
# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="infinium_methylationEPIC_sample_sheet_CPOS-221123-DS-15432.csv")
targets

# read in the raw intensity data into R using the read.metharray.exp function. This creates an RGChannelSet object that contains all the raw intensity data, from both the red and green colour channels, for each of the samples.
rgSet <- read.metharray.exp(targets=targets,force=TRUE) # read in the raw data from the IDAT files.

# sample metadata information
group_list<-data.frame(Sample_Name=targets$Sample_Name)
group_list$group<-group_list$Sample_Name %>% lapply(.,function(x){strsplit(x,"-")[[1]][4]}) %>% unlist
group_list[which(!(group_list$group %in% c("N","T"))),"group"]<-group_list[which(!(group_list$group %in% c("N","T"))),"Sample_Name"] %>% lapply(.,function(x){strsplit(x,"-")[[1]][3]}) %>% unlist
group_list$batch<-group_list$Sample_Name %>% lapply(.,function(x){targets[which(targets$Sample_Name==x),"Sample_Group"]}) %>% unlist
```

## process the data

```{r}
# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$Sample_Name

# RGChannelSet to MethylSet
MSet <- preprocessRaw(rgSet)
MSet

# MethylSet to RatioSet
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)

# RatioSet to GenomicRatioSet
gset <- mapToGenome(ratioSet)

# The functions getBeta, getM and getCN work on the GenomicRatioSet return respectively the Beta value matrix, M value matrix and a the Copy Number matrix.
beta <- getBeta(gset)
m <- getM(gset)
cn <- getCN(gset)
```

## quality control

### minfi QC

```{r}
qc <- getQC(MSet)
plotQC(qc)

# Control probes plot
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")
qcReport(rgSet, pdf= "E:/PhD_period/breast-cancer_multi-omics/DNA_methyl/qcReport.pdf")

# mds plot
library("wateRmelon")
beta.m=m[rowMeans(m)>0.005,]

pdf(file="densityBeanPlot.pdf")
par(oma=c(2,10,2,2))
densityBeanPlot(beta.m, sampGroups = group_list$group)
dev.off()
pdf(file="mdsPlot.pdf")
mdsPlot(beta.m, numPositions = 1000, sampGroups = group_list$group)
dev.off()
```

### ChAMP QC

```{r}
library(RColorBrewer)
library(dendextend)
champ.QC(beta,pheno=targets$Sample_Group)
champ.QC(beta,pheno=group_list$group)
```

### PCA

```{r}
library("missMDA")
library("FactoMineR")
library("factoextra")

dat=t(beta)
write.csv(dat,"allsample.beta_value.raw.csv",row.names = T,col.names=T)

res.comp <- imputePCA(dat, graph = FALSE) # Missing values are imputed by the mean of the variable: you should use the imputePCA function of the missMDA package
dat.pca <- PCA(t(res.comp$completeObs), graph = FALSE)

fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list$group, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

ggsave('all_samples_PCA.png')
```

## beta value normalization

```{r}
library("ChAMP")
myLoad.impute = champ.impute(beta = beta,pd=targets) # NA value is not allowed in beta matrix.
myNorm <- champ.norm(beta=myLoad.impute$beta,arraytype="EPIC",cores=3,plotBMIQ=FALSE)
dim(myNorm)
pD=myLoad.impute$pd

write.csv(myNorm,"allsample.beta_value.NAimputed_normalized.csv")
```

## differential analysis

```{r}
# Differential Methylation Probes
## ChAMP
myDMP <- champ.DMP(beta = myNorm, pheno = group_list$group)
head(myDMP[[1]])
## minfi
dmp <- dmpFinder(MSet, pheno=group_list, type="categorical")
dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]

# Differential Methylation Blocks
myBlock <- champ.Block(beta=myNorm,pheno=group_list$group,arraytype="EPIC")

# Differential Methylation Regions
myDMR <- champ.DMR(beta=myNorm,pheno=group_list$group,method="Bumphunter",arraytype="EPIC") # Bumphunter detected 4004 DMRs with P value <= 0.05
```


# RNA-seq data

The RNA-seq data were preprocessed via STAR alignment and stringtie transcript assembly and quantification.

## remove adaptors

Two softwares were used for trimming adapters thoroughly.

```bash
module load java
module load cutadapt

for i in $(cat ~/sample_list.RNA.unique.txt)
do
j_1=$(echo $i\_1*)
j_2=$(echo $i\_2*)
outfile=`basename $i`
t_1=~/clean_fastq/`basename $i`_1P.fastq.gz
t_2=~/clean_fastq/`basename $i`_2P.fastq.gz
o_1=~/fastq_clean/`basename $i`_R1.fastq.gz
o_2=~/fastq_clean/`basename $i`_R2.fastq.gz
java -jar ~/biosofts/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 16 $j_1 $j_2 -baseout ~/clean_fastq/$outfile.fastq.gz ILLUMINACLIP:~/biosofts/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:65 && \
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --minimum-length 20 -o $o_1 -p $o_2 $t_1 $t_2
done
```

## fastqc

```bash
for i in $(cat ~/sample_list.RNA.unique.txt)
do
file_1=~/fastq_clean/`basename $i`_R1.fastq.gz
file_2=~/fastq_clean/`basename $i`_R2.fastq.gz
~/biosofts/FastQC/fastqc -o ~/fastqc/ $file_1 -t 16
~/biosofts/FastQC/fastqc -o ~/fastqc/ $file_2 -t 16
done
wait

multiqc ~/fastqc/ -o ~/fastqc/ -i multiqc_clean
```

## STAR alignment

- The reference human genome version hg19/GRCh37 is downloaded from: 
1. fasta file: GRCh37.primary_assembly.genome.fa.gz from gencode Release 45 (GRCh37).
2. gtf file: gencode.v45lift37.basic.annotation.gtf.gz from gencode Release 45 (GRCh37).

- EBV genome: Epstein-Barr virus, strain: B95-8, NCBI RefSeq Complete Genome, link to the gff and fasta file: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/10376/

The fasta files and gtf files from both species were merged for alignment.

```bash
# build index
ref=~/references/hg19_EBV/GRCh37.EBV.genome.fa
gtf=~/references/hg19_EBV/EBV.genomic.gene.deRNA.gtf

~/biosofts/STAR-2.7.11a/bin/Linux_x86_64_static/STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ~/references/hg19_EBV/ --genomeFastaFiles $ref --sjdbGTFfile $gtf --sjdbOverhang 149
wait

# alignment
module load java
for i in $(cat ~/sample_list.RNA.unique.txt)
do
file_1=~/fastq_clean/`basename $i`_R1.fastq.gz
file_2=~/fastq_clean/`basename $i`_R2.fastq.gz
~/biosofts/STAR-2.7.11a/bin/Linux_x86_64_static/STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 10 \
--genomeDir ~/references/hg19_EBV/ --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outFilterMismatchNmax 2 --outSJfilterReads Unique \
--outSAMmultNmax 1 --outFileNamePrefix ~/alignments_STAR/`basename $i` --readFilesCommand gunzip -c --readFilesIn $file_1 $file_2
done
```

## stringtie transcript assembly

```bash
ref_gtf=~/references/hg19_EBV/GRCh37.EBV.genome.gtf
for i in $(cat ~/sample_list.RNA.unique.txt)
do
stringtie -p 10 -G $ref_gtf -o ~/stringtie/`basename $i`.gtf ~/alignments_STAR/`basename $i`Aligned.sortedByCoord.out.bam
done

# merge individual gtf files
gtf_list=~/stringtie/mergelist.txt
stringtie --merge -p 20 -G $ref_gtf -o /groups/cgsd/cyc001/BC_multiomics/RNAseq/stringtie/BC_stringtie_merged.gtf $gtf_list

# gffcompare
gffcompare -r $ref_gtf -G -o ~/stringtie/gffcompare/BC_merged ~/stringtie/BC_stringtie_merged.gtf
```

## stringtie requantification with merged gtf

```bash
ref_gtf=~/stringtie/BC_stringtie_merged.gtf
for i in $(cat ~/sample_list.RNA.unique.txt)
do
stringtie -e -B -p 16 -G $ref_gtf -o ~/requantification/`basename $i`/`basename $i`.gtf ~/alignments_STAR/`basename $i`Aligned.sortedByCoord.out.bam
done

# output gene and transcript count matrices
prepDE=~/biosofts/prepDE.py3
cd ~/requantification/
python $prepDE -i ~/requantification/sample_lst.txt # this will generate gene_count_matrix.csv and transcript_count_matrix.csv (raw data).
```

## screen for known genes only

```{r}
df.raw_read_count<-read.csv("~/requantification/gene_count_matrix.csv")
cts.ref<-df.raw_read_count[grepl("\\|",df.raw_read_count$gene_id),]
cts.ref$gene_id %<>% lapply(.,function(x){strsplit(x,"\\|")[[1]][2]}) %>% unlist
tmp.1<-lapply(colnames(tmp),function(x){gsub("^X","",gsub("\\.","-",x))}) %>% unlist
tmp<-rbind(tmp.1,tmp);rownames(tmp)[1]<-""
write.csv(tmp,"~/gene_count_matrix.known_genes.csv",col.names = F)
```

## qualimap QC

```bash
ref_gtf=~/references/hg19_EBV/GRCh37.EBV.genome.gtf
unset DISPLAY
for i in $(cat ~/sample_list.RNA.unique.txt)
do
qualimap rnaseq \
-outdir ~/qualimap/`basename $i` \
-a proportional \
-bam ~/alignments_STAR/`basename $i`Aligned.sortedByCoord.out.bam \
-gtf $ref_gtf \
--java-mem-size=8G --paired
done
```

## Failure modes during running the scripts

Some factors may result in failure in the process:
1. Version conflicts for the scripts of R and python function in the above tools and your installed R/python version. Adjusting a few command lines according to the error report will solve the problem.
2. When converting the EBV reference gff file into gtf file, the types of gene elements may need to be selected for correct recognision by the STAR tool.
