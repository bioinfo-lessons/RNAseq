
# RNAseq hands-on session

## 1.Download Raw Data fro SRA:

### GEO DATASET

Links to the samples used in the article: GSE81214 (GEO id)

- GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81214
- SRA link: https://www.ncbi.nlm.nih.gov/sra?term=SRP074494

### Install sra-tools from https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.current-ubuntu64.tar.gz
export PATH=/YOUR_PATH/sratoolkit.3.0.1-ubuntu64/bin/:$PATH
which fastq-dump
```

### Download accession list and use it in a loop for downloading data and tranform it to fastq:

for SRAid in $(cat SRR_Acc_List.txt); do echo $SRAid; fastq-dump $SRAid; done

The following files are produced:
* SRR3480371.fastq (CONTROL 1): 42059631 reads
* SRR3480372.fastq (CONTROL 2): 54401313 reads
* SRR3480373.fastq (TREATMENT 1): 52991800 reads
* SRR3480374.fastq (TREATMENT 2): 64635061 reads


## 2.Create a conda environment with all the required tools:

```
conda create --name RNAseq python=3.7
conda activate RNAseq
```

### A quality control tool for high throughput sequence data:
```
conda install -c bioconda fastqc
```

### Graph-based alignment of next generation sequencing reads to a population of genomes:
```
conda install -c bioconda hisat2
 ```
### Tools for dealing with SAM, BAM and CRAM files:
```
conda install -c bioconda samtools
```
### HTSeq is a Python library to facilitate processing and analysis of data from high-throughput sequencing (HTS) experiments:
```
conda install -c bioconda htseq
```
### A set of user-friendly tools for normalization and visualzation of deep-sequencing data:
```
conda install -c bioconda deeptools
```
### Integrative Genomics Viewer. Fast, efficient, scalable visualization tool for genomics data and annotations:
```
conda install -c bioconda igv
```
### Check that you have everything set:
```
conda list
```

**IN ORDER TO RUN THE EXERCISE IN CLASS WE ARE GOING TO USE A SUBSET OF THE INITIAL FASTQS**

## 3.Alignment of reads with HISAT2:

if reference indexing was not previously done:
```
hisat2-build --seed 123 -p 2 REF/chr19_GRCh38.fa REF/chr19_GRCh38 
```

**read alignment**

```
mkdir -p hisat2
hisat2 --new-summary --summary-file hisat2/SRR3480371.hisat2.summary  --rna-strandness R --seed 123 --phred33 -p 2 -k 1 -x REF/chr19_GRCh38 -U SRR3480371.subset.fastq -S hisat2/SRR3480371.sam

hisat2 --new-summary --summary-file hisat2/SRR3480372.hisat2.summary  --rna-strandness R --seed 123 --phred33 -p 2 -k 1 -x REF/chr19_GRCh38 -U SRR3480372.subset.fastq -S hisat2/SRR3480372.sam

hisat2 --new-summary --summary-file hisat2/SRR3480373.hisat2.summary  --rna-strandness R --seed 123 --phred33 -p 2 -k 1 -x REF/chr19_GRCh38 -U SRR3480373.subset.fastq -S hisat2/SRR3480373.sam

hisat2 --new-summary --summary-file hisat2/SRR3480374.hisat2.summary  --rna-strandness R --seed 123 --phred33 -p 2 -k 1 -x REF/chr19_GRCh38 -U SRR3480374.subset.fastq -S hisat2/SRR3480374.sam
```

## 4.SAMtools: view, sort and index:
```
cd hisat2
samtools view -bS SRR3480371.sam > SRR3480371.bam
samtools view -bS SRR3480372.sam > SRR3480372.bam
samtools view -bS SRR3480373.sam > SRR3480373.bam
samtools view -bS SRR3480374.sam > SRR3480374.bam

samtools sort SRR3480371.bam -o SRR3480371.sorted.bam
samtools sort SRR3480372.bam -o SRR3480372.sorted.bam
samtools sort SRR3480373.bam -o SRR3480373.sorted.bam
samtools sort SRR3480374.bam -o SRR3480374.sorted.bam 

samtools index SRR3480371.sorted.bam
samtools index SRR3480372.sorted.bam
samtools index SRR3480373.sorted.bam
samtools index SRR3480374.sorted.bam
```

## 5.HTseq-count to count reads:
```
cd ..
mkdir -p htseq
htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name hisat2/SRR3480371.sorted.bam REF/subset_gencode.v38.annotation.gtf > htseq/SRR3480371.htseq

htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name hisat2/SRR3480372.sorted.bam REF/subset_gencode.v38.annotation.gtf > htseq/SRR3480372.htseq

htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name hisat2/SRR3480373.sorted.bam REF/subset_gencode.v38.annotation.gtf > htseq/SRR3480373.htseq

htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name hisat2/SRR3480374.sorted.bam REF/subset_gencode.v38.annotation.gtf > htseq/SRR3480374.htseq
```
## 6. BigWig files for the Genome Browser
```
bamCoverage -b hisat2/SRR3480371.sorted.bam -o control_1.bw --normalizeUsing BPM
bamCoverage -b hisat2/SRR3480372.sorted.bam -o control_2.bw --normalizeUsing BPM
bamCoverage -b hisat2/SRR3480373.sorted.bam -o treatment_1.bw --normalizeUsing BPM
bamCoverage -b hisat2/SRR3480374.sorted.bam -o treatment_2.bw --normalizeUsing BPM
```
## 7. Open Igv load GRCh38 genome and search information for the following genes

* UPF1
* LPAR2
* SIRT2
* TPKC
* EML2
* JUND
