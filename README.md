# deepSNV
Single nucleotide variant calling in NGS with deep learning

deepSNV is a deep learning framework for calling SNVs in NGS data.  It can be run as either a model building and validation tool or used only to make predictions using the included pre-trained convoluted neural network or one of your choosing.  An example of building and training a CNN and making SNV predictions is as follows:

```
python3 deepsnv.py --sample_num 10 \
    --vcf_path my_vcf.vcf.gz \
    --genome_path my_refererence_genome.fa \
    --bam_path my_indexed_sorted_bam.sorted.bam \
    --len_path chromosome_sizes.txt \
    --epochs 50
```

The options are as follows:  
```--vfc_path```:  path to vcf file containing known SNVs   
```--genome_path```:  path to indexed reference genome file. This can be created for hg38, for example, via:
```
$ wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz  
$ gunzip hg38.fa.gz   
$ bwa index hg38.fa
```   

```--bam_path```:  path to indexed, sorted bam file
```--len_path```:  path to chromosome size tsv file (see hg19.sizes for example)
```--epochs```:  number of epochs to train model   
  
Output from the above command is a trained convoluted neural netowork and confusion matrix for held out testing data.




