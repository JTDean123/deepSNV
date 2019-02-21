# deepSNV
Single nucleotide variant calling in NGS with deep learning

deepSNV is a deep learning framework for calling SNVs in NGS data.  Technical details can be found in the deepSNV.pdf in this repository.  It can be run as either a model building and validation tool or used only to make predictions using the included pre-trained convoluted neural network or one of your choosing.  
  
## Usage  
An example of building and training a CNN and making SNV predictions is as follows:

```
$ python3 deepsnv.py --sample_num 10 \
    --vcf_path my_vcf.vcf.gz \
    --genome_path my_refererence_genome.fa \
    --bam_path my_indexed_sorted_bam.bam \
    --len_path chromosome_sizes.txt \
    --epochs 50
```

The options are as follows:  
`--vfc_path`:  path to vcf file containing known SNVs   
`--sample_num`: number of SNVs for training, randomly selected from vcf  
`--len_path`:  path to text file containing chromosome size info, see `chromosome_sizes.txt` for example  
`--genome_path`:  path to indexed reference genome file. This can be created for hg38, for example, via:  

```
$ wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz  
$ gunzip hg38.fa.gz   
$ bwa index hg38.fa
```   

```--bam_path```:  path to indexed, sorted bam file
```--len_path```:  path to chromosome size tsv file (see hg19.sizes for example)
```--epochs```:  number of epochs to train model   
  
Output from the above command is a trained convoluted neural netowork and confusion matrix for held out testing data.  Alternatively, predictions can be made using a pre-trained model as follows:  

```
$ python3 predictCNN.py --genome_path bams/wg.fa \
    --bam_path my_indexed_sorted_bam.bam \
    --preds_path predict.csv \
    --model_path deepSNV.h5

```  

Options are as follows:   
```--bam_path```:  path to indexed, sorted bam file to analyze  
```--preds_path```:  path to a csv containing chromosome and coordinate information for positions to evalate for SNV.  (see predict.csv for example)  
```--model_path```:  path to pre-trained keras CNN model   

Output from the above command is a csv of input SNV candidate coordinates and a label of 1 (SNV) or 0 (non-SNV).   


## Get deepSNV  
  
```$ git clone https://github.com/JTDean123/deepSNV.git```

## Requirements  
  
python3, pandas, numpy, sklean, keras (tensorflow backend), pysam
  
## Other  
  
This is a skeleton, bare bones first attempt at this idea.  I am sure there is lots of room for improvement, and I'd love to hear from you.




