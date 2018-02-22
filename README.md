# deepSNV
Single nucleotide variant calling in NGS with deep learning

deepSNV is a deep learning framework for calling SNVs in NGS data.  It can be run as either a model building and predictive tool or used only to make predictions using the included pre-trained convoluted neural network.  An example of building and training a CNN and making SNV predictions is as follows:

```
python3 deepsnv.py --sample_num 10 \
    --vcf_path my_vcf.vcf.gz \
    --genome_path my_refererence_genome.fa \
    --bam_path my_indexed_sorted_bam.sorted.bam \
    --len_path chromosome_sizes.txt \
    --epochs 50
```
