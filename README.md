# BWA aligning snakefile
Will align paired-end FASTQ files with `BWA mem` followed by mark duplicate with `samblaster`. Outputs a sorted bam. 

This is equivalent to the following bash script. 
```bash
bwa mem -M -t {threads} -R '@RG\tID:{rgname}\tLB:{rgname}\tSM:{rgname}\tPL:ILLUMINA' {REFERENCE} {fq1} {fq2} | samblaster -M --addMateTags | samtools sort -@ {threads} -O BAM -o {output_bam} -
```

# Configuring `sample_config.yaml`
```yaml
'sample':
    'T2664-Hiseq':
        'name': 'T2664-blood-wgs-ILLUMINA'
        'fq1': ['/path/to/T2664-Hiseq_0_R1.fastq.gz', '/path/to/T2664-Hiseq_1_R1.fastq.gz']
        'fq2': ['/path/to/T2664-Hiseq_0_R2.fastq.gz', '/path/to/T2664-Hiseq_1_R2.fastq.gz']
    'T2664-Novaseq':
        'name': 'T2664-blood-wgs-ILLUMINA'
        'fq1': ['/path/to/T2664-Novaseq_0_R1.fastq.gz']
        'fq2': ['/path/to/T2664-Novaseq_0_R2.fastq.gz']

```
In the above example, running this snakefile will give you two bam files, `dna_bam/T2664-Hiseq.bam` and `dna_bam/T2664-Novaseq.bam`. Both will have their `name` set as their `@RG SM:` read names. In the case of having more than one pairs of FASTQ files for a single sample (e.g. T2664-Hiseq), each pairs will be aligned, mark duplicated, and sorted independently and then will be merged. Once a merged bam file is made, then `GATK IndelRealigner` will realign each bam files around potential indel sites to produce the final bam file. 

# To run
Type in the following command in a directory with `sample_config.yaml` file. 
```
snakemake -p -s ~/scripts/bwa_align/highdepth_bwa_align.smk
```
