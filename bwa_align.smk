
configfile: 'sample_config.yaml'

REFERENCE = '/home/users/cjyoon/Projects/mosaic/reference/Homo_sapiens_assembly38.fasta'
SAMBLASTER = '/usr/local/bin/samblaster'
BWA = '/home/users/cjyoon/tools/bwa-0.7.12/bwa'
SAMTOOLS = '/home/users/cjyoon/anaconda3/bin/samtools'
SAMBAMBA = '/home/users/cjyoon/anaconda3/bin/sambamba'
rule all:
	input:
		expand("temp_dna/{sample}.fmarked.bam", sample=config['sample'])

rule bwa_align:
    params:
        name = lambda wildcards: config['sample'][wildcards.sample]['name'], 
        # rg = "@RG\tID:{sample}\tLB:{sample}\tSM:{sample}\tPL:ILLUMINA", 
    input:
        fq1 = lambda wildcards: config['sample'][wildcards.sample]['fq1'], 
        fq2 = lambda wildcards: config['sample'][wildcards.sample]['fq2']
    output:
        bam = "temp_dna/{sample}.fmarked.bam", 
    log:
        "logs/{sample}.bwa.log"
    threads: 12
    shell:
        "({BWA} mem -M -t {threads} -R '@RG\tID:{params.name}\tLB:{params.name}\tSM:{params.name}\tPL:ILLUMINA' {REFERENCE} {input.fq1} {input.fq2} |"
        "{SAMBLASTER}  -M --addMateTags | "
        "{SAMTOOLS} sort -@ {threads} -O BAM -o {output.bam} - ; {SAMBAMBA} index -t {threads} {output.bam}) "
        "&> {log}"

