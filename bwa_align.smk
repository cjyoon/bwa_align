
configfile: 'sample_config.yaml'
configfile: srcdir('path_config.yaml')

SAMBLASTER = config['SAMBLASTER']
BWA = config['BWA']
SAMTOOLS = config['SAMTOOLS']
SAMBAMBA = config['SAMBAMBA']
JAVA = config['JAVA']
GATK = config['GATK']

if config['assembly']=='GRCh38':
    chromosomes = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y']]
    REFERENCE = config['GRCh38_REFERENCE']
    KNOWNINDEL = config['GRCh38_KNOWNINDEL']

elif config['assembly']=='GRCh37':
    chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y']]
    REFERENCE = config['GRCh37_REFERENCE']
    KNOWNINDEL = config['GRCh37_KNOWNINDEL']

else:
    print(config['assembly'] + ' not supported, exiting...')
    sys.exit(1)


rule all:
    input:
        expand("dna_bam/{sample}.fmarked.realigned.bam", sample=config['sample'])

rule bwa_align:
    params:
        name = lambda wildcards: config['sample'][wildcards.sample]['name'], 
        # rg = "@RG\tID:{sample}\tLB:{sample}\tSM:{sample}\tPL:ILLUMINA", 
    input:
        fq1 = lambda wildcards: config['sample'][wildcards.sample]['fq1'], 
        fq2 = lambda wildcards: config['sample'][wildcards.sample]['fq2']
    output:
        bam = temp("tmp_bam/{sample}.fmarked.bam"), 
    log:
        "logs/{sample}.bwa.log"
    threads: 10
    shell:
        "({BWA} mem -M -t {threads} -R '@RG\tID:{params.name}\tLB:{params.name}\tSM:{params.name}\tPL:ILLUMINA' {REFERENCE} {input.fq1} {input.fq2} |"
        "{SAMBLASTER}  -M --addMateTags | "
        "{SAMTOOLS} sort -@ {threads} -O BAM -o {output.bam} - ; {SAMBAMBA} index -t {threads} {output.bam}) "
        "&> {log}"



rule indel_targetcreator:
    input:
        bam = "tmp_bam/{sample}.fmarked.bam", 
    output:
        intervals = 'realign_intervals/{sample}.intervals'
    threads: 1
    log:
        'logs/{sample}.indel_targetcreator.log'
    shell:
        "({JAVA} -Djava.io.tmpdir=tmp/  -jar {GATK} "
        "-T RealignerTargetCreator "
        "-R {REFERENCE} "
        "-known {KNOWNINDEL} "
        "-I {input.bam} "
        "-o {output.intervals}) &> {log}"

rule realign:
    input:
        bam = "tmp_bam/{sample}.fmarked.bam", 
        intervals = 'realign_intervals/{sample}.intervals'
    output:
        bam = "dna_bam/{sample}.fmarked.realigned.bam"
    threads: 1
    log:
        "logs/{sample}.realign.log"
    shell:
        "({JAVA} -Xmx8G -Djava.io.tmpdir=tmp/ -jar {GATK} "
        "-T IndelRealigner "
        "-R {REFERENCE} "
        "-targetIntervals {input.intervals} "
        "-known {KNOWNINDEL} " 
        "-I {input.bam} "
        "-o {output.bam}; "
        "{SAMBAMBA} index -t {threads} {output.bam}) &> {log}"
