"""Script to streamline high depth WGS sequence from FASTQ files to merged BAM files
2019.05.19 cjyoon
"""
import subprocess
import shlex

configfile: 'sample_config.yaml'

REFERENCE = '/home/users/cjyoon/Projects/mosaic/reference/Homo_sapiens_assembly38.fasta'
SAMBLASTER = '/usr/local/bin/samblaster' # 0.1.24
BWA = '/home/users/cjyoon/tools/bwa-0.7.12/bwa' # 0.7.12-r1039
SAMTOOLS = '/home/users/cjyoon/anaconda3/bin/samtools' # 1.9 
SAMBAMBA = '/home/users/cjyoon/anaconda3/bin/sambamba' # 0.6.6
JAVA = '/home/users/cjyoon/anaconda3/bin/java'
GATK = '/home/users/cjyoon/tools/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar'
KNOWNINDEL = '/home/users/cjyoon/reference/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf'
    
rule all:
    input:
        expand("dna_bam/{sample}.fmarked.realigned.bam", sample=config['sample'])

def bwa_align(fq1, fq2, reference, rgname, output_bam, threads, log):
    """function to run bwa alignment
    equivalent to the shell command: 
    bwa mem -M -t {threads} -R '@RG\tID:{rgname}\tLB:{rgname}\tSM:{rgname}\tPL:ILLUMINA' {REFERENCE} {fq1} {fq2} | samblaster -M --addMateTags | samtools sort -@ {threads} -O BAM -o {output_bam} - """
    bwa_cmd = f"{BWA} mem -M -t {threads} -R '@RG\tID:{rgname}\tLB:{rgname}\tSM:{rgname}\tPL:ILLUMINA' {reference} {fq1} {fq2}"
    bwa_process = subprocess.Popen(shlex.split(bwa_cmd), stdout=subprocess.PIPE)
    print(bwa_cmd)
    samblaster_cmd = f'{SAMBLASTER} -M --addMateTags'
    samblaster_process = subprocess.Popen(shlex.split(samblaster_cmd), stdin=bwa_process.stdout, stdout=subprocess.PIPE)
    bwa_process.stdout.close()
    print(samblaster_cmd)
    sort_cmd = f'{SAMTOOLS} sort -@ {threads} -O BAM -o {output_bam} -'
    sort_process = subprocess.Popen(shlex.split(sort_cmd), stdin=samblaster_process.stdout)
    samblaster_process.stdout.close()
    print(sort_cmd)

    output = sort_process.communicate()[0]

    return 0


rule bwa_align:
    params:
        name = lambda wildcards: config['sample'][wildcards.sample]['name'], 
        # rg = "@RG\tID:{sample}\tLB:{sample}\tSM:{sample}\tPL:ILLUMINA", 
    input:
        fq1 = lambda wildcards: config['sample'][wildcards.sample]['fq1'], 
        fq2 = lambda wildcards: config['sample'][wildcards.sample]['fq2']
    output:
        bam = "tmp_bam/{sample}.fmarked.bam", 
        bai = "tmp_bam/{sample}.fmarked.bam.bai", 
    log:
        "logs/{sample}.bwa.log"
    threads: 12
    run:
        bamfiles = []
        for i, (fq1, fq2) in enumerate(zip(input.fq1, input.fq2)):
            bamfile = re.sub(r'.bam$', f'.{i}.bam', output.bam)
            bamfiles.append(bamfile)
            print(i, fq1, fq2, bamfile)
            bwa_align(fq1, fq2, REFERENCE, params.name, bamfile, threads, log)

        bamfiles_string = ' '.join(bamfiles)
        merge_cmd = f'{SAMTOOLS} merge -@ {threads} {output.bam} {bamfiles_string}'
        print(merge_cmd)
        merge_execute = subprocess.Popen(shlex.split(merge_cmd))
        merge_execute.wait()

        index_cmd = f'{SAMBAMBA} index -t {threads} {output.bam}'
        index_execute = subprocess.Popen(shlex.split(index_cmd))
        index_execute.wait()
        # cleanup_cmd = f'rm -rf {bamfiles_string}'
        # cleanup_execute = subprocess.Popen(shlex.split(cleanup_cmd))
        # cleanup_execute.wait()


rule indel_targetcreator:
    input:
        bam = "tmp_bam/{sample}.fmarked.bam", 
    output:
        intervals = 'realign_intervals/{sample}.intervals'
    threads: 4
    log:
        'logs/{sample}.indel_targetcreator.log'
    shell:
        "({JAVA} -jar {GATK} "
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
    threads: 4
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
