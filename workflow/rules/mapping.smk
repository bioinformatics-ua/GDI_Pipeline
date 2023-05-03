
rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}_{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    threads: 20,
    log:
        "logs/trimmomatic/{sample}_{unit}.log",
    wrapper:
        "v1.23.5/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
#        r1=temp("results/trimmed/{sample}_{unit}.1.fastq.gz"),
#        r2=temp("results/trimmed/{sample}_{unit}.2.fastq.gz"),
#        r1_unpaired=temp("results/trimmed/{sample}_{unit}.1.unpaired.fastq.gz"),
#        r2_unpaired=temp("results/trimmed/{sample}_{unit}.2.unpaired.fastq.gz"),
        r1="results/trimmed/{sample}_{unit}.1.fastq.gz",
        r2="results/trimmed/{sample}_{unit}.2.fastq.gz",
        r1_unpaired="results/trimmed/{sample}_{unit}.1.unpaired.fastq.gz",
        r2_unpaired="results/trimmed/{sample}_{unit}.2.unpaired.fastq.gz",
        ## create a huge file, avoid this
        #trimlog="results/trimmed/{sample}_{unit}.trimlog.txt",
    threads: 10,
    params:
        **config["params"]["trimmomatic"]["pe"],
        #extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    resources:
        mem_mb=2048
    log:
        "logs/trimmomatic/{sample}_{unit}.log",
    wrapper:
        "v1.23.5/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        "results/mapped/{sample}_{unit}.sorted.bam",
    log:
        "logs/bwa_mem/{sample}_{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
        sort_extra="",
    threads: 6,
    wrapper:
        "v1.23.5/bio/bwa/mem"

#rule sort_bam_mapped_reads:
#    input:
#        "results/mapped/{sample}_{unit}.temp.bam",
#    output:
#        "results/mapped/{sample}_{unit}.sorted.bam",
#    log:
#        "logs/samtools/sort/{sample}_{unit}.log",
#    params:
#        extra="-m 2G",
#    threads: 2,
#    wrapper:
#        "v1.23.5/bio/samtools/sort"


rule mark_duplicates:
    input:
        bams="results/mapped/{sample}_{unit}.sorted.bam",
    output:
        bam="results/dedup/{sample}_{unit}.bam",
        metrics="results/qc/dedup/{sample}_{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}_{unit}.log",
    params:
        ## remove duplciates
        config["params"]["picard"]["MarkDuplicates"],
    resources:
        mem_mb=50000,
    threads: 12,
    wrapper:
        "0.74.0/bio/picard/markduplicates"
#        "v1.23.5/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        known_idx="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table="results/recal/{sample}_{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}_{unit}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        #spark_runner="",  # optional, local by default
    resources:
        mem_mb=50000,
    threads: 1,
    wrapper:
#        "0.74.0/bio/gatk/baserecalibrator"
        "v1.23.5/bio/gatk/baserecalibratorspark"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}_{unit}.grp",
    output:
        bam=protected("results/recal/{sample}_{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}_{unit}.log",
    params:
        extra=get_regions_param(),
        #spark_runner="",  # optional, local by default
    resources:
        mem_mb=50000,
    threads: 1,
    wrapper:
#        "v1.23.5/bio/gatk/applybqsr"
        "v1.23.5/bio/gatk/applybqsrspark"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "v1.23.5/bio/samtools/index"
