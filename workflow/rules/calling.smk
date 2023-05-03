if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("results/called/{sample}.{contig}.g.vcf.gz"),
        gvcf_tbi=protected("results/called/{sample}.{contig}.g.vcf.gz.tbi"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
    resources:
        mem_mb=50000,
    threads: 1
    wrapper:
        "v1.23.5/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
        gvcfs_tbi=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz.tbi", sample=samples.index
        ),
    output:
        gvcf=temp("results/called/all.{contig}.g.vcf.gz"),
        gvcf_tbi=temp("results/called/all.{contig}.g.vcf.gz.tbi"),
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    resources:
        mem_mb=50000,
    wrapper:
        "v1.23.5/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="results/called/all.{contig}.g.vcf.gz",
        gvcf_tbi="results/called/all.{contig}.g.vcf.gz.tbi",
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
        vcf_tbi=temp("results/genotyped/all.{contig}.vcf.gz.tbi"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    resources:
        mem_mb=50000,
    threads: 1
    wrapper:
        "v1.23.5/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
        vcfs_tbi=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz.tbi", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.gatk.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    resources:
        mem_mb=50000,
    threads: 6
    wrapper:
        "v1.23.5/bio/picard/mergevcfs"
