### split vcf in snp and indel
rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/all.gatk.vcf.gz",
    output:
        vcf=temp("results/filtered/all.gatk.{vartype}.vcf.gz"),
        vcf_tbi=temp("results/filtered/all.gatk.{vartype}.vcf.gz.tbi"),
    params:
        extra=get_vartype_arg,
    log:
        "logs/gatk/selectvariants/gatk.{vartype}.log",
    threads: 1
    resources:
        mem_mb=50000,
    wrapper:
        "v1.23.5/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.gatk.{vartype}.vcf.gz",
        vcf_tbi="results/filtered/all.gatk.{vartype}.vcf.gz.tbi",
    output:
        vcf=temp("results/filtered/all.gatk.{vartype}.hardfiltered.vcf.gz"),
        vcf_tbi=temp("results/filtered/all.gatk.{vartype}.hardfiltered.vcf.gz.tbi"),
        tranches=temp("results/filtered/calls/all.gatk.{vartype}.tranches"),
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    resources:
        mem_mb=10000,
    threads: 1
    wrapper:
        "v1.23.5/bio/gatk/variantfiltration"

## remove .snakemake/metadata/*
rule recalibrate_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.gatk.{vartype}.vcf.gz",
        vcf_tbi="results/filtered/all.gatk.{vartype}.vcf.gz.tbi",
        hapmap=config["filtering"]["vqsr_resources"]["path"] + "/" + config["filtering"]["vqsr_resources"]["hapmap"],
        dbsnp=config["filtering"]["vqsr_resources"]["path"] + "/" + config["filtering"]["vqsr_resources"]["dbsnp"],
        omni=config["filtering"]["vqsr_resources"]["path"] + "/" + config["filtering"]["vqsr_resources"]["omni"],
        g1k=config["filtering"]["vqsr_resources"]["path"] + "/" + config["filtering"]["vqsr_resources"]["g1k"],
        mills=config["filtering"]["vqsr_resources"]["path"] + "/" + config["filtering"]["vqsr_resources"]["mills"],
    output:
        vcf=temp("results/filtered/all.gatk.{vartype}.recal.vcf"),
        idx=temp("results/filtered/all.gatk.{vartype}.recal.vcf.idx"),
        tranches=temp("results/filtered/all.gatk.{vartype}.tranches"),
    params:
        mode=get_wildcard_vartype,
        resources=get_variant_recalibrator_resources,
        annotation=get_variant_recalibrator_annotation,
        extra=config["params"]["gatk"]["VariantRecalibrator"],
    log:
        "logs/gatk/variantrecalibrator/gatk.{vartype}.log",
    resources:
        mem_mb=10000,
    threads: 1
    wrapper:
        "v1.23.5/bio/gatk/variantrecalibrator"

## remove .sankemake/metadata/*
rule apply_vqsr:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.gatk.{vartype}.vcf.gz",
        vcf_tbi="results/filtered/all.gatk.{vartype}.vcf.gz.tbi",
        recal="results/filtered/all.gatk.{vartype}.recal.vcf",
        idx="results/filtered/all.gatk.{vartype}.recal.vcf.idx",
        tranches="results/filtered/all.gatk.{vartype}.tranches",
    output:
        vcf=temp("results/filtered/all.gatk.{vartype}.recalibrated.vcf.gz"),
        vcf_tbi=temp("results/filtered/all.gatk.{vartype}.recalibrated.vcf.gz.tbi"),
    params:
        mode=get_wildcard_vartype,
    log:
        "logs/gatk/apply_vqsr/gatk.{vartype}.log",
    resources:
        mem_mb=10000,
    threads: 1
    wrapper:
        "v1.23.5/bio/gatk/applyvqsr"



rule merge_calls_gatk:
    input:
        vcfs=expand(
            "results/filtered/all.{algorithm}.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated" if config["filtering"]["vqsr"] else "hardfiltered",
            algorithm=config["params"]["algorithm"],
        ),
        vcfs_tbi=expand(
            "results/filtered/all.{algorithm}.{vartype}.{filtertype}.vcf.gz.tbi",
            vartype=["snvs", "indels"],
            filtertype="recalibrated" if config["filtering"]["vqsr"] else "hardfiltered",
            algorithm=config["params"]["algorithm"],
        ),

    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    resources:
        mem_mb=50000,
    threads: 6
    wrapper:
        "v1.23.5/bio/picard/mergevcfs"
