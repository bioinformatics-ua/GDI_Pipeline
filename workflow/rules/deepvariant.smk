### If ther more than one bam for sample need to join them, Result: <sample_name>.bam
rule join_bam_files:
    input:
        get_sample_bams,
    output:
        bam="results/joined_bam/{sample}.bam",
##        bam=temp("results/joined_bam/{sample}.bam"),
    log:
        "logs/samtools_merge/{sample}.log",
    params:
        extra="",  # optional additional parameters as string
    benchmark:
        "benchmarks/samtools_merge_bam/{sample}.tsv"
    threads: 8
    wrapper:
        "0.74.0/bio/samtools/merge"

rule deepvariant_gvcf:
    input:
        bam="results/joined_bam/{sample}.bam",
        ref="resources/genome.fasta",
    output:
#        vcf=temp("results/deepvariant/gvcf_calls/{sample}.vcf.gz"),
#        vcf_tbi=temp("results/deepvariant/gvcf_calls/{sample}.vcf.gz.tbi"),
#        gvcf=temp("results/deepvariant/gvcf_calls/{sample}.g.vcf.gz"),
#        gvcf_tbi=temp("results/deepvariant/gvcf_calls/{sample}.g.vcf.gz.tbi"),
        vcf="results/deepvariant/gvcf_calls/{sample}.vcf.gz",
        vcf_tbi="results/deepvariant/gvcf_calls/{sample}.vcf.gz.tbi",
        gvcf="results/deepvariant/gvcf_calls/{sample}.g.vcf.gz",
        gvcf_tbi="results/deepvariant/gvcf_calls/{sample}.g.vcf.gz.tbi",

    params:
        model=config["params"]["library"],   # {wgs, wes, pacbio, hybrid}
        extra=""
    threads: 10
    benchmark:
        "benchmarks/deepvariant/gvcf_calls/{sample}.tsv"
    log:
        "logs/deepvariant/{sample}/stdout.log"
#    shell:
#        To Run in DeepVariant docker 
#        """ BIN_VERSION="1.4.0"; docker run -v `pwd`/resources:/resources -v `pwd`/results/joined_bam:/input -v `pwd`/results/deepvariant:/output google/deepvariant:"${{BIN_VERSION}}" /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/resources/genome.fasta --reads=/input/{wildcards.sample}.bam --output_vcf=/output/{wildcards.sample}.vcf --output_gvcf=/output/{wildcards.sample}.gvcf --num_shards=20 --logging_dir=/output/logs --dry_run=false """

    wrapper:
        "0.74.0/bio/deepvariant"


rule deepvariant_combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand(
            "results/{algorithm}/gvcf_calls/{sample}.g.vcf.gz", sample=samples.index,
            algorithm=config["params"]["algorithm"]
        ),
        gvcfs_tbi=expand(
            "results/{algorithm}/gvcf_calls/{sample}.g.vcf.gz.tbi", sample=samples.index,
            algorithm=config["params"]["algorithm"]
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    benchmark:
        "benchmarks/deepvariant/combine_calls/all.tsv"
    log:
        "logs/deepvariant/gvcf_calls/glnexus_cli.log",
    conda:
        "../envs/glnexus.yaml",
    shell:
        """ glnexus_cli –config DeepVariantWGS {input.gvcfs} | bcftools view - | bgzip -c > {output.vcf} """
##      """ docker run -v `pwd`:/data quay.io/mlin/glnexus:v1.3.1 /usr/local/bin/glnexus_cli --config DeepVariantWGS /data/HG004.g.vcf.gz /data/HG003.g.vcf.gz /data/HG002.g.vcf.gz | bcftools view - | bgzip -c > results/deepvariant/gvcf_calls/all.g.vcf.gz """
##      """ docker run -v `pwd`:/data quay.io/mlin/glnexus:v1.3.1 /usr/local/bin/glnexus_cli --config DeepVariantWGS /data/results/deepvariant/gvcf_calls/NA12877_200.g.vcf.gz /data/results/deepvariant/gvcf_calls/NA12878_200.g.vcf.gz | bcftools view - | bgzip -c > {output.gvcf} ""




