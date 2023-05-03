rule get_genome:
    output:
        "resources/genome_38.fasta.gz",
    log:
        "logs/get-genome.log",
    cache: True
    params:
        url=config['ref']['url'],
    shell:
        ## Link to selected reference genome (from https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use):
        "wget {params.url} -O {output}"


rule unzip_genome:
    input:
        "resources/genome_38.fasta.gz",
    output:
        "resources/genome.fasta",
    cache: True
    log:
        "logs/unzip-genome.log",
    shell:
        "gzip -cd {input} > {output}"


checkpoint genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.23.5/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        ## need to chr to the names replace "MT" by "M"
        """zcat {input} | sed 's/##contig=<ID=chr/##contig=<ID=/' | sed 's/^MT/M/g' | rbt vcf-fix-iupac-alleles | bcftools view -Ov | sed 's/##contig=<ID=/##contig=<ID=chr/' | awk '{{ if ( $0 ~ /#/ ) {{ print $0 }} else {{ printf "chr%s\\n", $0 }} }}' | bgzip > {output}"""


rule tabix_known_variants:
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "v1.23.5/bio/tabix/index"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v1.23.5/bio/bwa/index"

rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/vep/plugins"


