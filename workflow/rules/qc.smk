rule fastqc:
    input:
        unpack(get_trimmed_reads),
    output:
        html="results/qc/fastqc_trimmed/{sample}_{unit}.html",
        zip="results/qc/fastqc_trimmed/{sample}_{unit}_fastqc.zip",
    wrapper:
        "v1.23.5/bio/fastqc"


rule samtools_stats:
    input:
        "results/recal/{sample}_{unit}.bam",
    output:
        "results/qc/samtools-stats/{sample}_{unit}.txt",
    log:
        "logs/samtools-stats/{sample}_{unit}.log",
    wrapper:
        "v1.23.5/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            "results/qc/fastqc_trimmed/{u.sample}_{u.unit}_fastqc.zip",
            u=units.itertuples(),
        ),
        expand(
            "results/qc/samtools-stats/{u.sample}_{u.unit}.txt",
            u=units.itertuples(),
        ),
        expand(
            "results/qc/dedup/{u.sample}_{u.unit}.metrics.txt",
            u=units.itertuples(),
        ),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    params:
        use_input_files_only=False,
    log:
        "logs/multiqc.log",
    wrapper:
        "v1.23.5/bio/multiqc"
