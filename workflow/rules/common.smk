import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Helper functions #####

# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group.
     -H '@HD\tVN:1.6\tSO:coordinate'"""
    return r"-H '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}_{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}_{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}_{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{sample}_{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/dedup/{sample}_{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f

def get_wildcard_vartype(wildcards):
    """  only three possibilities for 'gatk --mode' INDEL/SNP/BOTH  """
    return "SNP" if wildcards.vartype == "snvs" else "INDEL"

def get_variant_recalibrator_resources(wildcards):
    """ return resources for gatk variant recalibrator """
    if (wildcards.vartype == "snvs"):
        return {
            "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni": {"known": False, "training": True, "truth": False, "prior": 12.0},
            "g1k": {"known": False, "training": True, "truth": False, "prior": 10.0},
            "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0},
        }
    else:
        return {
            "mills": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0},
        }

def get_variant_recalibrator_annotation(wildcards):
    """ return annotations for gatk variant recalibrator """
    if (wildcards.vartype == "snps"):
        return config["filtering"]["vqsr_annotation_filters_snp"].split(',') + \
           ['InbreedingCoeff'] if len(samples.index) > 10 else []
    else:
        return config["filtering"]["vqsr_annotation_filters_indels"].split(',') + \
           ['InbreedingCoeff'] if len(samples.index) > 10 else []

def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
       get_wildcard_vartype(wildcards)
    )


def get_filter(wildcards):
    """ return hard-filter parameters for SNP or INDEL """ 
    return config["filtering"]["hard"][wildcards.vartype]



