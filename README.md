# GDI_Pipeline


This Snakemake pipeline implements the GATK best-practices workflow for calling small germline variants with possibility of using GATK Hard/VQSR filter or using deepvariant variant caller.


Improvements from pipeline (https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling.git)

1) Update software versions;
2) Use GATK SPARK tools, when possible;
3) Set an human genome version more robust;
4) Add GATK VQSR filter version;
5) Add deep variant software call available;

## How to install

```
$ mamba create -c conda-forge -c bioconda --name call_variants snakemake snakedeploy
$ source activate call_variants
$ mkdir -p path/to/project-workdir
$ cd path/to/project-workdir
$ snakedeploy deploy-workflow https://github.com/bioinformatics-ua/GDI_Pipeline.git . --branch main
```

:warning: Important, If you intend to use DeepVariant variant caller check if your server/cluster has glibc >= 2.18 

```
$ ldd --version
```

Snakedeploy will create two folders `workflow` and `config`. The `workflow` as a Snakemake module, and `config` contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main Snakefile in the workflow subfolder.

## How to use

To configure this workflow, modify `config/config.yaml` according to your needs:

* To run GATK - hard filter -> `params: -> algorithm: "gatk"` and `filtering: -> vqsr: false`
* To run GATK - VQSR -> `params: -> algorithm: "gatk"`and `filtering: -> vqsr: true`
* To run DeepVariant -> `params: -> algorithm: "deepvariant"`

:warning: Important, if you want to use "GATK - VQSR" you need to download all the resources present in the file "download_VQSR_resources.sh" to a specific diretory and set that directory in `vqsr_resources: -> path: "<vqsr diretory here>"`.


## Sample and unit sheet

Add samples to `config/samples.tsv`. Only the column sample is mandatory, but any additional columns can be added.
For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. For each unit, define platform, and either one (column fq1) or two (columns fq1, fq2) FASTQ files (these can point to anywhere in your system).
The pipeline will jointly call all samples that are defined, following the GATK best practices.

If you want to try with some examples go to `fastq` diretory and run `download_some_samples_to_test.sh`. The `config/samples.tsv` and `config/units.tsv` config files are already set to these samples.

## How to run

Run only in your server:

```
$ snakemake --jobs 20 --use-conda
```

# Working in a cluster environment

## using SGE

Activate your model for Sun Grid Engine (sGE)

```
$ snakemake --profile sge --jobs 20 --use-conda
```

## using Slurm

Activate your model for Sun Grid Engine (sGE)

```
$ snakemake --profile slurm --jobs 20 --use-conda
```

## Activate a personal profile in a cluster for SnakeMake

To add a cluster environment for SnakeMake you need to install on in your home path:

This is the case for SGE:

```
$ cookiecutter https://github.com/Snakemake-Profiles/sge.git
## add the queue name in the next file
$ vi ~.config/snakemake/sge/cluster.yaml
```

