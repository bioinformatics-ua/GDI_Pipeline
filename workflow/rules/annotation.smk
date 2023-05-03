rule annotate_variants:
    input:
        calls="results/filtered/all.vcf.gz",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
    output:
        calls=report(
            "results/annotated/all.vcf.gz",
            caption="../report/vcf.rst",
            category="Calls",
        ),
        stats=report(
            "results/stats/all.stats.html",
            caption="../report/stats.rst",
            category="Calls",
        ),
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"],
    log:
        "logs/vep/annotate.log",
    threads: 1
    wrapper:
        "v1.23.5/bio/vep/annotate"
