rulename = "semibin"
rule semibin:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.bam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        semibin = OUTDIR /  "{key}/semibin",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/semibin.yaml"
    shell:
        """
        SemiBin2 multi_easy_bin -i {input.contigs} -b {input.bamfiles} -o {output.semibin} \
        --separator C -t {threads} --write-pre-reclustering-bins --self-supervised # &> {log}
        """
