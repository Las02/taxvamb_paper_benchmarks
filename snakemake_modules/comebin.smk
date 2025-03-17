# https://github.com/BigDataBiology/SemiBin
rulename = "comebin"
rule comebin:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.bam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/contigs.flt.fna",
    output: 
        comebin = OUTDIR /  "{key}/comebin",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/comebin.yaml"
    shell:
        """
        run_comebin.sh -a {input.contigs_decompressed} \
        -o {output.comebin} \
        -p {input.bamfiles} \
        -t {threads} &> {log}
        """
