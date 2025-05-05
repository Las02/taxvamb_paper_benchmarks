# rulename = "semibin"
# rule semibin:
#     input:
#         contigs = contigs_all,
#         bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
#     output:
#         semibin = directory(OUTDIR /  "{key}/semibin"),
#         semibin_bins = OUTDIR /  "{key}/semibin/bins",
#     threads: threads_fn(rulename)
#resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
#     log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
#     conda: THIS_FILE_DIR / "envs/semibin.yaml"
#     shell:
#         """
#         rm -rf {output.semibin}
#         SemiBin2 multi_easy_bin -i {input.contigs} -b {input.bamfiles} -o {output.semibin} \
#         --separator C -t {threads} --write-pre-reclustering-bins --self-supervised # &> {log}
#         """

rulename = "semibin_checkm"
rule semibin_checkm:
    input: 
        bin_dir = OUTDIR /  "{key}/semibin/bins",
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/semibin"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'gz' --database_path {params.database}
        """
