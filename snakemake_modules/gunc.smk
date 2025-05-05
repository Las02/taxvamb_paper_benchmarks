# Decompress semibin bins, as they are gzipped, and GUNC needs fasta format# Decompress semibin bins, as they are gzipped, and GUNC needs fasta format
rulename = "gzip_semibins"
rule gzip_semibins:
    input:
        semibin_bins = OUTDIR /  "{key}/semibin/bins",
    output: 
        touch_me = OUTDIR / "{key}/tmp/semibin_decompressed.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        gzip -dk {input.semibin_bins}/*.gz
        touch {output.touch_me}
        """

# name: [filepath, bin suffix]
all_bin_dirs = {
            # Binners
            "semibin": ["{key}/semibin/bins", "fa"],
            "vamb": [os.path.join(OUTDIR,"{key}", 'vamb_default_bins_filtered'), ".fna"],
            "metabat": ["{key}/metabat/metabat", "fa"],
            "metadecoder": ["{key}/metadecoder/clusters", "fasta"],
            # Main classifiers
            "taxvamb_gtdb_mmseqs": [OUTDIR / "{key}/gtdb_taxvamb_default_bins_filtered", "fna"],
            # metabuli
            # centrifuge
            # Main classifiers with --no predictor
            # metabuli, centrifuge, gtdb

            # Extra test classifiers
            "taxvamb_kalmari_mmseqs": [OUTDIR / "{key}/kalmari_taxvamb_default_bins_filtered", "fna"],
            "taxvamb_trembl_mmseqs": [OUTDIR / "{key}/trembl_taxvamb_default_bins_filtered", "fna"],
                }
bin_dir_names = all_bin_dirs.keys()

rulename = "gunc"
rule gunc:
    input:
        # semi_bin_decompressed = OUTDIR / "{key}/tmp/semibin_decompressed.done", # recomment me
        bin_dir = lambda wildcards: all_bin_dirs[wildcards.bin_dir][0]
    output: 
        dir = directory(OUTDIR / "{key}/gunc/{bin_dir}"),
        file = OUTDIR / "{key}/gunc/{bin_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    params:
        suffix = lambda wildcards: all_bin_dirs[wildcards.bin_dir][1],
        database = config.get("gunc_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        # Input dir should contain files in fasta format -- need to gzip for some
        rm -rf {output.dir}
        mkdir -p {output.dir}
        gunc run --input_dir {input.bin_dir} --file_suffix {params.suffix} --db_file {params.database} --threads {threads} --out_dir {output.dir} 
        """

rulename = "collect_gunc"
rule collect_gunc:
    input:
        expand(OUTDIR / "{{key}}/gunc/{bin_dir}",  bin_dir = bin_dir_names)
    output: 
        OUTDIR / "{key}/tmp/gunc.done"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename), gpu=gpu_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/gunc.yaml"
    shell:
        """
        touch  {output}
        """
