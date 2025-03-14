rulename = "metabuli"
rule metabuli:
    input:
        contigs = contigs_all,
    output:
        metabuli= OUTDIR /  "{key}/classifiers/metabuli",
    params: 
        database = "."
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/metabuli.yaml"
    shell:
        """
        metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key} --seq-mode 1 --threads {threads}
        """
# NOTE: does it work with gzipped fastafiles? 

# metabuli classify <i:FASTA/Q> <i:DBDIR> <o:OUTDIR> <Job ID> [options]
# - INPUT : FASTA/Q file of reads you want to classify. (gzip supported)
# - DBDIR : The directory of reference DB. 
# - OUTDIR : The directory where the result files will be generated.
# - Job ID: It will be the prefix of result files.  
#
# # Paired-end
# metabuli classify read_1.fna read_2.fna dbdir outdir jobid
#
# # Single-end
# metabuli classify --seq-mode 1 read.fna dbdir outdir jobid
#
# # Long-read 
# metabuli classify --seq-mode 3 read.fna dbdir outdir jobid
#
#   * Important parameters:
#    --threads : The number of threads used (all by default)
#    --max-ram : The maximum RAM usage. (128 GiB by default)
#    --min-score : The minimum score to be classified 
#    --min-sp-score : The minimum score to be classified at or below species rank. 
#    --taxonomy-path: Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
#    --accession-level : Set 1 to use accession level classification (0 by default). 
#                        It is available when the DB is also built with accession level taxonomy.

rulename = "mmseqs2"
rule mmseqs2:
    output:
        mmseqs2 = OUTDIR /  "{key}/classifiers/mmseqs2",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/mmseqs2.yaml"
    shell:
        """
        echo
        """

rulename = "centrifuge"
rule centrifuge:
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
    shell:
        """
        
        -k 1
        """

rulename = "kraken2"
rule kraken2:
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/kraken2.yaml"
    shell:
        """
        echo
        """

