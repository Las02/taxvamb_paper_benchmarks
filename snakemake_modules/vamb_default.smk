# Run vamb 
rulename = "run_vamb"
rule run_vamb:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.bam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'vamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'vamb_default','vae_clusters_unsplit.tsv'),
        compo = os.path.join(OUTDIR, '{key}','vamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    shell:
        """
        rmdir {output.directory}
        vamb bin contr_vamb --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} \
        -m {MIN_CONTIG_LEN} &> {log}
        """
