# This older version of metadecoder used in the paper only takes in sam files..
rule bam_to_sam:
    input:
        bamfile = OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.bam.sort"
    output:
        samfile = OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}" + rulename
    conda: THIS_FILE_DIR / "envs/samtools.yaml"
    shell: 
        """
        samtools view -h {input.bamfile} > {output.samfile} 2> {log}
        """

rule metadecoder:
    input:
        contigs = contigs_all,
        samfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        coverage_file = OUTDIR /  "{key}/metadecoder/coverage_file.coverage",
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/contigs.flt.fna",
        seed = OUTDIR /  "{key}/metadecoder/seed.seed",
        metadecoder = OUTDIR /  "{key}/metadecoder/clusters.metadecoder",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/meta_decoder.yaml"
    shell:
        """
        # Calculate coverage
        metadecoder coverage --threads {threads} -s {input.samfiles} -o {output.coverage_file}

        # Decompress fastafile 
        gzip --decompress {input.contigs} > {output.contigs_decompressed}

        # Get seed
        metadecoder seed --threads {threads} -f {output.contigs_decompressed} -o {output.seed}

        # Actually cluster 
        metadecoder cluster --threads {threads} -f {output.contigs_decompressed} -c {output.coverage_file} -s {output.seed} -o {output.metadecoder}
        """
# NOTE: Likely needs sam files as input >(
