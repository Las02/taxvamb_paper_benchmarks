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

rule decompress_fastafile:
    input:
        contigs = contigs_all,
    output:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/contigs.flt.fna",
    shell:
        """
        # Decompress fastafile 
        gzip --decompress -c {input.contigs} > {output.contigs_decompressed}
        echo decompression done
        """

rule metadecoder:
    input:
        samfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/contigs.flt.fna",
    output:
        coverage_file = OUTDIR /  "{key}/metadecoder/coverage_file.coverage",
        seed = OUTDIR /  "{key}/metadecoder/seed.seed",
    params:
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
        echo metadecoder coverage done

        # # Get seed
        metadecoder seed --threads {threads} -f {input.contigs_decompressed} -o {output.seed}
        echo metadecoder seed done

        # Actually cluster 
        metadecoder cluster -f {input.contigs_decompressed} -c {output.coverage_file} -s {output.seed} -o {params.metadecoder}
        echo metadecoder cluster done
        """
