########### METABULI ####################

rulename = "metabuli"
rule metabuli:
    input:
        contigs = contigs_all,
    output:
        metabuli= directory(OUTDIR /  "{key}/classifiers/metabuli"),
        metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_classifications.tsv",
        metabuli_report= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_report.tsv",
    params: 
        database = config.get("metabuli_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/metabuli.yaml"
    shell:
        """
        metabuli classify {input.contigs} {params.database} {output.metabuli} {wildcards.key}.metabuli --seq-mode 1 --threads {threads}
        """

rulename = "metabuli_taxconv"
rule metabuli_taxconv:
    input:
        metabuli= directory(OUTDIR /  "{key}/classifiers/metabuli"),
        metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_classifications.tsv",
        metabuli_report= OUTDIR /  "{key}/classifiers/metabuli/{key}.metabuli_report.tsv",
    output:
        metabuli_classification= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: "taxconv"
    shell:
        """
        taxconverter metabuli -c {input.metabuli_classification} -r {input.metabuli_report} -o {output.metabuli_classification}
        """


# Run taxvamb 
rulename = "run_taxvamb_metabuli"
rule run_taxvamb_metabuli:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_metabuli"
rule run_taxvamb_no_predictor_metabuli:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy= OUTDIR /  "{key}/classifiers/metabuli/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'metabuli_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','metabuli_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

# Run taxvamb 
rulename = "format_bins_default_vamb"
rule format_bins_default_vamb:
    input:
        contigs = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
        directory = directory(os.path.join(OUTDIR,"{key}", 'vamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'vamb_default','vae_clusters_split.tsv'),
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'vamb_default_bins')),
    params:
        create_fasta = SRC_DIR / "create_fasta.py"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
            rm -rf {output.directory} # clean up dir eg. for failed runs
            python {params.create_fasta} {input.contigs} {input.bins} 0 {output.directory} 
        """

rulename = "metabuli_checkm"
rule metabuli_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_default')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/metabuli_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database"),
        create_fasta = SRC_DIR / "create_fasta.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

rulename = "metabuli_checkm_no_predictor"
rule metabuli_checkm_no_predictor:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'metabuli_taxvamb_no_predictor')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/metabuli_taxvamb_no_predictor"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database"),
        create_fasta = SRC_DIR / "create_fasta.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

########### KRAKEN2 ####################

rulename = "kraken2"
rule kraken2:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",
    params: 
        database = config.get("kraken2_database")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/kraken2.yaml"
    shell:
        """
        kraken2 --minimum-hit-groups 3 --db {params.database} --threads {threads}  {input.contigs_decompressed} > {output.kraken2}
        """

rulename = "kraken_taxconv"
rule kraken_taxconv:
    input:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/kraken2_predictions.tsv",
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: "taxconv"
    shell:
        """
        taxconverter kraken2 -i {input.kraken2} -o {output.kraken2}
        """

# Run taxvamb 
rulename = "run_taxvamb_kraken"
rule run_taxvamb_kraken:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_default')),
        bins = os.path.join(OUTDIR,"{key}",'kraken_taxvamb_default','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','kraken_taxvamb_default/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """


# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_kraken"
rule run_taxvamb_no_predictor_kraken:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/kraken2/taxvamb_formatted_classifications",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'kraken_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','kraken_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """

rulename = "kraken_checkm_no_predictor"
rule kraken_checkm_no_predictor:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_no_predictor')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/kraken_taxvamb_no_predictor"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database"),
        create_fasta = SRC_DIR / "create_fasta.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

rulename = "kraken_checkm"
rule kraken_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'kraken_taxvamb_default')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/kraken_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database"),
        create_fasta = SRC_DIR / "create_fasta.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

########### CENTRIFUGE ####################

rulename = "centrifuge"
rule centrifuge:
    input:
        contigs_decompressed = OUTDIR /  "{key}/metadecoder/{key}_contigs.flt.fna",
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",
    params: 
        database_dir = config.get("centrifuge_database_dir"),
        database_name = config.get("centrifuge_database_name")
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
    shell:
        """
        centrifuge -x {params.database_dir}/{params.database_name} -k 1 -f {input.contigs_decompressed} --threads {threads} > {output.centrifuge}
        """

rulename = "centri_taxconv"
rule centri_taxconv:
    input:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/centrifuge_predictions.tsv",
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: "taxconv"
    shell:
        """
        taxconverter centrifuge -i {input.centrifuge} -o {output.centrifuge}
        """

# Run taxvamb_no_predictor 
rulename = "run_taxvamb_no_predictor_centrifuge"
rule run_taxvamb_no_predictor_centrifuge:
    input:
        contigs = contigs_all,
        bamfiles = lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.sort.bam", key=wildcards.key, id=sample_id[wildcards.key]),
        taxonomy = OUTDIR /  "{key}/classifiers/centrifuge/taxvamb_formatted_classifications.tsv",
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_no_predictor')),
        bins = os.path.join(OUTDIR,"{key}",'centrifuge_taxvamb_no_predictor','vaevae_clusters_split.tsv'),
        compo = os.path.join(OUTDIR, '{key}','centrifuge_taxvamb_no_predictor/composition.npz'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell:
        """
        rm -rf {output.directory}
        vamb bin taxvamb --cuda --no_predictor --taxonomy {input.taxonomy} --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles} # &> {log}
        """



rulename = "centrifuge_checkm_no_predictor"
rule centrifuge_checkm_no_predictor:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_no_predictor')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb_no_predictor"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """

rulename = "centrifuge_checkm"
rule centrifuge_checkm:
    input: 
        bin_dir = directory(os.path.join(OUTDIR,"{key}", 'centrifuge_taxvamb_default')),
    output:
        outdir = directory(OUTDIR /  "{key}/checkm2/centrifuge_taxvamb"),
    threads: threads_fn(rulename)
    params:
        database = config.get("checkm2_database")
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/checkm2.yaml"
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bin_dir} --output-directory {output.outdir} --extension 'fna' --database_path {params.database}
        """
