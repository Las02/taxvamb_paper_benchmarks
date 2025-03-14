rulename = "metabuli"
rule metabuli:
    output:
        metabuli= OUTDIR /  "{key}/classifiers/metabuli",
    conda: THIS_FILE_DIR / "envs/metabuli.yaml"
    shell:
        """
        echo
        """

rulename = "centrifuge"
rule centrifuge:
    output:
        centrifuge = OUTDIR /  "{key}/classifiers/centrifuge",
    conda: THIS_FILE_DIR / "envs/centrifuge.yaml"
    shell:
        """
        echo
        """

rulename = "kraken2"
rule kraken2:
    output:
        kraken2 = OUTDIR /  "{key}/classifiers/kraken2",
    conda: THIS_FILE_DIR / "envs/kraken2.yaml"
    shell:
        """
        echo
        """

