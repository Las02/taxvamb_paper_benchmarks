import pandas as pd
import collections
import os
from pathlib import Path
import sys

# Set the directory the snakefile exists in. This makes us able to call the pipeline with the relevant src files from other directories.
THIS_FILE_DIR = Path(workflow.basedir)

# If the configfile is not set explicit fall back to the default
CONFIG_PATH = THIS_FILE_DIR / "config/config.yaml"
# TODO should go to configfile
if CONFIG_PATH.exists():
    configfile: CONFIG_PATH

# Define the src directory for the files used in the snakemake workflow
SRC_DIR = THIS_FILE_DIR / "files_used_in_snakemake_workflow"  

# Get the output_directory defined by the user or fallback to current directory, which is the default way snakemake handles output directories
OUTDIR = Path("") if config.get("output_directory") is None else Path(config.get("output_directory"))

#### Setting parameters from the config file ####
##  For a more throughout description of what the different config options mean see the /config/config.yaml file

# Default resources used
default_walltime = config.get("default_walltime", "48:00:00")
default_threads = config.get("default_threads", 16)
default_mem_gb = config.get("default_mem_gb", 50)
use_minimap = config.get("use_minimap", True)

# Minimum contig length used
MIN_CONTIG_LEN = int(config.get("min_contig_len", 2000)) 

# Other options
CUDA = True if config.get("cuda") ==  "True" else False
 
## ----------- ##

# Assert that input files are actually passed to snakemake
if config.get("read_file") == None and config.get("read_assembly_dir") == None and config.get("bam_contig") == None:
    print("ERROR: read_file or read_assembly_dir or bam_contig not passed to snakemake as config. Define either. Eg. snakemake <arguments> --config read_file=<read file>. If in doubt refer to the README.md file")
    sys.exit()

# Set default paths for the SPades outputfiles - running the pipeline from allready assembled reads overwrite these values
contigs =  OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.fasta"
contigs_paths =  OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.paths"
assembly_graph = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/assembly_graph_after_simplification.gfa"

# Default paths for bamfiles 
bamfiles_before = OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam"
bamfiles = OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam"
contigs_all = OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz"


# Set default values for dictonaries containg information about the input information
# The way snakemake parses snakefiles means we have to define them even though they will always be present
sample_id = dict()               
sample_id_path= dict() 
sample_id_path_assembly = dict()
sample_id_contig = collections.defaultdict()

read_fw = ""
read_rv = ""
# If the read_file is defined the pipeline will also run SPades and assemble the reads
if config.get("read_file") != None:
    df = pd.read_csv(config["read_file"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    for id, (sample, read1, read2) in enumerate(zip(df["sample"], df.read1, df.read2)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]
    read_fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
    read_rv =  lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]

# If read_assembly dir is defined the pipeline will run user defined SPades output files
if config.get("read_assembly_dir") != None:
    df = pd.read_csv(config["read_assembly_dir"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    sample_id_path_assembly = collections.defaultdict(dict)
    for id, (sample, read1, read2, assembly) in enumerate(zip( df["sample"], df.read1, df.read2, df.assembly_dir)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]
        sample_id_path_assembly[sample][id] = [assembly]

    # Setting the output paths for the user defined SPades files
    contigs =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "contigs.fasta"
    assembly_graph  =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "assembly_graph_after_simplification.gfa"
    contigs_paths  =  lambda wildcards: Path(sample_id_path_assembly[wildcards.key][wildcards.id][0]) / "contigs.paths"

    read_fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
    read_rv =  lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]

if config.get("bam_contig") != None:
    df = pd.read_csv(config["bam_contig"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    contigs = collections.defaultdict(list)
    for id, (sample, bamfile, contig) in enumerate(zip(df["sample"], df.bamfile, df.contig)):
        id = f"sample{str(id)}"
        sample = sample
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [bamfile]
        sample_id_contig[sample] = contig
        contigs[sample].append(contig)
        bamfiles = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
        contigs_all = lambda wildcards: sample_id_contig[wildcards.key]

    # Check that all contigs are the same for each sample
    for sample in contigs.keys():
        for contig in contigs[sample]:
            if contigs[sample][0] != contig:
                print("Not all contigs are the same")
                sys.exit()
            


# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads) 
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime) 
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb) 

rulename = "all"
rule all:
    input:
        # metadecoder = expand(OUTDIR / "{key}/metadecoder/seed.seed", key=sample_id.keys()),
        # metabat = expand(OUTDIR /  "{key}/metabat/depht.txt", key=sample_id.keys()),
        # semibin = expand(OUTDIR /  "{key}/semibin", key=sample_id.keys()),
        comebin = expand(OUTDIR /  "{key}/comebin", key=sample_id.keys()),
        # vamb_default = expand( OUTDIR / "{key}" / 'vamb_default' / 'vae_clusters_unsplit.tsv', key=sample_id.keys()),
        # metabuli=   expand( OUTDIR /  "{key}/classifiers/metabuli",    key=sample_id.keys()),
        # mmseqs2 = expand(OUTDIR /  "{key}/classifiers/mmseqs2", key=sample_id.keys()),
        # kraken2 =  expand( OUTDIR /  "{key}/classifiers/kraken2",          key=sample_id.keys()),     
        # centrifuge = expand(OUTDIR /  "{key}/classifiers/centrifuge",           key=sample_id.keys()),

#### Rules general for all tools ####

# If only reads are passed run metaspades to assemble the reads
rulename = "spades"
rule spades:
    input:
       fw = read_fw, 
       rv = read_rv, 
    output:
       outdir = directory(OUTDIR / "{key}/assembly_mapping_output/spades_{id}"),
       outfile = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.fasta",
       graph = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/assembly_graph_after_simplification.gfa",
       graphinfo  = OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.paths",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda: THIS_FILE_DIR / "envs/spades_env.yaml"
    shell:
       "spades.py --meta "
       "-t {threads} -m 180 "
       "-o {output.outdir} -1 {input.fw} -2 {input.rv} " 
       "-t {threads} --memory {resources.mem_gb} &> {log} " 

# Rename the contigs to keep sample information for later use 
rulename = "rename_contigs"
rule rename_contigs:
    input:
        contigs,
    output:
        OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.renamed.fasta"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    shell:
        """
        sed 's/^>/>S{wildcards.id}C/' {input} > {output} 2> {log}
        """

# Cat the contigs together in one file to later map each pair of reads against all the contigs together
rulename="cat_contigs"
rule cat_contigs:
    input: lambda wildcards: expand(OUTDIR / "{key}/assembly_mapping_output/spades_{id}/contigs.renamed.fasta", key=wildcards.key, id=sample_id[wildcards.key]),
    output: OUTDIR / "{key}/assembly_mapping_output/contigs.flt.fna.gz"
    threads: threads_fn(rulename)
    params: script =  SRC_DIR / "concatenate.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
    conda: THIS_FILE_DIR / "envs/vamb.yaml"
    shell: 
        "python {params.script} {output} {input} --keepnames -m {MIN_CONTIG_LEN} &> {log} "  


if use_minimap:
    INDEX_SIZE = "12G" # get_config("index_size", "12G", r"[1-9]\d*[GM]$")
    # Index resulting contig-file with minimap2
    rule index:
        input:
            # contigs = os.path.join(OUTDIR,"contigs.flt.fna.gz")
            contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
        output:
            mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi")
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
        log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_" + rulename
        conda: 
            THIS_FILE_DIR / "envs/minimap2.yaml"
        shell:
            "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log}"

    # This rule creates a SAM header from a FASTA file.
    # We need it because minimap2 for truly unknowable reasons will write
    # SAM headers INTERSPERSED in the output SAM file, making it unparseable.
    # To work around this mind-boggling bug, we remove all header lines from
    # minimap2's SAM output by grepping, then re-add the header created in this
    # rule.
    rule dict:
        input:
            # contigs = os.path.join(OUTDIR,"contigs.flt.fna.gz")
            contigs = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
        output:
            dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict")  
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
        log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}" + rulename
        conda:
            THIS_FILE_DIR / "envs/samtools.yaml"
        shell:
            "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

    # Generate bam files 
    rule minimap:
        input:
            fw = read_fw,
            rv = read_rv,
            mmi = os.path.join(OUTDIR, "{key}", "contigs.flt.mmi"),
            dict = os.path.join(OUTDIR,"{key}", "contigs.flt.dict"),
            # fq = lambda wildcards: sample2path[wildcards.sample],
            # mmi = os.path.join(OUTDIR,"contigs.flt.mmi"),
            # dict = os.path.join(OUTDIR,"contigs.flt.dict")
        output:
            bam = bamfiles_before, 
            # bam = temp(os.path.join(OUTDIR,"mapped/{sample}.bam"))
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
        log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
        conda:
            THIS_FILE_DIR / "envs/minimap2.yaml"
        shell:
            # See comment over rule "dict" to understand what happens here
            "minimap2 -t {threads} -ax sr {input.mmi} {input.fw} {input.rv} -N 5"
            " | grep -v '^@'"
            " | cat {input.dict} - "
            " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
            " > {output.bam} 2> {log}"

else:
    # Run strobealign to get the abundances  
    rulename = "Strobealign_bam_default"
    rule Strobealign_bam_default:
            input: 
                fw = read_fw,
                rv = read_rv,
                contig = OUTDIR /"{key}/assembly_mapping_output/contigs.flt.fna.gz",
            output:
                bamfiles_before,
            threads: threads_fn(rulename)
            resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
            benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
            log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
            conda: THIS_FILE_DIR / "envs/strobe_env.yaml"
            shell:
                """
                strobealign -t {threads} {input.contig} {input.fw} {input.rv} > {output} 2> {log}
                """

# Sort the bam files and index them
rulename="sort"
rule sort:
    input:
        # OUTDIR / "{key}/assembly_mapping_output/mapped/{id}.bam",
        bamfiles,
    output:
        OUTDIR / "{key}/assembly_mapping_output/mapped_sorted/{id}.bam.sort",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", f"{str(OUTDIR)}/log/") + "{key}_{id}_" + rulename
    conda: THIS_FILE_DIR / "envs/samtools.yaml"
    shell:
        """
	samtools sort --threads {threads} {input} -o {output} 2> {log}
	# samtools index {output} 2>> {log}
	"""


## Include the specific rules for each tool
include: THIS_FILE_DIR / "snakemake_modules/vamb_default.smk"
include: THIS_FILE_DIR / "snakemake_modules/metabat.smk"
include: THIS_FILE_DIR / "snakemake_modules/comebin.smk"
include: THIS_FILE_DIR / "snakemake_modules/metadecoder.smk"
include: THIS_FILE_DIR / "snakemake_modules/semibin.smk"
include: THIS_FILE_DIR / "snakemake_modules/taxonomy_classifiers.smk"
