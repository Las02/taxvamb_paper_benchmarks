# /maps/projects/rasmussen/people/bxc755/conda_env/conda/bin/snakemake  -n --snakefile /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/taxvamb_paper_benchmarks/snakefile.smk --rerun-triggers mtime --nolock -c 40 -p --keep-going --software-deployment-method apptainer --use-conda --rerun-incomplete --keep-incomplete --config bam_contig=/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/taxvamb_paper_benchmarks/data_configs/esrum_airways_bamfile_contigs.tsv output_directory=/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways --directory /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways --keep-going
# Config file /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/taxvamb_paper_benchmarks/config/config.yaml is extended by additional config specified via the command line.
# The flag 'temp' used in rule metabat is only valid for outputs, not inputs.
# host: esrumcmpn05fl.unicph.domain
# Building DAG of jobs...
# Your conda installation is not configured to use strict channel priorities. This is however important for having robust and correct environments (for details, see https://conda-forge.org/docs/user/tipsandtricks.html). Please consider to configure strict priorities by executing 'conda config --set channel_priority strict'.
# Job stats:
# job        count
# -------  -------
# all            1
# comebin        1
# total          2
#
#
# [Mon Mar 17 23:53:52 2025]
# rule comebin:
#     input: /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/tmp_contigs_2kbp.fna.gz, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample0.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample1.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample2.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample3.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample4.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample5.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample6.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample7.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample8.bam.sort, /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample9.bam.sort
#     output: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin
#     log: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/log/Airways_comebin
#     jobid: 1
#     benchmark: benchmark/Airways_comebin
#     reason: Missing output files: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin
#     wildcards: key=Airways
#     threads: 16
#     resources: tmpdir=<TBD>, walltime=48:00:00, mem_gb=50
#
#
#
run_comebin.sh -a /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/contigs_2kbp.fna         -o /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin         -p /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample0.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample1.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample2.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample3.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample4.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample5.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample6.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample7.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample8.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample9.bam.sort         -t 16 
# run_comebin.sh -a /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/tmp_contigs_2kbp.fna.gz         -o /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin         -p /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample0.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample1.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample2.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample3.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample4.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample5.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample6.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample7.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample8.bam.sort /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample9.bam.sort         -t 16 

#
#
# [Mon Mar 17 23:53:52 2025]
# rule all:
#     input: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin
#     jobid: 0
#     reason: Input files updated by another job: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/comebin
#     resources: tmpdir=<TBD>
#
# Job stats:
# job        count
# -------  -------
# all            1
# comebin        1
# total          2
#
# Reasons:
#     (check individual jobs above for details)
#     input files updated by another job:
#         all
#     output files have to be generated:
#         comebin
#
# This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
