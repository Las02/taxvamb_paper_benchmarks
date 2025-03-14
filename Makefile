all:
	taxvamb_benchmark --help
emtpy_test_data-from_bamfiles:
	taxvamb_benchmark  --output ./test_results/empty_test_data-from_reads --bam_contig ./data_configs/bam_assembly.tsv --snakemake_arguments '-p  --software-deployment-method apptainer'
emtpy_test_data-from_reads:
	taxvamb_benchmark -n --output ./test_results/empty_test_data-from_reads --reads ./data_configs/reads.tsv --snakemake_arguments '-p' 
realsample_data-from_spades:
	taxvamb_benchmark --threads 40 --output ../realsample_data-from_spades  --reads_and_assembly_dir ./data_configs/spades_file.tsv --snakemake_arguments '-p' 
realsample_data-airways_from_bam:
	taxvamb_benchmark --threads 40 --output ../airways  --bam_contig ./data_configs/airways_bamfile_contigs.tsv --snakemake_arguments '-p --software-deployment-method apptainer' 
conda: 
	ls .snakemake/conda/ | grep .yaml | fzf  --preview 'cat .snakemake/conda/{}' | sed 's/.yaml//'

