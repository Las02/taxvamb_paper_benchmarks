all:
	taxvamb_benchmark --help
emtpy_test_data-from_bamfiles:
	taxvamb_benchmark -n --output ./test_results/empty_test_data-from_reads --bam_assembly ./data_configs/bam_assembly.tsv --snakemake_arguments '-p' 
emtpy_test_data-from_reads:
	taxvamb_benchmark -n --output ./test_results/empty_test_data-from_reads --reads ./data_configs/reads.tsv --snakemake_arguments '-p' 
# emtpy_test_data-from_spades:
# 	PlasMAAG -n --output output  --reads_and_assembly_dir reads_and_assembly_dir.tsv
