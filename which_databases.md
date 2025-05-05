# Kraken:
default kraken database:
```
This will download NCBI taxonomic information, as well as the complete genomes in RefSeq for the bacterial, archaeal, and viral domains, along with the human genome and a collection of known vectors (UniVec_Core)
```
https://github.com/DerrickWood/kraken2/issues/38  | use ftp instead of rsync

# Metabuli:
GTDB chooise for DB:
```
1. GTDB (101 GiB)
- GTDB 214.1 (Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5).
```
 
# Centrifuge:
  I followed this:
```
  * Use centrifuge-download to download genomes from NCBI. The following two commands download the NCBI taxonomy to taxonomy/ in the current directory, and all complete archaeal, bacterial and viral genomes to library/. Low-complexity regions in the genomes are masked after download (parameter -m) using blast+'s dustmasker. centrifuge-download outputs tab-separated sequence ID to taxonomy ID mappings to standard out, which are required by centrifuge-build.
centrifuge-download -o taxonomy taxonomy
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map
To build the index, first concatenate all downloaded sequences into a single file, and then run centrifuge-build:

cat library/*/*.fna > input-sequences.fna

## build centrifuge index with 4 threads
centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv
```

# MMSEQS
3 Default databases:
```
# mmseqs databases
Usage: mmseqs databases <name> <o:sequenceDB> <tmpDir> [options]

  Name                	Type      	Taxonomy	Url
- UniRef100           	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniRef90            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniRef50            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniProtKB           	Aminoacid 	     yes	https://www.uniprot.org/help/uniprotkb
- UniProtKB/TrEMBL    	Aminoacid 	     yes	https://www.uniprot.org/help/uniprotkb     <-----------------------------
- UniProtKB/Swiss-Prot	Aminoacid 	     yes	https://uniprot.org
- NR                  	Aminoacid 	     yes	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                  	Nucleotide	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- GTDB                  Aminoacid	     yes	https://gtdb.ecogenomic.org  <-----------------------------
- PDB                 	Aminoacid 	       -	https://www.rcsb.org
- PDB70               	Profile   	       -	https://github.com/soedinglab/hh-suite
- Pfam-A.full         	Profile   	       -	https://pfam.xfam.org
- Pfam-A.seed         	Profile   	       -	https://pfam.xfam.org
- Pfam-B              	Profile   	       -	https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released
- CDD                   Profile                -        https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
- eggNOG              	Profile   	       -	http://eggnog5.embl.de
- VOGDB                 Profile                -        https://vogdb.org
- dbCAN2              	Profile   	       -	http://bcb.unl.edu/dbCAN2
- SILVA                 Nucleotide           yes        https://www.arb-silva.de
- Resfinder           	Nucleotide	       -	https://cge.cbs.dtu.dk/services/ResFinder
- Kalamari            	Nucleotide	     yes	https://github.com/lskatz/Kalamari       <-----------------------------
```
