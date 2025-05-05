
Kraken downloaded from
https://benlangmead.github.io/aws-indexes/k2:
Standard	Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core	12/28/2024

# Kraken download is at : /people/bxc755/taxvamb_benchmarks/airways/slurm-779675.out

# Now i'm trying to see the results of using centrifuge
- Depending on how they look i'll see if I need to fix mmseqs
- We are running checkm2 in the faster way. - it'll also be easier to produce results 




## IMPORTANT
# TODO mmseq GTDB needs to be rerun for all !!
Plan: 
- Set over for one dataset at a time. 
- If it works. 
- Make sure to use GPU
- Esrum has 8 GPUS imma use 4 of them
- 4 jobs max, with avg. 128 cores

When done the only thing we should need is gunc on GTDB mmseqs







# Apr 9
sludge + human  bare brug hvad der ligger 
vamb + taxvamb rerun (mmseqs2) 

Rerun metabat with v.2.12.1
Is kraken version correct? 



# TODO gunc
Run on bio data... Make sure not to rerun anything !!



# New
- rerun create bins for taxvamb ::: and checkm2 for them, as they use vae_unsplit.tsv bins
- rerun vamb with 5.0.2

# Apr 5
Update input data, then run vamb + taxvamb (gtdb) on ALL datasets. 
- Make sure vamb and taxvamb are being rerun.

<!-- Run Centrifuge, Metabuli and Kraken2 on 1 dataset to get results which can be used to test input to taxvamb  -->
<!-- - Likely needs high memory usage -->
<!-- - Figure out which versions of the databases they are using (/ that i downloaded) -->
<!-- - # Test if download of Kraken2 always have problem..  -->

run GUNC On bins. - what does it take as input ? 
- Which datasets do we want to run GUNC on ?

## When semibin is done 
- Update it's input files, such that we can run checkm2 on it

## WHICH VERSIONS OF EACH DATABASE IS DOWNLOADED ? ? 

## Metabuli
1. GTDB (101 GiB)
- GTDB 214.1 (Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5).
(when running says:  DB name: GTDB214.1+humanT2T)
## Kraken NCBI (downloaded march 30) 
RefSeq Release 229, is the one as it's released 10 match

## Centrifuge NCBI  (end march downloaded)
RefSeq Release 229

## MMSeqs
# Trembl v 2025_01
```
`UniProt Release 2025_01

The UniProt consortium European Bioinformatics Institute (EBI), SIB Swiss 
Institute of Bioinformatics and Protein Information Resource (PIR), 
is pleased to announce UniProt Knowledgebase (UniProtKB) Release 
2025_01 (05-Feb-2025). UniProt (Universal Protein Resource) is a 
comprehensive catalog of information on proteins.

UniProtKB Release 2025_01 consists of 253,206,171 entries (UniProtKB/Swiss-Prot: 
572,970 entries and UniProtKB/TrEMBL: 252,633,201 entries)
UniRef100 Release 2025_01 consists of 453,950,711 entries
UniRef90 Release 2025_01 consists of 204,806,910 entries
UniRef50 Release 2025_01 consists of 69,290,910 entries
UniParc Release 2025_01 consists of 916,871,247 entries, where 852,312,545 are active and 64,558,702 inactive
UniProt databases can be accessed from the web at http://www.uniprot.org and 
downloaded from http://www.uniprot.org/downloads. Detailed release 
statistics for TrEMBL and Swiss-Prot sections of the UniProt Knowledgebase 
can be viewed at http://www.ebi.ac.uk/uniprot/TrEMBLstats/ and 
http://web.expasy.org/docs/relnotes/relstat.html respectively.
```
: UniProt Knowledgebase Release 2025_01 consists of:
UniProtKB/Swiss-Prot Release 2025_01 of 05-Feb-2025
UniProtKB/TrEMBL Release 2025_01 of 05-Feb-2025
 if notExists "${TMP_PATH}/uniprot_trembl.fasta.gz"; then
 │ downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${TMP_PATH}/ver
 │ downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz" "${T
 fi
# GTDB v220
Released Apr 24, 2024
 │   if notExists "${TMP_PATH}/download.done"; then
 │   │   downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION.txt" "${TMP_PATH}/
 │   │   downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_pr
# Kalmari (v3.7)
 │ "Kalamari")
 │ │   if notExists "${TMP_PATH}/kalamari.tsv"; then
 │ │   │   printf "3.7 %s\n" "$(date "+%s")" > "${TMP_PATH}/version"
 │ │   │   downloadFile "https://raw.githubusercontent.com/lskatz/Kalamari/18d71da740546ba4a5117682e1ae2a037379afe0/src/Kalamari_v3.7.tsv" "${TMP_PATH
 │ │   fi
 │ │   ACCESSIONS=""

# Questions
Minimap2: By using -F 3584, you are excluding reads that have any of these flags set. In other words, you are keeping only reads that are not supplementary, not duplicates, and not failing QC.
Read QC ? 

# TODO
## Minimap2 for alignment


## Tools to benchmark
# Taxvamb classifiers
- Taxometer -> taxvamb
All the taxonomic annotations were first refined with Taxometer 42 (v.b5fd0ea) with the default parameters (epochs 100, batch size 1024)

- Metabuli
For Metabuli, we used the metabuli classify command with –seq-mode 1
MMseqs2 and Metabuli were configured to use GTDB v207 as the reference database. 
https://github.com/steineggerlab/Metabuli
GTDB: 214.1 is availible

- Centrifuge
For Centrifuge, we used the centrifuge command with -k 1 flag. 
https://github.com/infphilo/centrifuge/
Centrifuge, Kraken2 and MetaMaps were configured to use NCBI identifiers. 
- Kraken2: https://github.com/DerrickWood/kraken2
For Kraken2, we used the kraken command with –minimum-hit-groups 3 flag. 
Centrifuge, Kraken2 and MetaMaps were configured to use NCBI identifiers. 

# Not conda
- MMSeqs2 : version not on conda
For MMseqs2, we used the mmseqs taxonomy command.  -> also uses prebuilt
MMseqs2 and Metabuli were configured to use GTDB v207 as the reference database. 
docker pull ghcr.io/soedinglab/mmseqs2:14-7e284 # correct version

```
We obtained the taxonomic annotations for contigs of all seven short-read and
two long-read datasets from MMseqs2 (v.7e2840) 33 , Metabuli (v.1.0.1) 79 , Centrifuge (v.1.0.4) 80 and Kraken2 (v.2.1.3) 81 . 
Centrifuge, Kraken2 and MetaMaps were configured to use NCBI identifiers. 
```

### VAMB based tools : Taxvamb + vamb should be latest version.
VAMB
VAMB (reclustering)
Metabuli (species) + TaxVAMB
Metabuli+TaxVAMB
Metabuli+TaxVAMB (reclustered)
Centrifuge+TaxVAMB
Centrifuge+TaxVAMB (reclustered)
Kraken2+TaxVAMB
Kraken2+TaxVAMB (reclustered)
MMSeqs+TaxVAMB
MMSeqs+TaxVAMB (reclustered)

### OTHER
#### Tested eg. it runs
Metadecoder: V: 1.0.19 as in paper: 
#### Not tested: as snakemake
Metabat : Works with docker img : # NOTE: we are using earlier version than in paper
Comebin :  Comebin (v.1.0.3)
#### Not done anything with
<!-- Metabuli --> # Classifier
SemiBin (NO reclustering)
SemiBin (WITH reclustering)

### Should we benchmark ? :: are included in list
AVAMB
AVAMB (reclustering)
Stacked VAE (Centrifuge)
Stacked VAE (Kraken)
Stacked VAE (MMseqs)
Stacked VAE (Metabuli)

# IMPORTANT consideration
- right now we filter contig size.. what should be done? 

## Other reviewer comments: KRAKEN with different databases
Figure out how to change databases with kraken :: 

# Metadecode
  # You may need to install pip3 before. #
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
  python3 get-pip.py

  # You may need to install or upgrade setuptools and wheel using pip3 before. #
  pip3 install --upgrade setuptools wheel

  # Download and install MetaDecoder version 1.2.0 #
  pip3 install -U https://github.com/liu-congcong/MetaDecoder/releases/download/v1.2.0/metadecoder-1.2.0-py3-none-any.whl # new
  https://github.com/liu-congcong/MetaDecoder/releases/download/v1.0.19/metadecoder-1.0.19-py3-none-any.whl # our version

# Testing TaxVAMB:
  @Lasse Schnell Danielsen
   maybe a good place to start is then to re-run TaxVamb from v 5.0.1 versus TaxVamb from commit 5f2cd7 on CAM

# fejl i paper
mentions MetaMaps: 482
-segmode 1 written, should be --segmode 1

# Ting i taxconverter
- does it work for mmseqs2 -- i dont think so

