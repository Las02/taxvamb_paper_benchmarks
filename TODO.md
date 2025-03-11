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
- Centrifuge
For Centrifuge, we used the centrifuge command with -k 1 flag. 
https://github.com/infphilo/centrifuge/
Centrifuge, Kraken2 and MetaMaps were configured to use NCBI identifiers. 
- Kraken2: https://github.com/DerrickWood/kraken2
For Kraken2, we used the kraken command with –minimum-hit-groups 3 flag. 
Centrifuge, Kraken2 and MetaMaps were configured to use NCBI identifiers. 
- MMSeqs2
For MMseqs2, we used the mmseqs taxonomy command. 
MMseqs2 and Metabuli were configured to use GTDB v207 as the reference database. 
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
#### Tested
#### Not tested: as snakemake
Metabat : Works with docker img : # NOTE: we are using earlier version than in paper
Metadecoder: V: 1.0.19 as in paper: 
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

