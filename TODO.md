# TODO
## Tools to benchmark
### VAMB based tools
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
#### Not done anything with
Comebin : 
Metabuli
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
