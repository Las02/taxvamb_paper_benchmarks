
# I want the centrifuge mmseqs results for the 3 other. And the default vamb ones
# In addition to the Airways ones
# Pehabs more?  
# Atleast Airways ~ 82 for taxvamb is expected
# I get simmilar now -> could it be the bug is not fixed ? 

# Check results of my runn
latest and commit: 35788c910 give simmilar


Hi all,

On the CAMI datasets the newest version of TaxVamb has significantly worse performance than an older commit for some of the classifiers (and the results from the preprint)
I also ran the other binners on the Airways dataset to verify whether it was a dataset problem. See below for results:
```
# binbencher (v0.3.0, same as in preprint) (genomes, recall more than 0.9 and precision more than 0.95, strain level)
semibin 95
metabat 38
metadecoder 91
comebin 103

Default vamb: 76 <--- Better (73 before)
Kraken + taxvamb(v.5.0.2): 91  <-- Better (88 before)
Centrifuge + taxvamb(v.5.0.2): 96 <-- Better (84 before)

Metabuli(GTDB v214.1) + taxvamb(v.5.0.2): 87 <---- worse (118 before) 
Metabuli (GTDB v220 + Refseq vira) + taxvamb:(v.5.0.2) 80
Metabuli + taxvamb (4.1.4.dev141+g21bdcdd) 92

MMSEQS (GTDB) + taxvamb(v.5.0.2): 60 <---- worse (86 Before)
MMSEQS results prepared for benchmarks by jakob (earlier version of mmseqs and earlier database)  + taxvamb (v.5.0.2): 67  (86 Before)
MMSEQS results prepared for benchmarks by jakob (earlier version of mmseqs and earlier database)  + d35788c910 => 87 bins (checkm2)
MMSEQS results prepared for benchmarks by jakob (earlier version of mmseqs and earlier database)  + 4.1.4.dev141+g21bdcdd => 91 (checkm2)


earlier commit 4.1.4.dev142+g6c97ece
kraken: 95
centrifuge: 95
metabuli: 88
mmseqs (gtdb): 89


```

I also completed some of the runs for the 3 new biological samples (Black Sea, Forrest Soil and Apple tree).
But based on the preliminary results TaxVamb (v.5.0.2) does not seem to perform better than both 'normal vamb" or semibin. 
But hopefully the performance will improve if we can figure out the bug. Let me know if you can figure out what is wrong, i'll be happy to start running the results whenever, as I have set up a pipeline for it. 

```
                     HQ-bins  MQ-bins
# PRJNA638805 (black sea)
kraken_taxvamb_default    33   224
run_taxvamb_centrifuge    36   311
metadecoder               65   250
metabuli_taxvamb_default  81   346
metabat                   97   288
mmseqs_gtdb_taxvamb       130  502
default_vamb              158  462
semibin                   169  538

# Older commit
run_taxvamb_gtdb/quality_report.tsv     187     511
metabuli_taxvamb_default/quality_report.tsv     102     375

# PRJNA1003562 (apple tree)
metabat                   5   38
kraken_taxvamb_default    7   35
run_taxvamb_centrifuge    7   33
metadecoder               9   42
mmseqs_gtdb_taxvamb       15  80
metabuli_taxvamb_default  17  63
default_vamb              33  79
semibin                   35  95

# PRJNA783873 (forest soil)
metadecoder               23  223
metabat                   25  207
kraken_taxvamb_default    28  168
run_taxvamb_centrifuge    32  237
metabuli_taxvamb_default  33  177
mmseqs_gtdb_taxvamb       37  284
default_vamb              39  210
semibin                   70  478

# Older commit
metabuli_taxvamb_default/quality_report.tsv     32      155

```






(Matching Fig.2 Airways, Genomes in preprint) Program versions: https://docs.google.com/spreadsheets/d/1aDR8dOS99-jp6NL44WY0938LA75t-KBupQotfYhiNEI/edit?usp=sharing
```

Julia code for benchmarking: 
```
level = 0 # at strain level
reference = PATH TO REFERENCE
bins = PATH TO BINS

ref = Reference(reference)
bins = open(i -> Binning(i, ref), bins)
nc_strains = n_recovered(bins, 0.90, 0.95; level=level, assembly=false)
print(nc_strains)


--------







To get a better overview of what could be going on (such that i could compare results with the results from the preprint), I also ran the other binners in addition to the TaxVamb on the CAMI Airways dataset. But again, here I seem to get worse-than-expected results. Here the other binners seems to have a performance nearly identical to the ones described in the preprint (and seen in the TaxVamb fig2 sheets document). On the other hand taxvamb has better performance using Kraken and Centrifuge - but worse using MMSEQS and Metabuili. Note these are results for Vamb and Taxvamb without reclustering.
I tested TaxVamb with an earlier version and here i get simmilar resutls to before, which leads to me to believe that a bug was introduced somewhere.


```
# binbencher (v0.3.0, same as in preprint) (genomes, recall more than 0.9 and precision more than 0.95, strain level)
-----
semibin 95
metabat 38
metadecoder 91
comebin 103

Default vamb: 76 <--- Better (73 before)
Kraken + taxvamb: 91  <-- Better (88 before)
Centrifuge + taxvamb: 96 <-- Better (84 before)

# These much worse
Metabuli(GTDB v214.1) + taxvamb: 87 (118 before) 
Metabuli (GTDB v220 + Refseq vira) + taxvamb: 80

MMSEQS (GTDB) + taxvamb: 60 (86 Before)
MMSEQS results prepared for benchmarks by jakob (earlier version of mmseqs and earlier database)  + taxvamb: 67  (86 Before)
d35788c910 => 87 bins (checkm2)
4.1.4.dev141+g21bdcdd => 91 (checkm2)

(Matching Fig.2 Airways, Genomes in preprint) Program versions: https://docs.google.com/spreadsheets/d/1aDR8dOS99-jp6NL44WY0938LA75t-KBupQotfYhiNEI/edit?usp=sharing
```

Julia code for benchmarking: 
```
level = 0 # at strain level
reference = PATH TO REFERENCE
bins = PATH TO BINS

ref = Reference(reference)
bins = open(i -> Binning(i, ref), bins)
nc_strains = n_recovered(bins, 0.90, 0.95; level=level, assembly=false)
print(nc_strains)
```
---------

It therefore seems likely that either
1) I have made a mistake somewhere 
2) There is a bug somewhere, introduced in a newer version of taxVAMB.
3) We are now using different versions of the classifiers + taxonomy databases which causes us to lose performance, to test this I tested the current version of taxvamb with mmseqs output files Jakob prepared earlier and got 67 bins on Airways (vs 60 on current) - but still far away from the ones reported in the paper. Although they are reported with reclustering - but does it really have such a large effect? 

To assess whether I have made a mistake, here is an minimal example of getting my results using Kraken2 on Airways. 
(The data is also on C2 at /home/projects/ku_00197/data/vambnew)
I ran these results on ESRUM and I have included full paths to the files used. 

```
## Environment
# Vamb and Centrifuge install
mamba create -n test_env -c bioconda -c default python==3.10.2 kraken2==2.14-0
mamba activate test_env
pip install vamb==5.0.2

# Taxconverter install
git clone https://github.com/RasmussenLab/taxconverter
cd taxconverter
pip install -e .
cd data
unzip linage.zip # Notice for me the unzip version on ESRUM did not work and I just unzipped in on my pc where it worked even though the versions "should" be the same.
cd ../..

# Running the tools
kraken2 --minimum-hit-groups 3 --db /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/kraken_database_download --threads 127  /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/metadecoder/Airways_contigs.flt.fna > kraken2_predictions.tsv

taxconverter kraken2 -i kraken2_predictions.tsv -o taxvamb_formatted_classifications

vamb bin taxvamb  --taxonomy taxvamb_formatted_classifications --outdir kraken_taxvamb_default --fasta /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/tmp_contigs_2kbp.fna.gz -p 128 --bamfiles /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/27.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/26.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/23.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/4.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/10.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/11.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/12.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/9.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/7.bam /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Airways/bam/8.bam


(returns 93)
```
My results for my runs can be found in (with checkm2 results in /checkm2 dir):
```
/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/PRJNA783873
/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/PRJNA638805
/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/PRJNA1003562
/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways
```

Script for calculating HQ and MQ 
```
import pandas as pd
from pathlib import Path
def HQ_mq(dataframe):
    return len(dataframe.query("Completeness > 90.0 and Contamination < 5"))
def MQ_mq(dataframe):
    return len(dataframe.query("Completeness >= 50.0 and Contamination < 10"))
def get_MQ_HQ(file: Path, header: bool = True):
    df = pd.read_csv(file, sep="\t", low_memory=False)
    HQ = HQ_mq(df)
    MQ = MQ_mq(df)
    if header:
        print("file\tHQ\tMQ")
    print(f"{file}\t{HQ}\t{MQ}")

```
I'll be leaving for the next 3 weeks, But when i'm back i'll be glad to run the rest of the benchmarks, as I have allready set up a pipeline for doing so :)


```







/maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample0.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample1.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample2.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample3.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample4.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample5.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample6.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample7.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample8.sort.bam /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Airways/assembly_mapping_output/mapped_sorted/sample9.sort.bam 


## HQ bins:
# PRJNA1003562
mmseqs_gtdb_taxvamb/                 15      80
run_taxvamb_centrifuge/                7       33
mmseqs_kalmari_taxvamb/              17      60
default_vamb/                        33      79
metadecoder/                           9       42
mmseqs_trembl_taxvamb/               17      60
metabat/                              5       38
semibin/                              35      95

# PRJNA638805
mmseqs_gtdb_taxvamb/              130     502
trembl_taxvamb/                    136     439
mmseqs_kalmari_taxvamb/           104     425
semibin/                          169     538
metadecoder/                       65      250
metabat/                             97      288

# PRJNA783873
metabat/                    25      207
mmseqs_gtdb_taxvamb/         37      284
semibin/                    70      478
metadecoder/                23      223

# Airways
We have all of t



## MMseqs results on Airways
gtdb_taxvamb/quality_report.tsv 145     251
metabat/quality_report.tsv      22      58
metadecoder/quality_report.tsv  162     233
default_vamb/quality_report.tsv 159     247
semibin/quality_report.tsv      173     296






metabat/vamb_format_bins.tsv
comebin/comebin_res/vamb_format_bins.tsv
semibin/vamb_format_bins.tsv
metadecoder/vamb_format_bins.tsv






I checked the performance of vamb - and it is simmilar for all of the CAMI datasets.



## Oral
Method	Bins
centrifuge_taxvamb_no_predictor	149
metabuli_taxvamb_default	158
metabuli_taxvamb_no_predictor	155
gtdb_taxvamb_default	121
vamb_default	149

# Skin
Method	Cluster Path	Result
centrifuge_taxvamb_no_predictor	vaevae_clusters_split.tsv	126
centrifuge_taxvamb_default	vaevae_clusters_split.tsv	125
metabuli_taxvamb_default	vaevae_clusters_split.tsv	121
metabuli_taxvamb_no_predictor	vaevae_clusters_split.tsv	118
gtdb_taxvamb_default	vaevae_clusters_split.tsv	65

# Urogenital
Method	Cluster Path	Result
centrifuge_taxvamb_no_predictor	vaevae_clusters_split.tsv	114
centrifuge_taxvamb_default	vaevae_clusters_split.tsv	113
metabuli_taxvamb_default	vaevae_clusters_split.tsv	112
metabuli_taxvamb_no_predictor	vaevae_clusters_split.tsv	116
gtdb_taxvamb_default	vaevae_clusters_split.tsv	80



Running:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false centrifuge_taxvamb_no_predictor/vaevae_clusters_split.tsv false
cwd: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Oral
Ran:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false centrifuge_taxvamb_no_predictor/vaevae_clusters_split.tsv false
The version of BinBencherBackend is: 0.3.0
149
Running:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false metabuli_taxvamb_default/vaevae_clusters_split.tsv false
cwd: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Oral
Ran:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false metabuli_taxvamb_default/vaevae_clusters_split.tsv false
The version of BinBencherBackend is: 0.3.0
158
Running:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false metabuli_taxvamb_no_predictor/vaevae_clusters_split.tsv false
cwd: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Oral
Ran:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false metabuli_taxvamb_no_predictor/vaevae_clusters_split.tsv false
The version of BinBencherBackend is: 0.3.0
155
Running:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false gtdb_taxvamb_default/vaevae_clusters_split.tsv false
cwd: /maps/projects/rasmussen/people/bxc755/taxvamb_benchmarks/airways/Oral
Ran:
/home/bxc755/rasmussen/scratch/ptracker/ptracker/bin/julia-1.10.3/bin/julia /maps/projects/rasmussen/scratch/ptracker/Benchmark_vamb_cli/Binbench.jl /maps/projects/ra
smussen/scratch/ptracker/Benchmark_vamb_cli/data/vambnew/Oral/reference.json false gtdb_taxvamb_default/vaevae_clusters_split.tsv false
The version of BinBencherBackend is: 0.3.0
121

Vamb default 149









# PRJNA638805 (black sea)
kraken_taxvamb_default    33   224
run_taxvamb_centrifuge    36   311
metadecoder               65   250
metabuli_taxvamb_default  81   346
metabat                   97   288
mmseqs_gtdb_taxvamb       130  502
default_vamb              158  462
semibin                   169  538

# PRJNA1003562 (apple tree)
metabat                   5   38
kraken_taxvamb_default    7   35
run_taxvamb_centrifuge    7   33
metadecoder               9   42
mmseqs_gtdb_taxvamb       15  80
run_taxvamb_gtdb          15  80
metabuli_taxvamb_default  17  63
trembl_taxvamb            17  60
default_vamb              33  79
semibin                   35  95

# PRJNA783873 (forest soil)
metadecoder               23  223
metabat                   25  207
kraken_taxvamb_default    28  168
run_taxvamb_centrifuge    32  237
metabuli_taxvamb_default  33  177
mmseqs_gtdb_taxvamb       37  284
default_vamb              39  210
semibin                   70  478
