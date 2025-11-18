<img src="logo.png" alt="image" width="800" height=auto>

## The official implementation of [ViralQC: A Tool for Assessing Completeness and Contamination of Predicted Viral Contigs]

![GitHub License](https://img.shields.io/github/license/ChengPENG-wolf/ViralQC)

## Contents

- [1. Introduction](#1-introduction)
- [2. Setup Environment](#2-setup-environment)
- [3. Quick Start](#3-quick-start)
- [4. Output explanation](#4-output-explanation)
- [5. Citation](#5-citation)
- [6. Contact](#6-contact)

## 1. Overview

ViralQC is a python library for quality assessment of assembled viral contigs or bins. ViralQC contains two primary modules. The contamination detection module identifies and removes non-viral regions within viral contigs. The completeness module evaluates the expected genome length and calculates the completeness of the viral contigs.

## 2. Setup environment

*Note*: we suggest you install all the packages using [mamba](https://github.com/mamba-org/mamba) or [conda](https://docs.conda.io/en/latest/miniconda.html).

```
# clone the repository to the local
git clone https://github.com/ChengPENG-wolf/ViralQC.git
cd ViralQC

# install and activate environment for ViralQC
conda env create -f viralqc.yaml -n viralqc
conda activate viralqc
# if the above commands fail, try manually installing the environment using the commands in installation.txt

# download and setup the database
python viralqc.py download_database --db /path/to/db
# manually download through Baidu Disk
https://pan.baidu.com/s/19C1EmYJfu3K3jc9WcJQ3Dg?pwd=r9mn
```

## 3. Quick start

**Run Contamination Detection:**

*Note*: The Contamination Detection module should be run on GPUs.

```
python viralqc.py contamination [-i INPUT_FA] [-o OUTPUT_PTH] [-d DATABASE_PATH] [-t THREADS]
```

**Run Completeness Estimation:**

```
python viralqc.py completeness [-i INPUT_FA] [-o OUTPUT_PTH] [-d DATABASE_PATH] [-t THREADS]
```

**Options**

```
  -i INPUT_FA
                        The name of your input file (FASTA format)
  -o OUTPUT_PTH
                        The path of the output directory
  -d DATABASE
                        Model directory
  -t THREADS
                        Number of threads
  -b BIN
                        Run in bin mode (Completeness estimation)
```

## 4. Output explanation

#### 1. `OUTPUT_PTH`/contamination_result.csv:

```
seq_name   seq_len   total_genes   virus_genes   non-virus_genes   virus_length   non-virus_length   region_types       region_coords_bp
--------   -------   -----------   -----------   ---------------   ------------   ----------------   ----------------   ------------------
contig_1   20865     44            32            5                 14196          6669               non-virus;virus;   1-6669;6670-20865;
contig_2   11336     16            10            4                 9881           1455               virus;non-virus;   1-9881;9882-11336;
…
```

This tabular file lists only the inputs with contamination:

- `seq_name`: The identifier of the sequence in the input FASTA file.
- `seq_len`: The length of the sequence.
- `total_genes`: The total number of genes in the sequence.
- `virus_genes`: The number of virus genes in the sequence.
- `non-virus_genes`: The number of non-virus genes in the sequence.
- `virus_length`: The Length of virus regions (bp).
- `virus_length`: The Length of non-virus regions (bp).
- `region_types`: The Type of regions (virus or non-virus) in the sequence (separated by `;`).
- `region_coords_bp`: The start and end coordinates (bp) for each region in the sequence. The coordinates correspond to the order of region_types (separated by `;`).

#### 2. `OUTPUT_PTH`/extracted_virus.fasta:

This FASTA file contains all the virus sequences after removing contamination.

#### 3. `OUTPUT_PTH`/completeness_result.csv:

```
seq_name   seq_len   num_protein   dtr     expected_length   completeness   confidence   significance   ref_id        ref_taxonomy                                              match   mutation   insertion   deletion   translocation   duplication   outlier
--------   -------   -----------   -----   ---------------   ------------   ----------   ------------   -----------   -------------------------------------------------------   -----   --------   ---------   --------   -------------   -----------   -------
contig_1   31816     55            False   111506            28.5           high         19.21          NC_019529.1   Duplodnaviria;...;Demerecviridae;...;Vibrio phage pVp-1   11      8          3           5          0               0             33
contig_2   11324     15            False   12754             88.8           high         17.82          NC_076111.1   Varidnaviria;...;...;Terrapene box turtle adintovirus     6       1          3           1          2               1             6
…
```

This tabular file lists only the inputs with contamination:

- `seq_name`: The identifier of the sequence in the input FASTA file.
- `seq_len`: The length of the sequence.
- `num_protein`: The number of proteins in the sequence.
- `dtr`: Whether the sequence contains direct terminal repeat.
- `expected_length`: The expected length of its complete genome (bp).
- `completeness`: The completeness of the sequence (%).
- `confidence`: Indicate how reliable the prediction result is based on the significance score.
- `significance`: A score indicating the statistical significance of the match between the sequence and the reference genome. A score of 5 corresponds to an e-value of 1e-5.
- `ref_id`: The identifier of the reference genome (e.g., GenBank accession number) that the sequence most closely matches.
- `ref_taxonomy`: The taxonomy lineage of the reference genome.
- `match`: The number of shared proteins.
- `mutation`: The number of proteins that change their homologous family due to the accumulation of nucleotide variations.
- `insertion`: The number of acquisitions of external proteins.
- `deletion`: The number of losses of proteins.
- `translocation`: The number of rearrangements within the genome that change the proteins' positions.
- `duplication`: The number of proteins with more than one copy.
- `outlier`: The number of proteins out of the aligned region.

## 5. Citation

```
@article{peng2025viralqc,
  title={ViralQC: A Tool for Assessing Completeness and Contamination of Predicted Viral Contigs},
  author={Peng, Cheng and Shang, Jiayu and Guan, Jiaojiao and Sun, Yanni},
  journal={arXiv preprint arXiv:2504.05790},
  year={2025}
}
```

## 6. Contact

If you have any questions, please post an issue or email us: cpeng29-c@my.cityu.edu.hk
