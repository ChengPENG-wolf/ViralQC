<img src="logo.png" alt="image" width="800" height=auto>

## The official implementation of [ViralQC: A Tool for Assessing Completeness and Contamination of Predicted Viral Contigs]

![GitHub License](https://img.shields.io/github/license/ChengPENG-wolf/ViralQC)

## Contents

- [1. Introduction](#1-introduction)
- [2. Setup Environment](#2-setup-environment)
- [3. Quick Start](#3-quick-start)

## 1. Overview

ViralQC is a python library for quality assessment of assembled viral contigs or bins. ViralQC contains two primary modules. The contamination detection module identifies and removes non-viral regions within viral contigs. The completeness module evaluates the expected genome length and calculates the completeness of the viral contigs.

## 2. Setup environment

*Note*: we suggest you install all the packages using [mamba](https://github.com/mamba-org/mamba) or [conda](https://docs.conda.io/en/latest/miniconda.html).

```
# clone the repository to the local
git clone [https://github.com/ChengPENG-wolf/ViraLM.git](https://github.com/ChengPENG-wolf/ViralQC.git)
cd ViralQC

# install and activate environment for ViraLM
conda env create -f viralqc.yaml -n viralqc
conda activate viralqc

# download and setup the database
python viralqc.py download_database --db /path/to/db
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
