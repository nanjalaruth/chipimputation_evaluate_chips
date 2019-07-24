# Chip imputation evaluation Workflow h3abionet/chipimputation_evaluate_chips

[![Build Status](https://travis-ci.org/h3abionet/chipimputation.svg?branch=master)](https://travis-ci.org/h3abionet/chipimputation)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/chipimputation.svg)](https://hub.docker.com/r/h3abionet/chipimputation)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is to evaluate the imputation performance and accuracy of different arrays starting from sequence data. 
It masks non tag variants for each array, and then impute to a reference panel using Minimac.  
It is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.   
It comes with singularity containers making installation trivial and results highly reproducible.


### Documentation
The evaluate_chips pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and Configuration](docs/installation.md)
2. [Configuration for other clusters](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)

### Setup (native cluster)

#### Headnode
  - [Nextflow](https://www.nextflow.io/) (can be installed as local user)
   - NXF_HOME needs to be set, and must be in the PATH
   - Note that we've experienced problems running Nextflow when NXF_HOME is on an NFS mount.
   - The Nextflow script also needs to be invoked in a non-NFS folder
  - Java 1.8+

#### Compute nodes

- The compute nodes need to have singularity installed.
- The compute nodes need access to shared storage for input, references, output
- The following commands need to be available in PATH on the compute nodes, in case of unavailabitity of singularity.

  - `minimac4` from [MINIMAC4](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
  - `vcftools` from [VCFtools](https://vcftools.github.io/index.html)
  - `bcftools`from [bcftools](https://samtools.github.io/bcftools/bcftools.html)
  - `bgzip` from [htslib](http://www.htslib.org)
  - `eagle` from [Eagle](https://data.broadinstitute.org/alkesgroup/Eagle/)
  - `python2.7`
  - `R` with the following packages ...


