<div align="center">
   
  <h1>methyl-seek 🔬</h1>
  
  **_Bisulphite-sequencing Methylation Pipeline_**

  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8387343.svg)](https://doi.org/10.5281/zenodo.8387343) [![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/OpenOmics/methyl-seek?color=blue&include_prereleases)](https://github.com/OpenOmics/methyl-seek/releases) [![Docker Pulls](https://img.shields.io/docker/pulls/skchronicles/trimgalore)](https://hub.docker.com/repository/docker/skchronicles/trimgalore)<br>[![docs](https://github.com/OpenOmics/methyl-seek/workflows/docs/badge.svg)](https://github.com/OpenOmics/methyl-seek/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/methyl-seek?color=brightgreen)](https://github.com/OpenOmics/methyl-seek/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/methyl-seek)](https://github.com/OpenOmics/methyl-seek/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, methyl-seek. Its long-term goals: to accurately call CpG sites within whole-genome and cell-free DNA fragments, to perform deconvolution for CpG calls within cell-free DNA fragments, and to boldly identify differentially methylated regions like no pipeline before!
  </i>
</div>

## Overview
Welcome to methyl-seek's documentation! This guide is the main source of documentation for users who are getting started with the [methlyation pipeline](https://openomics.github.io/methyl-seek/).

The **`./methyl-seek`** pipeline is composed of several interrelated pipelines to set up and run different types of analysis. Each of the available pipelines performs different functions:

 * [<code>methyl-seek <b>run</b></code>](https://openomics.github.io/methyl-seek/usage/run/): Run the genome-seek pipeline with your input files.
 * [<code>methyl-seek <b>unlock</b></code>](https://openomics.github.io/methyl-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>methyl-seek <b>cache</b></code>](https://openomics.github.io/methyl-seek/usage/cache/): Cache software containers locally.

**methyl-seek** is a comprehensive bisulphite-sequencing-based DNA methylation pipeline. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As [inputs](usage/run.md), it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/methyl-seek/usage/run/) section of each available sub-command.

For more information about issues or troubleshooting a problem, please check out our [FAQ](https://openomics.github.io/methyl-seek/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/methyl-seek/issues).


## Dependencies
**Requires:** `singularity>=3.5`  `snakemake<8.0`

At the current moment, the pipeline uses a mixture of environment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/u/skchronicles). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/methyl-seek.git
# Change your working directory
cd methyl-seek/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./methyl-seek -h
```

## Contribute

This site is a living document, created for and by members like you. **methyl-seek** is maintained by the members of OpenOmics and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull requests to our [GitHub repository](https://github.com/OpenOmics/methyl-seek).


## Cite

If you use this software, please cite it as below:  

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
@software{neelam_redekar_2023_8387344,
  author       = {Neelam Redekar and Tom Hill and Skyler Kuhn},
  title        = {OpenOmics/methyl-seek: v1.0.0},
  month        = sep,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.8387343},
  url          = {https://doi.org/10.5281/zenodo.8387343}
}
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Neelam Redekar, Tom Hill, & Skyler Kuhn. (2023). OpenOmics/methyl-seek: v1.0.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8387343
```

</details>

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.8387343).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
