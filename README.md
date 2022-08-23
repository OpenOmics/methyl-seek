# methyl-seek 🔬  [![docs](https://github.com/OpenOmics/methyl-seek/workflows/docs/badge.svg)](https://github.com/OpenOmics/methyl-seek/actions) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/methyl-seek?color=brightgreen)](https://github.com/OpenOmics/methyl-seek/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/methyl-seek)](https://github.com/OpenOmics/methyl-seek/blob/main/LICENSE)

> **_Bisulphite-sequencing methylation pipeline_**. This is the home of the pipeline, methyl-seek. Its long-term goals: to accurately call CpG sites within whole genome and cell-free DNA fragments, to perform deconvolution for CpG calls within cell-free DNA fragments, and to boldly identify differentially methylated regions like no pipeline before!

---
## Overview
Welcome to methyl-seek's documentation! This guide is the main source of documentation for users that are getting started with the [methlyation pipeline](https://openomics.github.io/methyl-seek/).

The **`./methyl-seek`** pipeline is composed several inter-related pipelines to setup and run different types of analysis. Each of the available pipelines perform different functions:

 * [<code>methyl-seek <b>run</b></code>](https://openomics.github.io/methyl-seek/usage/run/): Identify CpG sites from whole genome or cell-free DNA bisulphite sequencing data.
 * [<code>methyl-seek <b>dmr</b></code>](https://openomics.github.io/methyl-seek/usage/run/): Determine differentially methylated regions populated with CpG sites.
 * [<code>methyl-seek <b>dcv</b></code>](https://openomics.github.io/methyl-seek/usage/run/): Identify cells/tissue of origin for cell-free DNA.

**methyl-seek** is a comprehensive bisulphite-sequencing based methylation pipeline. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As [inputs](usage/run.md) , it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/methyl-seek/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/methyl-seek/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/methyl-seek/issues).

## Contribute

This site is a living document, created for and by members like you. **methyl-seek** is maintained by the members of OpenOmics and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/methyl-seek).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). Snakemake-a scalable bioinformatics workflow engine. Bioinformatics 34(20): 3600.</sup>  
