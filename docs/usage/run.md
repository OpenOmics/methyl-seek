# <code>methyl-seek</code>

## 1. Overview
The `methyl-seek` executable is composed of several inter-related pipelines intended to analyze human whole genome and cell-free DNA bisulphite-sequencing data to accurately locate CpG methylation sites, identify differentially methylated regions, and perform CpG deconvolution to determine cell/tissue of origin for cell-free DNA fragments. The expected input is paired-end Illumina FASTQ files, and output varies with the execution mode chosen.

#### Execution modes

* [<code><b>run</b></code>](./run.md): Performs quality control and CpG methylation calling from whole genome or cell-free DNA bisulphite sequencing data.
* [<code><b>dmr</b></code>](./dmr.md): Identifies differentially methylated regions (DMRs) populated with CpG sites between experimental groups.
* [<code><b>dcv</b></code>](./deconvolution.md): Identifies cells/tissue of origin for the cell-free DNA fragments based on their CpG methylation profiles.

#### Programs

- FastQC is used to assess the sequencing quality. FastQC is ran twice, before and after adapter trimming. It generates a set of basic statistics to identify problems that can arise during sequencing or library preparation. FastQC will summarize per base and per read QC metrics such as quality scores and GC content. It will also summarize the distribution of sequence lengths and will report the presence of adapter sequences.

- TrimGalor
{is used to remove adapter sequences, perform quality trimming, and remove very short sequences that would otherwise multi-map all over the genome prior to alignment.}

- Kraken2 and FastQ Screen are used to screen for various sources of contamination. During the process of sample collection to library preparation, there is a risk for introducing wanted sources of DNA. FastQ Screen compares your sequencing data to a set of different reference genomes to determine if there is contamination. It allows a user to see if the composition of your library matches what you expect. Also, if there are high levels of microbial contamination, Kraken can provide an estimation of the taxonomic composition. Kraken can be used in conjunction with Krona to produce interactive reports.

- BBmerge **(need text here!)**

- Bismark is used to map bisulfite treated sequencing reads to bisulphite marked reference genome of interest (default: human hg38 genome, hg38) and perform cytosine methylation calls in CpG, CHG and CHH context.

- MultiQC is used to aggregate the results of all above mentioned tools into a single interactive report.

- bsseq (R package) is used for  **(need text here!)**

- **(add other tools here .....)**

#### Flowchart
**(_Add workflow flowchart here !!!)**


## 2. Download repository

Download methyl-seek repository from github: `https://github.com/OpenOmics/RNA-seek/archive/refs/heads/main.zip`

```
mkdir ~/project
cd ~/project

wget https://github.com/OpenOmics/RNA-seek/archive/refs/heads/main.zip
unzip main.zip
cd methyl-seek-main
```

## 3. Initial setup files

#### FASTQ inputs

Pipeline accepts paired-end Illumina FASTQ files named as below: `{sample}.R1.fastq.gz` and `{sample}.R2.fastq.gz`

#### Sample groups

Sample group information is provided in `samples.txt`, which is a tab-delimited, two-column formatted with header labels as `samples` and `group`. First column includes `{sample}` names matching to those in FASTQ file names. Second column includes group labels corresponding to each sample (e.g.: `group1`, `group2` or `test`, `control`).

| samples | group  |
| ------- | ------ |
| S1      | group1 |
| S2      | group1 |
| S3      | group2 |
| S4      | group2 |
| S5      | group3 |
| S6      | group3 |

For DMR analyses, user should provide `group1` and `group2` labels, that should match with sample group information. file.

#### Config file

Edit the `config.yaml` file as below:

- Set `rawdata_dir` to the absolute path of the directory containing all your fastqs.
- Set `result_dir` to the absolute path of the working directory containing `methyl-seek` pipeline, where all result files will be stored.
- Set `samples` to be the absolute path of your `samples.txt` located within `result_dir`.

For example,

```
samples: "~/project/methyl-seek-main/samples.txt"
rawdata_dir: "~/project/fastq/"
result_dir: "~/project/methyl-seek-main/"
```
See below instructions for adding reference to `config.yaml` file.

## 4. Custom reference

By default, bisulphite marked human hg38 reference genome is provide with this repo. To create a custom bisulphite marked reference for different genome of interest, simply provide genome fasta file (for e.g.`hg19.genome.fa`), and update the information in `config.yaml` as below:

```
bisulphite_ref: "~/project/bisulphite_genome"
bisulphite_fa: "~/project/bisulphite_genome/hg19.genome.fa"
species: "hg19"
```

## 5. Execution

Still in progress.. how to use execution modes?
```bash
# Step 1.) Dry-run the pipeline
sinteractive -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

cd ~/project/methyl-seek-main/

## run : generate CpG reports
sbatch pipeline_launch.sh run npr ./

## dcv : perform CpG deconvolution 
sbatch pipeline_launch.sh dcv npr ./

## dmr : perform CpG deconvolution
sbatch pipeline_launch.sh dmr npr ./

# Step 2.) To launch pipeline
cd ~/project/methyl-seek-main/

## run : generate CpG reports
sbatch pipeline_submit.sh run process ./

## dcv : perform CpG deconvolution 
sbatch pipeline_submit.sh dcv process ./

## dmr : perform CpG deconvolution 
sbatch pipeline_submit.sh dmr process ./ group1 group2

```
