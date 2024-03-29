# <code>methyl-seek</code>

## 1. Overview
The `methyl-seek` executable is composed of several inter-related pipelines intended to analyze human whole genome and cell-free DNA bisulphite-sequencing data to accurately locate CpG methylation sites, identify differentially methylated regions, and perform CpG deconvolution to determine cell/tissue of origin for cell-free DNA fragments. The expected input is paired-end Illumina FASTQ files, and output varies with the execution mode chosen.

#### Execution modes

* [<code><b>run</b></code>](./run.md): Performs quality control and CpG methylation calling from whole genome or cell-free DNA bisulphite sequencing data.
* [<code><b>dmr</b></code>](./dmr.md): Identifies differentially methylated regions (DMRs) populated with CpG sites between experimental groups.
* [<code><b>dcv</b></code>](./deconvolution.md): Identifies cells/tissue of origin for the cell-free DNA fragments based on their CpG methylation profiles.

#### Programs

- `FastQC` is used to assess the sequencing quality. `FastQC` is ran twice, before and after adapter trimming. It generates a set of basic statistics to identify problems that can arise during sequencing or library preparation. FastQC will summarize per base and per read QC metrics such as quality scores and GC content. It will also summarize the distribution of sequence lengths and will report the presence of adapter sequences.

- `TrimGalor` is used to clip 10 bp ends from both 5' and 3' ends of each sequences, and eliminate any reads shorter than 50 bp length, and runs `FASTQC` on trimmed reads. `BBmerge` is used to get insert sizes of paired read sequences.

- `Kraken2` and `FQScreen` are used to screen for various sources of contamination. During the process of sample collection to library preparation, there is a risk for introducing wanted sources of DNA. FastQ Screen compares your sequencing data to a set of different reference genomes to determine if there is contamination. It allows a user to see if the composition of your library matches what you expect. Also, if there are high levels of microbial contamination, Kraken can provide an estimation of the taxonomic composition. Kraken can be used in conjunction with Krona to produce interactive reports.

- `Bismark` is used to map bisulfite treated sequencing reads to bisulphite marked reference genome of interest (default: human hg38 genome, hg38) and perform cytosine methylation calls in CpG, CHG and CHH context.

- `MultiQC` is used to aggregate the results of all above mentioned tools into a single interactive report.

- (dmr mode) `bsseq` (R package) is used to perform differential methylation calls using CpGs that exist in 25% of the population with coverage of at least 2. Differentially methylated CpGs are used to determine differentially methylated windows that are at least 300 bp long and contain 3 significant CpGs (with significance value of (p<0.05). Comb-P is used to merge adjacent differentially methylated windows (at most 6-bp apart) into differentially methylated regions (DMR).

- (dcv mode) `meth_atlas` (array-based deconvolution) algorithm is used to compare the CpG methylation profiles with known methylation profiles of 25 tissues types (determined using methylation array) to deconvole the source of DNA.

- (dcv mode) `UMX` (WGBS-based deconvolution) algorithm is used to compare the CpG methylation profiles with known methylation profiles of 32 tissues types (determined using whole genome bisulfite sequencing data) to deconvole the source of DNA.

#### Flowchart
**(_Add workflow flowchart here !!!)**


## 2. Download repository

Download methyl-seek repository from github: `https://github.com/OpenOmics/RNA-seek/archive/refs/heads/main.zip`

```
mkdir ~/project
cd ~/project

wget https://github.com/OpenOmics/methyl-seek/archive/refs/heads/main.zip
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

For DMR analyses, user should provide `group1` and `group2` labels, that should match with sample group information. This information will be used to create the `contrasts.txt` file, which is a tab-delimited, four-column formatted with header labels as `sample`, `group`,`comparisons`, and `path`. First column includes `{sample}` names matching to those in FASTQ file names. Second column includes group labels corresponding to each sample (e.g.: `group1`, `group2` or `test`, `control`). Third column includes group labels that are to be compared for DMR analyses. Fourth column includes path to CpG output files generated from bismark (run mode) in `~/project/methyl-seek-main/CpG/{Sample}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz`

Make sure should provide `group` and `comparison` labels match with labels in `samples.txt` file.
Make sure `contrasts.txt` file is stored in main result directory(for e.g.: `~/project/methyl-seek-main/`)

```
cat contrasts.txt
```

| samples | group  |   comparison   |                                         path                                     |
| ------- | ------ | -------------- | ---------------------------------------------------------------------------------|
| S1      | group1 | group1vsgroup2 | ~/project/methyl-seek-main/CpG/S1.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |
| S2      | group1 | group1vsgroup2 | ~/project/methyl-seek-main/CpG/S2.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |
| S3      | group2 | group1vsgroup2 | ~/project/methyl-seek-main/CpG/S2.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |
| S4      | group2 | group1vsgroup2 | ~/project/methyl-seek-main/CpG/S4.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |
| S5      | group3 | group1vsgroup3 | ~/project/methyl-seek-main/CpG/S5.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |
| S6      | group3 | group1vsgroup3 | ~/project/methyl-seek-main/CpG/S6.bismark_bt2_pe.deduplicated.CpG_report.txt.gz  |

#### Config file

Edit the `preconfig.yaml` file as below:

- Set `rawdata_dir` to the absolute path of the directory containing all your fastqs.
- Set `result_dir` to the absolute path of the working directory containing `methyl-seek` pipeline, where all result files will be stored.
- Set `samples` to be the absolute path of your `samples.txt` located within `result_dir`.

For example,

```
samples: "~/project/methyl-seek-main/samples.txt"
rawdata_dir: "~/project/fastq/"
result_dir: "~/project/methyl-seek-main/"
```
See below instructions for adding reference to `preconfig.yaml` file.

## 4. Custom reference

By default, bisulphite marked human hg38 reference genome is provide with this repo. To create a custom bisulphite marked reference for different genome of interest, simply provide genome fasta file (for e.g.`hg19.genome.fa`), and update the information in `preconfig.yaml` as below:

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

## run : generate CpG reports
sbatch ~/project/methyl-seek-main/pipeline_launch.sh run npr ~/project/methyl-seek-main/

## dcv : perform CpG deconvolution
sbatch ~/project/methyl-seek-main/pipeline_launch.sh dcv npr ~/project/methyl-seek-main/

## dmr : perform CpG deconvolution
sbatch ~/project/methyl-seek-main/pipeline_launch.sh dmr npr ~/project/methyl-seek-main/ group1 group2

# Step 2.) To launch pipeline
module purge
module load singularity snakemake

## run : generate CpG reports
sbatch ~/project/methyl-seek-main/pipeline_submit.sh run process ~/project/methyl-seek-main/

## dcv : perform CpG deconvolution
sbatch ~/project/methyl-seek-main/pipeline_submit.sh dcv process ~/project/methyl-seek-main/

## dmr : perform CpG deconvolution
sbatch ~/project/methyl-seek-main/pipeline_submit.sh dmr process ~/project/methyl-seek-main/ group1 group2

```
