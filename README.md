
#### Setting up the working environment

Before running the pipeline, certain packages are required to be installed within a custom conda environment.

```
module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda create --name meth
conda activate meth
mamba install -yc bioconda bwameth methyldackel
conda deactivate meth
```

#### Setting up the working files

Alter the config.yaml so the rawdata_dir is the absolute path of the directory containing all your fastqs.
Alter the result_dir so it is the absolute path of the working directory containing your snakemake pipeline, where results will be stored.
Alter samples in config.yaml to be the absolute path of your samples.txt. Check this is correct. The samples file should have the following tab-delimited format:

```
sample  group comp
S1  GA  GAvsGB
S2  GA  GAvsGB
S3  GB  GAvsGB
S4  GB  GAvsGB
S5  GC  GAvsGC
S6  GC  GAvsGC
S1  GA  GAvsGC
S2  GA  GAvsGC
```

Where GA, GB & GC are the groups these 6 samples belong to.

Within MethylSnake.sh, alter the R variable to the absolute path of your working directory.

The pipeline is divided into 4 steps:

 * bismark - which performs quality control steps and maps the data using Bismark, before extracting CpG profiles using MethylDackel.
 * bismark - which performs quality control steps and maps the data using BWA-Meth instead of Bismark, before extracting CpG profiles using MethylDackel.
 * dmr - which uses the previously generated CpG profiles to identify differentially methylated regions between groups.
 * dcv - which uses the previously generated CpG to deconvolute the data and identify which tissues samples belong to based on methylation profiles

#### Dry run of the pipeline

To perform a dry run a step of the pipeline, choose a step (e.g. bismark) and submit:

```
sh MethylSnake.sh bismark npr
```

#### Actual run of the pipeline

Once everything seems to work, to perform a full run of a step of the pipeline, submit:

```
sbatch --partition=norm --gres=lscratch:500 --time=10-00:00:00 --mail-type=BEGIN,END,FAIL MethylSnake.sh bismark process
```
