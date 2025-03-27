# <code>methyl-seek <b>install</b></code>

## 1. About
 
The `methyl-seek` executable is composed of several inter-related sub commands. Please see `methyl-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>methyl-seek <b>install</b></code> sub command in more detail.

This page is still under construction 👷, more information is coming soon!

<!--
The `methyl-seek` executable is composed of several inter-related sub commands. Please see `methyl-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>methyl-seek <b>install</b></code> sub command in more detail.

With minimal configuration, the **`install`** sub command enables you to download the pipeline's resource bundle locally. This is necessary when setting up the pipeline on a new target system or cluster. 

The pipeline uses a set of reference files to process the data. These reference files are required and need to be available on the local file system prior to execution. This command can be used to download any required reference files of the pipeline. 

Since most resource bundles are very large; we recommend using multiple threads for pulling reference files concurrently. The resource bundle can be very large so please ensure you have sufficent disk space prior to running this sub command.

**Please Note:** The resource bundle requires about X GB of available disk space. If you are running the pipeline on the Biowulf cluster, you do *NOT* need to download the pipeline's resource bundle. It is already accessible to all HPC users. This sub command is for users running the pipeline outside of the Biowulf cluster.

Downloading the resource bundle is fast and easy! In its most basic form, <code>methyl-seek <b>install</b></code> only has *one required input*.

## 2. Synopsis
```text
$ methyl-seek install [--help] [--dry-run] \
     [--force] [--threads] \
     --ref-path REF_PATH
```

The synopsis for each command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide a output directory for the reference file download  via the `--ref-path` argument. Once the download of the resource bundle  has completed, a new child directory called methyl-seek will be created. This new directory will contain all of the pipeline's required reference files. The path to this new directory can be passed to the `--resource-bundle` option of the <code>methyl-seek <b>run</b></code> subcomand. This allow users outside of Biowulf to run the pipeline.

Use you can always use the `-h` option for information on a specific command.

### 2.1 Required Arguments

`--ref-path REF_PATH` 
 
> **Path where the resource bundle will be downloaded.**  
> *type: path*
> 
> Any resouces defined in the 'config/install.json' will be pulled onto the local filesystem. After the files have been downloaded, a new directory with the name `methyl-seek` will be created. It contains all the required reference files of the pipeline. The path to this new directory can be passed to the run sub command's `--resource-bundle` option. Please see the run sub command for more information.
> 
> ***Example:*** `--ref-path /data/$USER/refs`

### 2.2 Options

Each of the following arguments are optional and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

---  
  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what remote resources would be pulled. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--force`            
> **Force downloads all files.**  
> *type: boolean flag*
> 
> By default, any files that do not exist locally are pulled; however if a previous instance of an install did not exit gracefully, it may be necessary to forcefully re-download all the files.
>
> ***Example:*** `--force`

---  
  `--threads`            
> **Number of threads to use for concurrent file downloads.**  
> *type: int*  
> *default: 2*  
> 
> Max number of threads to use for concurrent file downloads.
>
> ***Example:*** `--threads 12`

## 3. Example
```bash 
# Step 0.) Grab an interactive node,
# do not run on head node! 
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=24gb  --cpus-per-task=12 --pty bash
module purge
module load singularity snakemake

# Step 1.) Dry-run download of the resource bundle
  methyl-seek install --ref-path /data/$USER/refs \
             --force \
             --dry-run \
             --threads 12

# Step 2.) Download the resource bundle,
# This command will NOT automatically submit
# a job to the cluster. As so, we recommend 
# submitting this next command to the cluster
# as a job. Download speeds will vary so it 
# is best to set the wall time to 2 days.
methyl-seek install --ref-path /data/$USER/refs \
           --force \
           --threads 12

# Checkout the downloaded files
cd /data/$USER/refs
tree methyl-seek
```
-->