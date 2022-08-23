# <code>methyl-seek</code>

## 1. About
The `methyl-seek` executable is composed of several inter-related pipelines.

This part of the documentation describes options and concepts for ...

## 2. Example
```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
TODO: ADD STEPS HERE...
```
