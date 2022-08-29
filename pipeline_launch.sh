#! /bin/bash
set -e

###################
#
# Launching shell script for NHLBI processing of methylation data
#
###################
module load python/3.7
module load snakemake/5.13.0
module load R

cd $SLURM_SUBMIT_DIR

##
## Test commandline arguments
##
if [ $1 != "run" ] && [ $1 != "dcv" ] && [ $1 != "dmr" ] ; then
    echo " "
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: run or dcv or dmr"
    echo " "
    exit
fi

if [ $2 != "npr" ] && [ $2 != "process" ] ; then
    echo " "
    echo "Invalid commandline option: $2"
    echo "Valid commandline options include: npr or process"
    echo " "
    exit
fi

R=$3
echo $R

mkdir -p $R/snakejobs
mkdir -p $R/reports

touch $R/contrasts.txt

if [ $1 == "dmr" ]
then
  Rscript $R/scripts/create_contrasts_file.R $R $4 $5
  #### Rscript $R/scripts/create_contrasts_file.R $R group1 group2
fi

cp $R/preconfig.yaml $R/config.yaml
echo "Mode: \"$1\"" >> $R/config.yaml

# and be sure this directory and all subdirectories are writable
#chmod -fR g+rwx .
#chmod +rwx .

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

echo $SCRIPT
echo $SCRIPTPATH

##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out"

if [ $2 == "npr" ]
then
    snakemake --unlock --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
    snakemake -npr --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
fi

if [ $2 == "process" ]
then
    snakemake --unlock --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
    snakemake --latency-wait 120  -s $R/Snakefile -d $R --printshellcmds --configfile $R/config.yaml --cluster-config $R/cluster.json --keep-going --restart-times 1 --cluster "$CLUSTER_OPTS" -j 500 --rerun-incomplete --stats $R/reports/snakemake.stats | tee -a $R/reports/snakemake.log
fi
