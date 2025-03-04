# Standard library imports 
from os.path import join
# 3rd party imports from pypi
from snakemake.io import expand
import pandas as pd
# Local imports
from scripts.common import (
    allocated,
    provided
)


# TODO: Add new scripts
# Differential methylation related rules
rule mvp_plots:
  input:
    bed=expand(join(dmr_dir, "bsseq", "{{group}}_{chr}_betas_pval.bed"), chr=CHRS),
  output:
    miami=join(dmr_dir, "combP", "{group}.miami.png"),
    violin=join(dmr_dir, "combP", "{group}.violin.png"),
  params:
    rname="mvp_plots",
    outdir=join(dmr_dir, "combP"),
    fdr_cutoff="0.05",
    highlight="20",
    bsseq_dir=join(dmr_dir, "bsseq"),
  resources:
    mem       = allocated("mem",       "mvp_plots", cluster),
    gres      = allocated("gres",      "mvp_plots", cluster),
    time      = allocated("time",      "mvp_plots", cluster),
    partition = allocated("partition", "mvp_plots", cluster),
  threads:
    int(allocated("threads", "mvp_plots", cluster)),
  shell:
    """
    module load R/4.1
    mkdir -p {params.outdir}
    Rscript {params.script_dir}/mvp_plots.R {params.bsseq_dir} {params.fdr_cutoff} {params.highlight} {output.miami} {output.violin}
    """


rule manhatten:
  input:
    expand(join(dmr_dir, "bsseq", "{group}_{chr}_betas_pval.bed"),chr=CHRS,group=COMPARISONS),
  output:
    LIST=temp(join(dmr_dir,"figs", "{group}.list")),
    MAN=join(dmr_dir,"figs", "{group}_manhatten_beta.png"),
  params:
    dir=join(dmr_dir, "figs"),
    rname="manhatten",
    script_dir=join(dmr_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "manhatten", cluster),
    gres      = allocated("gres",      "manhatten", cluster),
    time      = allocated("time",      "manhatten", cluster),
    partition = allocated("partition", "manhatten", cluster),
  threads:
    int(allocated("threads", "manhatten", cluster)),
  shell:
    """
    mkdir -p {params.dir}
    ls {params.dir}/{wildcards.group}_*_betas_pval.bed > {output.LIST}
    module load R
    Rscript {params.script_dir}/bsseq_Manhatten.R {output.LIST} {output.MAN}
    """