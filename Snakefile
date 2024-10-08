#####################################################################################################
# Bisulphite sequencing (methyl-seek) analysis workflow
#
# Execution mode: run and dmr
# Last Modified: August, 2022
#
#####################################################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards
from os import listdir
import pandas as pd
import json


# Helper functions 
def allocated(resource, rule, lookup, default="__default__"):
    """Pulls resource information for a given rule. If a rule does not have any information 
    for a given resource type, then it will pull from the default. Information is pulled from
    definitions in the cluster.json (which is used a job submission). This ensures that any 
    resources used at runtime mirror the resources that were allocated.
    :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
    :param rule <str>: rule to lookup its information
    :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
    :param default <str>: default information to use if rule information cannot be found
    :return allocation <str>: 
        allocation information for a given resource type for a given rule
    """

    try: 
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]
    
    return allocation


# Global workflow variables
sample_file = config["samples"]
rawdata_dir = config["rawdata_dir"]
working_dir = config["result_dir"]

# References
hg38_fa                = config["hg38_fa"]
phage_fa               = config["phage_fa"]
hg38_gtf               = config["hg38_gtf"]
hg38_rRNA_intervals    = config["hg38_rRNA_intervals"]
hg38_bed_ref           = config["hg38_bed_ref"]
hg38_refFlat           = config["hg38_refFlat"]
bisulphite_genome_path = config["bisulphite_ref"]
phage_genome_path      = config["phage_ref"]
bisulphite_fa          = config["bisulphite_fa"]
species                = config["species"]
REF_ATLAS              = config["REF_ATLAS"]
CpG_MAP_TABLE          = config["CpG_MAP_TABLE"]
wgbsDir                = config["wgbsDir"]
mode                   = config["Mode"]

REF_MARKERS            = config["REF_MARKERS"]
REF_IDS                = config["REF_IDS"]
GOLD_MARKERS           = config["GOLD_MARKERS"]
HG38TOHG19             = config["HG38TOHG19"]


if species == "hg38":
    atlas_bed = config["atlas_bed_hg38"]
    atlas_tsv = config["atlas_tsv_hg38"]

if species == "hg19":
    atlas_bed = config["atlas_bed_hg19"]
    atlas_tsv = config["atlas_tsv_hg19"]

# sample list
df      = pd.read_csv(sample_file, header=0, sep='\t')
SAMPLES = list(set(df['samples'].tolist()))
GROUPS  = list(set(df['group'].tolist()))

# Read in resource information,
# containing information about 
# threads, mem, walltimes, etc.
with open(join('cluster.json')) as fh:
    cluster = json.load(fh)

CHRS = [
  'chr1', 'chr2', 'chr3', 
  'chr4', 'chr5', 'chr6', 
  'chr7', 'chr8', 'chr9',
  'chr10', 'chr11', 'chr12',
  'chr13', 'chr14', 'chr15', 
  'chr16', 'chr17', 'chr18', 
  'chr19', 'chr20', 'chr21', 
  'chr22', 'chr23', 'chrX'
]

RN = ['R1', 'R2']


dmr_dir = working_dir
if mode == "dmr":
  # contrasts list
  df2 = pd.read_csv(working_dir + "/contrasts.txt", header=0, sep='\t')
  GROUPS=list(set(df2['comparisons'].tolist()))
  GRPSAMPLES=list(set(df2['samples'].tolist()))
  dmr_dir = join(working_dir, "dmr", GROUPS[0])


def output_from_modes():
    outputs = list()

    if mode == "run":
        outputs.append(join(working_dir, "multiqc_report.html"))
    if mode == "dmr":
        outputs.append(join(working_dir, "multiqc_report.html"))

        #homer files
        outputs.append(expand(join(dmr_dir, "homer/{group}_{chr}.homerOutput.txt"),group=GROUPS,chr=CHRS))
        outputs.append(expand(join(dmr_dir, "homer/{group}_{chr}.homerOutput2.txt"),group=GROUPS,chr=CHRS))

        # lm methylation sites from bismark alignments
        outputs.append(expand(join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),chr=CHRS,group=GROUPS))
        #outputs.append(expand(join(dmr_dir, "bsseq/{group}_{chr}_betas_wtFDR_unfiltered.bed"),chr=CHRS,group=GROUPS))
        outputs.append(expand(join(dmr_dir, "combP/{group}_{chr}.regions-p.bed.gz"),chr=CHRS,group=GROUPS))

        # summary figures of analysis
        outputs.append(expand(join(dmr_dir,"figs/{group}_manhatten_beta.png"),group=GROUPS))

        # combp summary of results
        outputs.append(expand(join(dmr_dir, "combP/{group}.regions-p.bed"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "figs/{group}_GREATfig1.png"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "figs/{group}_manhatten.png"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "combP/{group}.GREATprocesses.txt"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "combP/{group}.GREATfunction.txt"),group=GROUPS))

        #MVP mvp_plots
        outputs.append(expand(join(dmr_dir, "combP/{group}.miami.png"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "combP/{group}.violin.png"),group=GROUPS))
        #outputs.append(expand(join(dmr_dir, "combP/{group}.pie-cpg.png"),group=GROUPS))
        #outputs.append(expand(join(dmr_dir, "combP/{group}.pie-genic.png"),group=GROUPS))

    if mode == "dcv":
        outputs.append(join(working_dir, "multiqc_report.html"))
        outputs.append(expand(join(working_dir, "CpG_CSV/{samples}.csv"),samples=SAMPLES))
        outputs.append(expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES))
        outputs.append(expand(join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),samples=SAMPLES))
        outputs.append(join(working_dir, "deconvolution_CSV/total.csv"))
        outputs.append(join(working_dir, "deconvolution_CSV/total_deconv_output.csv"))
        outputs.append(join(working_dir, "deconvolution_CSV/total_deconv_plot.png"))
        #outputs.append(expand(join(working_dir, "bismarkAlign/{samples}.bismark_pe.deduplicated.bam"),samples=SAMPLES))
        outputs.append(expand(join(working_dir,"CpG/{samples}/{samples}.cfDNAmeInput.bedGraph"),samples=SAMPLES))
        outputs.append(expand(join(working_dir,"CpG/{samples}/{samples}.cfDNAmeDeconvolution.tsv"),samples=SAMPLES))
        outputs.append(expand(join(working_dir,"UXM/{samples}.pat.gz"), samples=SAMPLES))
        outputs.append(expand(join(working_dir,"UXM/{samples}_deconv.250.csv"), samples=SAMPLES))
    return(outputs)


rule All:
    input:
      # Creating data links:
      expand(join(working_dir, "raw/{samples}.{rn}.fastq.gz"), samples=SAMPLES, rn=RN),

      # Checking data quality:
      expand(join(working_dir, "rawQC/{samples}.{rn}_fastqc.html"), samples=SAMPLES, rn=RN),

      # Quality trimming output:
      expand(join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),samples=SAMPLES),
      expand(join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),samples=SAMPLES),

      # kraken output
      #expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.out.txt"),samples=SAMPLES),
      #expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.taxa.txt"),samples=SAMPLES),
      #expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.krona.html"),samples=SAMPLES),

      #FQscreen output
      #expand(join(working_dir, "FQscreen","{samples}_val_1_screen.txt"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen","{samples}_val_1_screen.png"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen","{samples}_val_2_screen.txt"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen","{samples}_val_2_screen.png"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen2","{samples}_val_1_screen.txt"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen2","{samples}_val_1_screen.png"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen2","{samples}_val_2_screen.txt"),samples=SAMPLES),
      #expand(join(working_dir, "FQscreen2","{samples}_val_2_screen.png"), samples=SAMPLES),

      # bisulphite genome preparation
      join(bisulphite_genome_path, species, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(bisulphite_genome_path, species, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),

      # bismark align to human reference genomes
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.flagstat"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),samples=SAMPLES),

      # get alignment statistics
      expand(join(working_dir, "bismarkAlign/{samples}.RnaSeqMetrics.txt"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.flagstat.concord.txt"),samples=SAMPLES),
      #expand(join(working_dir, "rseqc/{samples}.inner_distance_freq.txt"),samples=SAMPLES),
      #expand(join(working_dir, "rseqc/{samples}.strand.info"),samples=SAMPLES),
      #expand(join(working_dir, "rseqc/{samples}.Rdist.info"),samples=SAMPLES),
      #expand(join(working_dir, "trimGalore/{samples}_insert_sizes.txt"),samples=SAMPLES),
      expand(join(working_dir, "CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),samples=SAMPLES),
      expand(join(working_dir, "CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),samples=SAMPLES),


      join(working_dir, "bismark_summary_report.txt"),
      join(working_dir, "bismark_summary_report.html"),
      output_from_modes()

      #########  generate multiqc output
      #join(working_dir, "multiqc_report.html"),

      ######## Following lines are applicable only when "dcv" mode is enabled ###############
      # Deconvolution output
      #expand(join(working_dir, "CpG_CSV/{samples}.csv"),samples=SAMPLES),
      #expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
      #expand(join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),samples=SAMPLES),
      #join(working_dir, "deconvolution_CSV/total.csv"),
      #join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
      #join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),

## Copy raw data to working directory
rule raw_data_links:
    input:
      join(rawdata_dir, "{samples}.{rn}.fastq.gz")
    output:
      join(working_dir, "raw/{samples}.{rn}.fastq.gz")
    params:
      rname="raw_data_links",
      dir=directory(join(working_dir, "raw")),
    resources:
      mem       = allocated("mem",       "raw_data_links", cluster),
      gres      = allocated("gres",      "raw_data_links", cluster),
      time      = allocated("time",      "raw_data_links", cluster),
      partition = allocated("partition", "raw_data_links", cluster),
    threads:
      int(allocated("threads", "raw_data_links", cluster)),
    shell:
      """
      mkdir -p {params.dir}
      ln -s {input} {output}
      """


## Run fastqc on raw data to visually assess quality
rule raw_fastqc:
    input:
      join(working_dir, "raw/{samples}.{rn}.fastq.gz")
    output:
      join(working_dir, "rawQC/{samples}.{rn}_fastqc.html")
    params:
      rname="raw_fastqc",
      dir=directory(join(working_dir, "rawQC")),
      batch='--cpus-per-task=2 --mem=8g --time=8:00:00',
    resources:
      mem       = allocated("mem",       "raw_fastqc", cluster),
      gres      = allocated("gres",      "raw_fastqc", cluster),
      time      = allocated("time",      "raw_fastqc", cluster),
      partition = allocated("partition", "raw_fastqc", cluster),
    threads:
      int(allocated("threads", "raw_fastqc", cluster)),
    shell:
      """
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
      """


## trimming/filtering with trimGalore
rule trimGalore:
    input:
      F1=join(working_dir, "raw/{samples}.R1.fastq.gz"),
      F2=join(working_dir, "raw/{samples}.R2.fastq.gz"),
    output:
      join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      join(working_dir, "trimGalore/{samples}_val_2.fq.gz")
    params:
      rname="trimGalore",
      dir=directory(join(working_dir, "trimGalore")),
      tag='{samples}',
      fastqcdir=directory(join(working_dir, "postTrimQC")),
      command="--fastqc --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --length 50 --gzip",
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=25g --time=10:00:00',
      workdir = working_dir,
    resources:
      mem       = allocated("mem",       "trimGalore", cluster),
      gres      = allocated("gres",      "trimGalore", cluster),
      time      = allocated("time",      "trimGalore", cluster),
      partition = allocated("partition", "trimGalore", cluster),
    threads:
      int(allocated("threads", "trimGalore", cluster)),
    shell:
      """
      module load trimgalore/0.6.7
      mkdir -p {params.dir}
      trim_galore --paired --cores {threads} {params.command} --basename {params.tag} --output_dir {params.dir} --fastqc_args "--outdir {params.fastqcdir}"  {input.F1} {input.F2}
      """


rule bbmerge:
    input:
      R1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      R2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      join(working_dir, "trimGalore/{samples}_insert_sizes.txt"),
    params:
      rname='pl:bbmerge',
      script_dir=join(working_dir,"scripts"),
    resources:
      mem       = allocated("mem",       "bbmerge", cluster),
      gres      = allocated("gres",      "bbmerge", cluster),
      time      = allocated("time",      "bbmerge", cluster),
      partition = allocated("partition", "bbmerge", cluster),
    threads:
      int(allocated("threads", "bbmerge", cluster)),
    shell: 
      """
      # Get encoding of Phred Quality Scores
      module load python
      encoding=$(python {params.script_dir}/phred_encoding.py {input.R1})
      echo "Detected Phred+${{encoding}} ASCII encoding"

      module load bbtools/38.87
      bbtools bbmerge-auto in1={input.R1} in2={input.R2} qin=${{encoding}} \
      ihist={output} k=62 extend2=200 rem ecct -Xmx900G
      """


rule fastq_screen:
    input:
      file1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      file2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      out1=join(working_dir,"FQscreen","{samples}_val_1_screen.txt"),
      out2=join(working_dir,"FQscreen","{samples}_val_1_screen.png"),
      out3=join(working_dir,"FQscreen","{samples}_val_2_screen.txt"),
      out4=join(working_dir,"FQscreen","{samples}_val_2_screen.png"),
      out5=join(working_dir,"FQscreen2","{samples}_val_1_screen.txt"),
      out6=join(working_dir,"FQscreen2","{samples}_val_1_screen.png"),
      out7=join(working_dir,"FQscreen2","{samples}_val_2_screen.txt"),
      out8=join(working_dir,"FQscreen2","{samples}_val_2_screen.png")
    params:
      rname='pl:fqscreen',
      outdir = join(working_dir,"FQscreen"),
      outdir2 = join(working_dir,"FQscreen2"),
      fastq_screen_config="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf",
      fastq_screen_config2="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen_2.conf",
    resources:
      mem       = allocated("mem",       "fastq_screen", cluster),
      gres      = allocated("gres",      "fastq_screen", cluster),
      time      = allocated("time",      "fastq_screen", cluster),
      partition = allocated("partition", "fastq_screen", cluster),
    threads:
      int(allocated("threads", "fastq_screen", cluster)),
    shell:
      """
      module load fastq_screen/0.14.1
      module load bowtie/2-2.3.4
      module load perl/5.24.3
      fastq_screen --conf {params.fastq_screen_config} --outdir {params.outdir} \
        --threads {threads} --subset 1000000 \
        --aligner bowtie2 --force {input.file1} {input.file2}
      fastq_screen --conf {params.fastq_screen_config2} --outdir {params.outdir2} \
        --threads {threads} --subset 1000000 \
        --aligner bowtie2 --force {input.file1} {input.file2}
      """


rule kraken_pe:
    input:
        fq1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
        fq2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
        krakenout = join(working_dir, "kraken","{samples}.trim.kraken_bacteria.out.txt"),
        krakentaxa = join(working_dir, "kraken","{samples}.trim.kraken_bacteria.taxa.txt"),
        kronahtml = join(working_dir, "kraken","{samples}.trim.kraken_bacteria.krona.html"),
    params:
        rname='pl:kraken',
        dir=join(working_dir,"kraken"),
        #bacdb="/fdb/kraken/20170202_bacteria",
        bacdb="/fdb/kraken/20210223_standard_kraken2",
        prefix="{samples}",
    resources:
      mem       = allocated("mem",       "kraken_pe", cluster),
      gres      = allocated("gres",      "kraken_pe", cluster),
      time      = allocated("time",      "kraken_pe", cluster),
      partition = allocated("partition", "kraken_pe", cluster),
    threads:
      int(allocated("threads", "kraken_pe", cluster)),
    shell:
      """
      module load kraken
      module load kronatools/2.7
      mkdir -p {params.dir}
      cd /lscratch/$SLURM_JOBID;
      cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;

      kdb_base=$(basename {params.bacdb})
      kraken2 --db /lscratch/$SLURM_JOBID/${{kdb_base}} \
        --threads {threads} --report {output.krakentaxa} \
        --output {output.krakenout} \
        --gzip-compressed \
        --paired {input.fq1} {input.fq2}
      # Generate Krona Report
      cut -f2,3 {output.krakenout} | ktImportTaxonomy - -o {output.kronahtml}

      #kdb_base=$(basename {params.bacdb})
      #kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload--paired {input.fq1} {input.fq2}
      #kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
      #cut -f 2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
      """


######## NEED TO PUT CONDITION HERE:
# If user provided different genome version, then run this "prep_bisulphite_genome" rule:
#
##########
rule prep_bisulphite_genome:
    input:
      bisulphite_fa
    output:
      join(bisulphite_genome_path, species, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(bisulphite_genome_path, species, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),
    params:
      rname="prep_bisulphite_genome",
      dir=directory(join(bisulphite_genome_path, species)),
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=20g --time=2:00:00',
    resources:
      mem       = allocated("mem",       "prep_bisulphite_genome", cluster),
      gres      = allocated("gres",      "prep_bisulphite_genome", cluster),
      time      = allocated("time",      "prep_bisulphite_genome", cluster),
      partition = allocated("partition", "prep_bisulphite_genome", cluster),
    threads:
      int(allocated("threads", "prep_bisulphite_genome", cluster)),
    shell:
      """
      module load bismark/0.23.0
      mkdir -p {params.dir}
      cd {params.dir}
      #cp reference_genome {params.dir}/genome.fa
      bismark_genome_preparation --verbose --parallel {threads} {params.dir} #--single_fasta
      """


rule bismark_align:
    input:
      F1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      F2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      B1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),
      B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.flagstat"),
      FQ1=temp(join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_unmapped_reads_1.fq.gz")),
      FQ2=temp(join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_unmapped_reads_2.fq.gz")),
      FQ3=temp(join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_ambiguous_reads_1.fq.gz")),
      FQ4=temp(join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_ambiguous_reads_2.fq.gz")),
    params:
      rname="bismark_align",
      dir=directory(join(working_dir, "bismarkAlign")),
      genome_dir=directory(join(bisulphite_genome_path, species)),
      command="--bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 4 --score_min L,-0.6,-0.6",
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=100g --time=10:00:00',
      outbam=join(working_dir, "bismarkAlign/{samples}_val_1_bismark_bt2_pe.bam"),
      R1=join(working_dir,"bismarkAlign/{samples}_val_1_bismark_bt2_PE_report.txt"),
      R2=join(working_dir,"bismarkAlign/{samples}.bismark_bt2_PE_report.txt"),
    resources:
      mem       = allocated("mem",       "bismark_align", cluster),
      gres      = allocated("gres",      "bismark_align", cluster),
      time      = allocated("time",      "bismark_align", cluster),
      partition = allocated("partition", "bismark_align", cluster),
    threads:
      int(allocated("threads", "bismark_align", cluster)),
    shell:
      """
      module load bismark/0.23.0 samtools
      mkdir -p {params.dir}
      bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
      mv {params.outbam} {output.B1}
      mv {params.R1} {params.R2}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """


rule bismark_dedup:
    input:
      F1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),
    output:
      T1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.deduplicated.bam")),
      B1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam")),
      B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),
      C1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),
    params:
      rname="bismark_dedup",
      dir=directory(join(working_dir, "bismarkAlign")),
      prefix=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe"),
      genome=bisulphite_fa,
    resources:
      mem       = allocated("mem",       "bismark_dedup", cluster),
      gres      = allocated("gres",      "bismark_dedup", cluster),
      time      = allocated("time",      "bismark_dedup", cluster),
      partition = allocated("partition", "bismark_dedup", cluster),
    threads:
      int(allocated("threads", "bismark_dedup", cluster)),
    shell:
      """
      module load bismark/0.23.0
      module load samtools/1.15
      cd {params.dir}
      samtools view -hb {input.F1} | samtools sort -n -@ {threads} -O BAM -o {output.T1}
      deduplicate_bismark --paired --bam --outfile {params.prefix} {output.T1}
      samtools view -T {params.genome} -@ 16 -h -C {output.B1} > {output.C1}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """


rule bismark_extract:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    cov=join(working_dir, "CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
    bed=join(working_dir,"CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
  params:
    rname='pl:bismark_extract',
    bismark_index=join(bisulphite_genome_path,species),
    outdir=join(working_dir, "CpG/{samples}"),
    genome=bisulphite_fa,
  resources:
    mem       = lambda wc, attempt: allocated("mem",       "bismark_extract", cluster) if attempt == 1 else "900g",
    gres      = allocated("gres",      "bismark_extract", cluster),
    time      = allocated("time",      "bismark_extract", cluster),
    partition = lambda wc, attempt: allocated("partition", "bismark_extract", cluster) if attempt == 1 else "largemem",
  threads:
    int(allocated("threads", "bismark_extract", cluster)),
  shell:
    """
    mkdir -p {params.outdir}
    module load bismark samtools bowtie
    bismark_methylation_extractor --paired-end --no_overlap --multicore 8 --gzip --report --bedGraph --counts --buffer_size '50%' --no_header --cytosine_report --output {params.outdir} --scaffolds --genome_folder {params.bismark_index} {input.bam}
    """


rule prep_bisulphite_phage_genome:
    input:
      phage_fa
    output:
      join(phage_genome_path, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(phage_genome_path, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),
    params:
      rname="prep_phage_genome",
      dir=directory(phage_genome_path),
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=20g --time=2:00:00',
    resources:
      mem       = allocated("mem",       "prep_bisulphite_phage_genome", cluster),
      gres      = allocated("gres",      "prep_bisulphite_phage_genome", cluster),
      time      = allocated("time",      "prep_bisulphite_phage_genome", cluster),
      partition = allocated("partition", "prep_bisulphite_phage_genome", cluster),
    threads:
      int(allocated("threads", "prep_bisulphite_phage_genome", cluster)),
    shell:
      """
      module load bismark/0.23.0
      mkdir -p {params.dir}
      cd {params.dir}
      cp reference_genome {params.dir}/genome.fa
      bismark_genome_preparation --verbose --parallel {threads} {params.dir} #--single_fasta
      """


rule bismark_phage:
  input:
    F1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
    F2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
  output:
    join(working_dir, "bismark_phage/{samples}_val_1_bismark_bt2_pe.bam"),
  params:
    rname="bismark_phage",
    dir=directory(join(working_dir, "bismark_phage")),
    genome_dir=directory(phage_genome_path),
    command="--bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 2 --score_min L,-0.6,-0.6",
    batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=100g --time=10:00:00',
  resources:
    mem       = allocated("mem",       "bismark_phage", cluster),
    gres      = allocated("gres",      "bismark_phage", cluster),
    time      = allocated("time",      "bismark_phage", cluster),
    partition = allocated("partition", "bismark_phage", cluster),
  threads:
    int(allocated("threads", "bismark_phage", cluster)),
  shell:
    """
    module load bismark/0.23.0
    mkdir -p {params.dir}
    bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
    """


rule rseqc:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    out1=join(working_dir,"rseqc/{samples}.strand.info"),
    out2=join(working_dir,"rseqc/{samples}.Rdist.info")
  params:
    bedref=hg38_bed_ref,
    dir=join(working_dir,"rseqc"),
    rname="pl:rseqc",
  resources:
    mem       = allocated("mem",       "rseqc", cluster),
    gres      = allocated("gres",      "rseqc", cluster),
    time      = allocated("time",      "rseqc", cluster),
    partition = allocated("partition", "rseqc", cluster),
  threads:
    int(allocated("threads", "rseqc", cluster)),
  shell:
    """
    module load rseqc/4.0.0
    mkdir -p {params.dir}
    infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
    read_distribution.py -i {input.file1} -r {params.bedref} > {output.out2}
    """


rule inner_distance:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    innerdists=join(working_dir,"rseqc/{samples}.inner_distance_freq.txt"),
  params:
    prefix=join(working_dir,"rseqc/{samples}"),
    dir=join(working_dir,"rseqc"),
    genemodel=hg38_bed_ref,
    rname="pl:inner_distance",
  resources:
    mem       = allocated("mem",       "inner_distance", cluster),
    gres      = allocated("gres",      "inner_distance", cluster),
    time      = allocated("time",      "inner_distance", cluster),
    partition = allocated("partition", "inner_distance", cluster),
  threads:
    int(allocated("threads", "inner_distance", cluster)),
  shell:
    """
    module load rseqc/4.0.0
    mkdir -p {params.dir}
    inner_distance.py -i {input.bam} -r {params.genemodel} -k 10000000 -o {params.prefix}
    """

rule stats:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    outstar1=join(working_dir,"bismarkAlign/{samples}.RnaSeqMetrics.txt"),
    outstar2=join(working_dir,"bismarkAlign/{samples}.flagstat.concord.txt"),
  params:
    rname='pl:stats',
    refflat=hg38_refFlat,
    rrnalist=hg38_rRNA_intervals,
    picardstrand="SECOND_READ_TRANSCRIPTION_STRAND",
    statscript=join("workflow", "scripts", "bam_count_concord_stats.py"),
  resources:
    mem       = allocated("mem",       "stats", cluster),
    gres      = allocated("gres",      "stats", cluster),
    time      = allocated("time",      "stats", cluster),
    partition = allocated("partition", "stats", cluster),
  threads:
    int(allocated("threads", "stats", cluster)),
  shell:
    """
    module load python/3.8 samtools/1.15 picard/2.27.3
    java -Xmx110g -jar ${{PICARDJARPATH}}/picard.jar CollectRnaSeqMetrics REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist} STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    samtools flagstat {input.file1} > {output.outstar2};
    ## python3 {params.statscript} {input.file1} >> {output.outstar2} ## does require "Rstat" not available in biowulf
    """


rule multiqc:
  input:
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),samples=SAMPLES),
  output:
    join(working_dir, "bismark_summary_report.txt"),
    join(working_dir, "bismark_summary_report.html"),
    join(working_dir, "multiqc_report.html"),
  params:
    rname="multiqc",
    dir=working_dir,
    bis_dir=directory(join(working_dir,"bismarkAlign")),
    script_dir=join(working_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "multiqc", cluster),
    gres      = allocated("gres",      "multiqc", cluster),
    time      = allocated("time",      "multiqc", cluster),
    partition = allocated("partition", "multiqc", cluster),
  threads:
    int(allocated("threads", "multiqc", cluster)),
  shell:
    """
    module load multiqc/1.9 bismark
    cd {params.bis_dir}
    bismark2report
    bismark2summary
    mv bismark_summary_report.html {params.dir}/
    mv bismark_summary_report.txt {params.dir}/
    cd {params.dir}
    multiqc --ignore '*/.singularity/*' -f --interactive .
    """

############### 
# Deconvolution rules begin here
# run following rules in dcv mode only
#################
rule get_CpG:
  input:
    join(working_dir, "CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
  output:
    join(working_dir, "CpG_CSV/{samples}.csv"),
  params:
    rname="get_CpG",
    cutoff=5,
    script_dir=join(working_dir,"scripts"),
    dir1=join(working_dir,"CpG_CSV"),
    dir2=join(working_dir,"deconvolution_CSV"),
  resources:
    mem       = allocated("mem",       "get_CpG", cluster),
    gres      = allocated("gres",      "get_CpG", cluster),
    time      = allocated("time",      "get_CpG", cluster),
    partition = allocated("partition", "get_CpG", cluster),
  threads:
    int(allocated("threads", "get_CpG", cluster)),
  shell:
    """
    mkdir -p {params.dir1}
    mkdir -p {params.dir2}
    module load R
    Rscript {params.script_dir}/get_methy.R {input} {wildcards.samples} {params.cutoff} {output}
    """


rule get_overlap_meth:
  input:
    join(working_dir, "CpG_CSV/{samples}.csv"),
  output:
    join(working_dir, "deconvolution_CSV/{samples}.csv"),
  params:
    rname="get_overlap_meth",
    map_table=CpG_MAP_TABLE,
  resources:
    mem       = allocated("mem",       "get_overlap_meth", cluster),
    gres      = allocated("gres",      "get_overlap_meth", cluster),
    time      = allocated("time",      "get_overlap_meth", cluster),
    partition = allocated("partition", "get_overlap_meth", cluster),
  threads:
    int(allocated("threads", "get_overlap_meth", cluster)),
  run:
    df_ref=pd.read_csv(params.map_table,sep='\t',header=None)
    df_ref.columns=['chromosome','start','end','cgid']
    df_ref=df_ref.loc[(df_ref['chromosome'].isin(CHRS)),]
    dfm=pd.read_csv(input[0])
    dfm=pd.merge(df_ref,dfm,on=['chromosome','start','end'],how='inner')
    dfm=dfm.drop(labels=['chromosome','start','end'],axis=1)
    dfm=dfm.set_index('cgid')
    dfm.to_csv(output[0])


rule run_deconv:
  input:
    join(working_dir, "deconvolution_CSV/{samples}.csv"),
  output:
    join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),
  params:
    script_dir=join(working_dir,"scripts"),
    dir=join(working_dir,"deconvolution_CSV"),
    rname="run_deconv",
    ref=REF_ATLAS,
  resources:
    mem       = allocated("mem",       "run_deconv", cluster),
    gres      = allocated("gres",      "run_deconv", cluster),
    time      = allocated("time",      "run_deconv", cluster),
    partition = allocated("partition", "run_deconv", cluster),
  threads:
    int(allocated("threads", "run_deconv", cluster)),
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve2.py --atlas_path {params.ref} --residuals {input}  > {output}  2>&1
    """


rule merge_tables:
  input:
    expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
  output:
    join(working_dir, "deconvolution_CSV/total.csv"),
  params:
    rname="merge_tables",
  resources:
    mem       = allocated("mem",       "merge_tables", cluster),
    gres      = allocated("gres",      "merge_tables", cluster),
    time      = allocated("time",      "merge_tables", cluster),
    partition = allocated("partition", "merge_tables", cluster),
  threads:
    int(allocated("threads", "merge_tables", cluster)),
  run:
    dfm=pd.read_csv(input[0])
    for f in input[1:]:
    	df=pd.read_csv(f)
    	dfm=pd.merge(dfm,df,on='cgid',how='outer')
    dfm.to_csv(output[0],index=False)


rule run_deconv_merged:
  input:
    join(working_dir, "deconvolution_CSV/total.csv"),
  output:
    join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
    join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),
  params:
    rname="run_deconv_merged",
    ref=REF_ATLAS,
    dir=join(working_dir,"deconvolution_CSV"),
    script_dir=join(working_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "run_deconv_merged", cluster),
    gres      = allocated("gres",      "run_deconv_merged", cluster),
    time      = allocated("time",      "run_deconv_merged", cluster),
    partition = allocated("partition", "run_deconv_merged", cluster),
  threads:
    int(allocated("threads", "run_deconv_merged", cluster)),
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve2.py --atlas_path {params.ref} --residuals {input}
    """


###############CFDNAME RULES######################
rule sorting_CpG_bedgraph:
  input:
    bed=join(working_dir,"CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
  output:
    bed=temp(join(working_dir,"CpG/{samples}/{samples}.mapped_autosomal_CpG.bedGraph.tmp")),
  params:
    rname="pl:format1",
  resources:
    mem       = allocated("mem",       "sorting_CpG_bedgraph", cluster),
    gres      = allocated("gres",      "sorting_CpG_bedgraph", cluster),
    time      = allocated("time",      "sorting_CpG_bedgraph", cluster),
    partition = allocated("partition", "sorting_CpG_bedgraph", cluster),
  threads:
    int(allocated("threads", "sorting_CpG_bedgraph", cluster)),
  shell:
    """
    module load bedtools
    zcat {input.bed} | tail -n +2 | bedtools sort -i - > {output.bed}
    """


rule liftover_bedgraph:
  input:
    bed=join(working_dir,"CpG/{samples}/{samples}.mapped_autosomal_CpG.bedGraph.tmp"),
  output:
    graph=temp(join(working_dir,"CpG/{samples}/{samples}.liftover.bedGraph.tmp")),
  params:
    rname="pl:format2",
    lift_file=HG38TOHG19,
  resources:
    mem       = allocated("mem",       "liftover_bedgraph", cluster),
    gres      = allocated("gres",      "liftover_bedgraph", cluster),
    time      = allocated("time",      "liftover_bedgraph", cluster),
    partition = allocated("partition", "liftover_bedgraph", cluster),
  threads:
    int(allocated("threads", "liftover_bedgraph", cluster)),
  shell:
    """
    module load bedtools crossmap
    crossmap bed {params.lift_file} {input.bed} {output.graph}
    """


rule extract_signature_beds:
  input:
    graph=join(working_dir,"CpG/{samples}/{samples}.liftover.bedGraph.tmp"),
  output:
    sort=temp(join(working_dir,"CpG/{samples}/{samples}.sorted.bedGraph")),
  params:
    rname="pl:format3",
    markers=GOLD_MARKERS,
  resources:
    mem       = allocated("mem",       "extract_signature_beds", cluster),
    gres      = allocated("gres",      "extract_signature_beds", cluster),
    time      = allocated("time",      "extract_signature_beds", cluster),
    partition = allocated("partition", "extract_signature_beds", cluster),
  threads:
    int(allocated("threads", "extract_signature_beds", cluster)),
  shell:
    """
    module load bedtools
    bedtools sort -i {input.graph} | bedtools intersect -wo -a {params.markers} -b stdin -sorted | awk '$6-$5==1 {{print $0}}' | awk 'NF{{NF-=1}};1' > {output.sort}
    """


rule aggregate_over_regions:
  input:
    sort=join(working_dir,"CpG/{samples}/{samples}.sorted.bedGraph"),
  output:
    tsv=join(working_dir,"CpG/{samples}/{samples}.cfDNAmeInput.bedGraph"),
  params:
    rname="pl:format4",
    script_dir=join(working_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "aggregate_over_regions", cluster),
    gres      = allocated("gres",      "aggregate_over_regions", cluster),
    time      = allocated("time",      "aggregate_over_regions", cluster),
    partition = allocated("partition", "aggregate_over_regions", cluster),
  threads:
    int(allocated("threads", "aggregate_over_regions", cluster)),
  shell:
    """
    module load R
    Rscript {params.script_dir}/aggregate_over_regions.R {input.sort} {output.tsv}
    """


rule cfDNAme:
  input:
    bed=join(working_dir,"CpG/{samples}/{samples}.cfDNAmeInput.bedGraph"),
  output:
    tsv=join(working_dir,"CpG/{samples}/{samples}.cfDNAmeDeconvolution.tsv"),
  params:
    rname="pl:cfDNAme",
    script_dir=join(working_dir,"scripts"),
    sampleName="{samples}",
    reference_markers=REF_MARKERS,
    reference_IDs=REF_IDS,
  resources:
    mem       = allocated("mem",       "cfDNAme", cluster),
    gres      = allocated("gres",      "cfDNAme", cluster),
    time      = allocated("time",      "cfDNAme", cluster),
    partition = allocated("partition", "cfDNAme", cluster),
  threads:
    int(allocated("threads", "cfDNAme", cluster)),
  shell:
    """
    module load R
    Rscript ${params.script_dir}/tissues_of_origin_v2.R \
    {input.bed}  \
    {params.reference_markers} \
    {output.tsv} \
    {params.reference_IDs} \
    {params.sampleName} \
    TRUE \
    FALSE \
    TRUE \
    colon_5/hsc_5/hsc_2/spleen_1/spleen_2/spleen_3/spleen_4/
    """


############### 
# Differential methylation rules begin here
# run following rules in dmr mode only
#################
rule bsseq_bismark:
  input:
    bizfile=join(working_dir,"contrasts.txt"),
  output:
    bed1=join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),
  params:
    rname="bsseq_bismark",
    chr='{chr}',
    dir=directory(join(dmr_dir, "bsseq")),
    cov="2",
    sample_prop="0.25",
    script_dir=join(working_dir, "dmr","scripts"),
  resources:
    mem       = allocated("mem",       "bsseq_bismark", cluster),
    gres      = allocated("gres",      "bsseq_bismark", cluster),
    time      = allocated("time",      "bsseq_bismark", cluster),
    partition = allocated("partition", "bsseq_bismark", cluster),
  threads:
    int(allocated("threads", "bsseq_bismark", cluster)),
  shell:
    """
    module load R
    mkdir -p {params.dir}
    Rscript {params.script_dir}/bsseq_lm.R {params.chr} {input.bizfile} {output.bed1} {params.sample_prop} {params.cov}
    """


rule combP:
  input:
    join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),
  output:
    RP=join(dmr_dir, "combP/{group}_{chr}.regions-p.bed.gz"),
    homerInput=temp(join(dmr_dir, "homer/{group}_{chr}.homerInput.txt")),
    homerOutput=join(dmr_dir, "homer/{group}_{chr}.homerOutput.txt"),
    homerOutput2=join(dmr_dir, "homer/{group}_{chr}.homerOutput2.txt"),
    homerAnn=temp(join(dmr_dir,"homer/{group}_{chr}.homer.annStats.txt")),
  params:
    rname="CombP",
    groups='{group}_{chr}',
    dir=join(dmr_dir, "combP"),
    annStat=join(dmr_dir,"scripts/blank.homer.annStats.txt"),
    homer_dir=join(dmr_dir,"homer"),
    dist="300",
    step="60",
  resources:
    mem       = allocated("mem",       "combP", cluster),
    gres      = allocated("gres",      "combP", cluster),
    time      = allocated("time",      "combP", cluster),
    partition = allocated("partition", "combP", cluster),
  threads:
    int(allocated("threads", "combP", cluster)),
  shell:
    """
    mkdir -p {params.dir}
    module load combp bedtools
    comb-p pipeline -c 4 --dist {params.dist} --step {params.step} --seed 0.01 -p {params.dir}/{params.groups} --region-filter-p 0.05 --region-filter-n 3 {input}

    module load homer
    mkdir -p {params.homer_dir}
    cp {params.annStat} {output.homerAnn}
    cat {output.RP} | sed "1,1d" | awk '{{print$1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$4"|"$5"|"$6"|"$7"\t""*"}}' > {output.homerInput}
    annotatePeaks.pl {output.homerInput} hg38 -annStats {output.homerAnn} > {output.homerOutput}
    awk 'NR==FNR{{a[$4]=$5;next}}NR!=FNR{{c=$1; if(c in a){{print $0"\t"a[c]}}}}' {output.homerInput} {output.homerOutput} > {output.homerOutput2}
    """


rule mvp_plots:
  input:
    dir=join(dmr_dir, "bsseq"),
    bed=join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),
  output:
    miami=join(dmr_dir, "combP/{group}.miami.png"),
    violin=join(dmr_dir, "combP/{group}.violin.png"),
    #pie_cpg=join(dmr_dir, "combP/{group}.pie-cpg.png"),
    #pie_genic=join(dmr_dir, "combP/{group}.pie-genic.png"),
  params:
    rname="mvp_plots",
    groups='{group}_{chr}',
    dir=join(dmr_dir, "combP"),
    fdr_cutoff="0.05",
    highlight="20",
  resources:
    mem       = allocated("mem",       "mvp_plots", cluster),
    gres      = allocated("gres",      "mvp_plots", cluster),
    time      = allocated("time",      "mvp_plots", cluster),
    partition = allocated("partition", "mvp_plots", cluster),
  threads:
    int(allocated("threads", "mvp_plots", cluster)),
  shell:
    """
    mkdir -p {params.dir}
    module load R/4.1
    Rscript {params.script_dir}/mvp_plots.R {input.dir} {params.fdr_cutoff} {params.highlight} {output.miami} {output.violin}
    ## {output.pie_cpg} {output.pie_genic}
    """


rule combP_figs:
  input:
    expand(join(dmr_dir, "combP/{group}_{chr}.regions-p.bed.gz"),group=GROUPS,chr=CHRS),
  output:
    p1=join(dmr_dir, "combP/{group}.regions-p.bed"),
    f1=join(dmr_dir, "figs/{group}_GREATfig1.png"),
    f3=join(dmr_dir, "figs/{group}_manhatten.png"),
    o1=join(dmr_dir, "combP/{group}.GREATprocesses.txt"),
    o2=join(dmr_dir, "combP/{group}.GREATfunction.txt"),
  params:
    rname="combP_figs",
    dir=join(dmr_dir, "figs"),
    prefix=join(dmr_dir, "combP/{group}_chr"),
    script_dir=join(dmr_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "combP_figs", cluster),
    gres      = allocated("gres",      "combP_figs", cluster),
    time      = allocated("time",      "combP_figs", cluster),
    partition = allocated("partition", "combP_figs", cluster),
  threads:
    int(allocated("threads", "combP_figs", cluster)),
  shell:
    """
    module load R
    mkdir -p {params.dir}
    zcat {params.prefix}*.regions-p.bed.gz | grep -v "^#" > {output.p1}
    Rscript {params.script_dir}/combp_Manhatten.R {output.p1} {output.f3}
    Rscript {params.script_dir}/combp_GREAT.R {output.p1} {output.f1} {output.o1} {output.o2}
    """


rule manhatten:
  input:
    expand(join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),chr=CHRS,group=GROUPS),
  output:
    LIST=temp(join(dmr_dir,"figs/{group}.list")),
    MAN=join(dmr_dir,"figs/{group}_manhatten_beta.png"),
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


rule bamsort:
  input:
    bam=join(working_dir,"bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    bam=temp(join(working_dir,"bismarkAlign/bams/{samples}.bam")),
    bai=temp(join(working_dir,"bismarkAlign/bams/{samples}.bam.bai")),
  params:
    rname="pl:bamsort",
    reference="/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa",
  resources:
    mem       = allocated("mem",       "bamsort", cluster),
    gres      = allocated("gres",      "bamsort", cluster),
    time      = allocated("time",      "bamsort", cluster),
    partition = allocated("partition", "bamsort", cluster),
  threads:
    int(allocated("threads", "bamsort", cluster)),
  shell:
    """
    module load samtools
    samtools sort -@ 12 --reference /data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa {input.bam} > {output.bam}
    samtools index -@ 12 {output.bam}
    """


rule wgbstools:
  input:
    bam=join(working_dir,"bismarkAlign/bams/{samples}.bam"),
    bai=join(working_dir,"bismarkAlign/bams/{samples}.bam.bai"),
  output:
    pat=join(working_dir,"UXM/{samples}.pat.gz"),
  params:
    rname="pl:wgbstools",
    ref=species,
    outdir=join(working_dir,"UXM"),
    script_dir=join(working_dir,"scripts"),
  resources:
    mem       = allocated("mem",       "wgbstools", cluster),
    gres      = allocated("gres",      "wgbstools", cluster),
    time      = allocated("time",      "wgbstools", cluster),
    partition = allocated("partition", "wgbstools", cluster),
  threads:
    int(allocated("threads", "wgbstools", cluster)),
  shell:
    """
    mkdir -p {params.outdir}
    module load samtools bedtools bamtools
    /data/NHLBI_IDSS/references/UXM/wgbstools bam2pat --genome hg38  -L /data/NHLBI_IDSS/references/UXM/supplemental/Atlas.U25.l4.hg38.full.bed --out_dir {params.outdir} -@ 12 {input.bam}
    """


rule UXM:
  input:
    pat=join(working_dir,"UXM/{samples}.pat.gz"),
  output:
    pat=join(working_dir,"UXM/{samples}_deconv.250.csv"),
  params:
    rname="pl:UXM",
    atlas="/data/NHLBI_IDSS/references/UXM/supplemental/Atlas.U250.l4.hg38.full.tsv",
  resources:
    mem       = allocated("mem",       "UXM", cluster),
    gres      = allocated("gres",      "UXM", cluster),
    time      = allocated("time",      "UXM", cluster),
    partition = allocated("partition", "UXM", cluster),
  threads:
    int(allocated("threads", "UXM", cluster)),
  shell:
    """
    module load samtools bedtools bamtools
    /data/NHLBI_IDSS/references/UXM/uxm deconv {input.pat} -o {output.pat} --atlas {params.atlas} --ignore Colon-Fibro Dermal-Fibro Gallbladder Bone-Osteob
    """
