#####################################################################################################
# Bisulphite sequencing (methylseq) analysis workflow
#
# This Snakefile is designed to perform QC of raw bisulphite sequencing data,
# align reads to human (hg38) bisulfite genome using either Bismark (bowtie2) or
# bwa_meth, summarizes alignment stats, extracts methylated loci from the genome.
#
# Created: December 21, 2021
# Contact details: Tom Hill (tom.hill@nih.gov), Neelam Redekar (neelam.redekar@nih.gov),
# Skyler Kuhn (skyler.kuhn@nih.gov), Asya Khleborodova (asya.khleborodova@nih.gov)
#
# Last Modified: March 8, 2022
#
#####################################################################################################

from os.path import join
from snakemake.io import expand, glob_wildcards
from snakemake.utils import R
from os import listdir
import pandas as pd

##
## Locations of working directories and reference genomes for analysis
##
sample_file= config["samples"]
rawdata_dir= config["rawdata_dir"]
working_dir= config["result_dir"]
hg38_fa= config["hg38_fa"]
phage_fa= config["phage_fa"]
hg38_gtf= config["hg39_gtf"] ################## New edition
hg38_rRNA_intervals= config["rRNA_intervals"]  ################## New edition
hg38_bed_ref= config["bed_ref"]  ################## New edition
hg39_refFlat= config["hg39_refFlat"] ################## New edition
bisulphite_genome_path= config["bisulphite_ref"]
phage_genome_path= config["phage_ref"]
bisulphite_fa= config["bisulphite_fa"]
species= config["species"]

REF_ATLAS=config["REF_ATLAS"]
CpG_MAP_TABLE=config["CpG_MAP_TABLE"]

##
## Read in the masterkey file for 3 tab-delimited columns of samples, groups and comparison
## Each sample can be in the file multiple times if used in multiple comparisons, but will only be mapped/process once.
##
## e.g.
##
##sample	group	comp
##S1	GA	GAvsGB
##S2	GA	GAvsGB
##S3	GB	GAvsGB
##S4	GB	GAvsGB
##S5	GC	GAvsGC
##S6	GC	GAvsGC
##S1	GA	GAvsGC
##S2	GA	GAvsGC

## The file requires these headings as they are used in multiple rules later on.

## Here we read in the samples file generated and begin processing the data, printing out the samples and group comparisons

df = pd.read_csv(sample_file, header=0, sep='\t')

SAMPLES=list(set(df['sample'].tolist()))
GROUPS=list(set(df['comp'].tolist()))

print(SAMPLES)
print(len(SAMPLES))
print(GROUPS)
print(len(GROUPS))

CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX']

RN = ['R1', 'R2']

rule All:
    input:
      # Creating data links:
      expand(join(working_dir, "raw/{samples}.{rn}.fastq.gz"), samples=SAMPLES, rn=RN),

      # Checking data quality:
      expand(join(working_dir, "rawQC/{samples}.{rn}_fastqc.html"), samples=SAMPLES, rn=RN),
      expand(join(working_dir, "trimQC/{samples}.{rn}.pe_fastqc.html"), samples=SAMPLES, rn=RN),

      # Quality trimming output:
      expand(join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),samples=SAMPLES),
      expand(join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),samples=SAMPLES),

      # kraken output
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.out.txt"),samples=SAMPLES),
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.taxa.txt"),samples=SAMPLES),
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.krona.html"),samples=SAMPLES),

      #FQscreen output
      expand(join(working_dir,"FQscreen","{samples}.R1.trim_screen.txt"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen","{samples}.R1.trim_screen.png"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen","{samples}.R2.trim_screen.txt"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen","{samples}.R2.trim_screen.png"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen2","{samples}.R1.trim_screen.txt"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen2","{samples}.R1.trim_screen.png"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen2","{samples}.R2.trim_screen.txt"),samples=SAMPLES),
      expand(join(working_dir,"FQscreen2","{samples}.R2.trim_screen.png"), samples=SAMPLES),

      # bisulphite genome preparation
      join(bisulphite_genome_path, species, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(bisulphite_genome_path, species, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),

      # bismark align to human reference genomes
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.flagstat"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),samples=SAMPLES),

      # get alignment statistics
      expand(join(working_dir,"bismarkAlign/{samples}.RnaSeqMetrics.txt"),samples=SAMPLES),
      expand(join(working_dir,"bismarkAlign/{samples}.flagstat.concord.txt"),samples=SAMPLES),
      expand(join(working_dir,"rseqc/{samples}.inner_distance_freq.txt"),samples=SAMPLES),
      expand(join(working_dir,"rseqc/{samples}.strand.info"),samples=SAMPLES),
      expand(join(working_dir,"rseqc/{samples}.Rdist.info")samples=SAMPLES),
      expand(join(working_dir,"QualiMap/{samples}/qualimapReport.html"),samples=SAMPLES),
      expand(join(working_dir,"QualiMap/{samples}/genome_results.txt"),samples=SAMPLES),
      expand(join(working_dir, "preseq/{sample}.ccurve"),samples=SAMPLES),
      expand(join(working_dir, "trimGalore/{samples}_insert_sizes.txt"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bai"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.star.duplic"), samples=SAMPLES),

      # extract CpG profile with methyldackel
      expand(join(working_dir, "CpG/{samples}.bedGraph"),samples=SAMPLES),

      # generate multiqc output
      "multiqc_report.html",

      # Deconvolution output
      expand(join(working_dir, "CpG_CSV/{samples}.csv"),samples=SAMPLES),
      expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
      expand(join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),samples=SAMPLES),
      join(working_dir, "deconvolution_CSV/total.csv"),
      join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
      join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),


## Copy raw data to working directory
rule raw_data_links:
    input:
      join(rawdata_dir, "{samples}.{rn}.fastq.gz")
    output:
      join(working_dir, "raw/{samples}.{rn}.fastq.gz")
    params:
      rname="raw_data_links",
      dir=directory(join(working_dir, "raw")),
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
    threads:
      2
    shell:
      """
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
      """

## Trim raw data
rule trimmomatic:
    input:
      F1=join(working_dir, "raw/{samples}.R1.fastq.gz"),
      F2=join(working_dir, "raw/{samples}.R2.fastq.gz"),
    output:
      PE1=temp(join(working_dir, "trimmed_reads/{samples}.R1.pe.fastq.gz")),
      UPE1=temp(join(working_dir, "trimmed_reads/{samples}.R1.ue.fastq.gz")),
      PE2=temp(join(working_dir, "trimmed_reads/{samples}.R2.pe.fastq.gz")),
      UPE2=temp(join(working_dir, "trimmed_reads/{samples}.R2.ue.fastq.gz"))
    params:
      rname="trimmomatic",
      dir=directory(join(working_dir, "trimmed_reads")),
      batch='--cpus-per-task=8 --partition=norm --gres=lscratch:180 --mem=25g --time=20:00:00',
      command='ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50'
    threads:
      8
    shell:
      """
      module load trimmomatic/0.39
      mkdir -p {params.dir}
      java -jar $TRIMMOJAR PE -phred33 -threads {threads} {input.F1} {input.F2} {output.PE1} {output.UPE1} {output.PE2} {output.UPE2} {params.command}
      """

## Second round of trimming/filtering
rule trimGalore:
    input:
      F1=join(working_dir, "trimmed_reads/{samples}.R1.pe.fastq.gz"),
      F2=join(working_dir, "trimmed_reads/{samples}.R2.pe.fastq.gz"),
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
    threads:
      16
    shell:
      """
      module load trimgalore/0.6.7
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      trim_galore --paired --cores {threads} {params.command} --basename {params.tag} --output_dir {params.dir} --fastqc_args "--outdir {params.fastqcdir}"  {input.F1} {input.F2}
      """

## Run fastqc on filtered/trimmed data to visually assess quality for R1
rule trim_fastqc:
    input:
      join(working_dir, "trimmed_reads/{samples}.{rn}.pe.fastq.gz"),
    output:
      join(working_dir, "trimQC/{samples}.{rn}.pe_fastqc.html"),
    params:
      rname="trim_fastqc",
      dir=directory(join(working_dir, "trimQC")),
      batch='--cpus-per-task=2 --mem=8g --time=8:00:00',
    threads:
      2
    shell:
      """
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
      """

################## New edition - started
rule bbmerge:
    input:
      R1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      R2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      join(working_dir, "trimGalore/{samples}_insert_sizes.txt"),
    params:
      rname='pl:bbmerge',
    threads: 4
    shell: """
      # Get encoding of Phred Quality Scores
      module load python
      encoding=$(python phred_encoding.py {input.R1})
      echo "Detected Phred+${{encoding}} ASCII encoding"

      module load bbtools/38.87
      bbtools bbmerge-auto in1={input.R1} in2={input.R2} qin=${{encoding}} \
      ihist={output} k=62 extend2=200 rem ecct -Xmx64G
          """

rule fastq_screen:
    input:
      file1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      file2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      out1=join(working_dir,"FQscreen","{samples}.R1.trim_screen.txt"),
      out2=join(working_dir,"FQscreen","{samples}.R1.trim_screen.png"),
      out3=join(working_dir,"FQscreen","{samples}.R2.trim_screen.txt"),
      out4=join(working_dir,"FQscreen","{samples}.R2.trim_screen.png"),
      out5=join(working_dir,"FQscreen2","{samples}.R1.trim_screen.txt"),
      out6=join(working_dir,"FQscreen2","{samples}.R1.trim_screen.png"),
      out7=join(working_dir,"FQscreen2","{samples}.R2.trim_screen.txt"),
      out8=join(working_dir,"FQscreen2","{samples}.R2.trim_screen.png")
    params:
      rname='pl:fqscreen',
      outdir = join(working_dir,"FQscreen"),
      outdir2 = join(working_dir,"FQscreen2"),
      fastq_screen_config="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf",
      fastq_screen_config2="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen_2.conf",
    threads: 24
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
        outdir=join(working_dir,"kraken"),
        bacdb="/fdb/kraken/20170202_bacteria"
    threads: 24
    shell:
      """
      module load kraken/1.1
      module load kronatools/2.7
      if [ ! -d {params.dir} ];then mkdir {params.dir};fi

      cd /lscratch/$SLURM_JOBID;
      cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;

      kdb_base=$(basename {params.bacdb})
      kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload--paired {input.F1} {input.F2}
      kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
      cut -f 2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
      """

################## New edition - ended


## prepare bi-sulphite genome
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
    threads:
      16
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
    params:
      rname="bismark_align",
      dir=directory(join(working_dir, "bismarkAlign")),
      genome_dir=directory(join(bisulphite_genome_path, species)),
      command="--bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 4 --score_min L,-0.6,-0.6",
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=100g --time=10:00:00',
      outbam=join(working_dir, "bismarkAlign/{samples}_val_1_bismark_bt2_pe.bam"),
    threads:
      16
    shell:
      """
      module load bismark/0.23.0 samtools
      mkdir -p {params.dir}
      bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
      mv {params.outbam} {output}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """

rule bismark_dedup:
    input:
      F1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),
    output:
      T1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.deduplicated.bam")),
      B1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam")),
      B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),
    params:
      rname="bismark_dedup",
      dir=directory(join(working_dir, "bismarkAlign")),
    threads:
      16
    shell:
      """
      module load bismark/0.23.0
      module load samtools/1.15
      cd {params.dir}
      deduplicate_bismark --paired --bam --outfile {output.B1} {input.F1}
      samtools view -hb {output.T1} | samtools sort -@ {threads} -O BAM -o {output.B1}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
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
    threads:
      16
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
  threads:
    16
  shell:
    """
    module load bismark/0.23.0
    mkdir -p {params.dir}
    bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
    """

################## New edition - start
rule picard:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),
    bai=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bai"),
    metrics=join(working_dir, "bismarkAlign/{samples}.star.duplic")
  params:
    rname='pl:picard',
    sampleName="{sample}",
  threads: 6
  shell:
    """
    module load samtools/1.15
    module load picard/2.26.9

    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups \
    I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.bam \
    TMP_DIR=/lscratch/$SLURM_JOBID RGID=id VALIDATION_STRINGENCY=SILENT RGLB=library RGPL=illumina RGPU=machine RGSM=sample;

    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates \
    I=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.bam \
    O=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bam \
    TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.metrics};

    mv /lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bam {output.bam};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bai {output.bai};
    sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.metrics};

    """

rule preseq:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),
  output:
    ccurve = join(working_dir, "preseq/{sample}.ccurve"),
  params:
    rname = "pl:preseq",
    dir = directory(join(working_dir, "preseq")),
  shell:
    """
    module load preseq/3.1.2
    mkdir -p {params.dir}
    preseq c_curve -B -o {output.ccurve} {input.bam}
    """

rule qualibam:
  input:
    bamfile=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),
  output:
    report=join(working_dir,"QualiMap/{samples}/qualimapReport.html"),
    results=join(working_dir,"QualiMap/{samples}/genome_results.txt"),
  params:
    rname='pl:qualibam',
    outdir=join(working_dir,"QualiMap/{samples}"),
    gtfFile=hg38_gtf,
  threads: 16
  shell:
    """
    module load qualimap/2.2.1
    mkdir -p {params.outdir}
    unset DISPLAY;
    qualimap bamqc -bam {input.bamfile} --feature-file {params.gtfFile} -outdir {params.outdir} -nt {threads} --java-mem-size=120G
    """

rule rseqc:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    out1=join(working_dir,"rseqc/{samples}.strand.info"),
    out2=join(working_dir,"rseqc/{samples}.Rdist.info")
  params:
    bedref=hg38_bed_ref
    dir=join(working_dir,"rseqc")
    rname="pl:rseqc",
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
    refflat=hg39_refFlat,
    rrnalist=hg38_rRNA_intervals,
    picardstrand="SECOND_READ_TRANSCRIPTION_STRAND",
    statscript=join("workflow", "scripts", "bam_count_concord_stats.py"),
  shell:
    """
    module load python/3.8 samtools/1.15 picard/2.26.9
    java -Xmx110g -jar ${{PICARDJARPATH}}/picard.jar CollectRnaSeqMetrics REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist} STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    samtools flagstat {input.file1} > {output.outstar2};
    ## python3 {params.statscript} {input.file1} >> {output.outstar2} ## does require "Rstat" not available in biowulf
    """

################## New edition - ended

rule extract_CpG_bismark:
    input:
      F1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
    output:
      B1=join(working_dir, "CpG/{samples}.bedGraph"),
    params:
      rname="extract_CpG",
      dir=directory(join(working_dir, "CpG")),
      genome=hg38_fa,
      prefix=join(working_dir,"CpG/{samples}"),
    threads:
      16
    shell:
      """
      module load python
      module load samtools
      mkdir -p {params.dir}
      source /data/$USER/conda/etc/profile.d/conda.sh
      conda activate meth
      module load samtools/1.9
      MethylDackel mbias -@ {threads} {params.genome} {input.F1} {params.prefix}
      MethylDackel extract -o {params.prefix} -@ {threads} {params.genome} {input.F1}
      """

rule cleanup_bams:
  input:
    B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
    G1=join(working_dir, "CpG/{samples}.bedGraph"),
  output:
    C2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),
  params:
    rname="cleanup_bams",
    genome=hg38_fa,
    FQ1=join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_unmapped_reads_1.fq.gz"),
    FQ2=join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_unmapped_reads_2.fq.gz"),
    FQ3=join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_ambiguous_reads_1.fq.gz"),
    FQ4=join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_ambiguous_reads_2.fq.gz"),
  threads:
    8
  shell:
    """
      module load samtools
      samtools -h -C -@ {threads} -T {params.genome} {input.B2} > {output.C2}
      rm {params.FQ1}
      rm {params.FQ2}
      rm {params.FQ3}
      rm {params.FQ4}
    """

rule multiqc:
  input:
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),samples=SAMPLES),
  output:
    "multiqc_report.html",
  params:
    dir=working_dir,
    bis_dir=directory(join(working_dir,"bismarkAlign")),
    script_dir=join(working_dir,"scripts"),
  shell:
    """
    module load multiqc/1.9 bismark
    cd {params.bis_dir}
    bismark2report
    bismark2summary
    cd {params.dir}
    multiqc --ignore '*/.singularity/*' -f --interactive .
    """

############### Deconvolution rules begin here
rule get_CpG:
	input:
		join(working_dir, "CpG/{samples}.bedGraph"),
	output:
		join(working_dir, "CpG_CSV/{samples}.csv"),
	params:
		rname="get_CpG",
		cutoff=5,
		script_dir=join(working_dir,"scripts"),
		dir1=join(working_dir,"CpG_CSV"),
		dir2=join(working_dir,"deconvolution_CSV"),
	shell:
		"""
		mkdir -p {params.dir1}
		mkdir -p {params.dir2}
		module load R
		Rscript {params.script_dir}/get_methy.R {input} {wildcards.samples} {params.cutoff} {output}
		"""

rule get_overlap_meth:
  input:
    join(working_dir, "deconvolution_CSV/{samples}.meth.csv"),
  output:
    join(working_dir, "deconvolution_CSV/{samples}.csv"),
  params:
    rname="get_overlap_meth",
    map_table=CpG_MAP_TABLE,
  shell:
    """
    df_ref=pd.read_csv(params.map_table,sep='\t',header=None)
    df_ref.columns=['chromosome','start','end','cgid']
    df_ref=df_ref.loc[(df_ref['chromosome'].isin(CHRS)),]
    dfm=pd.read_csv(input[0])
    dfm=pd.merge(df_ref,dfm,on=['chromosome','start','end'],how='inner')
    dfm=dfm.drop(labels=['chromosome','start','end'],axis=1)
    dfm=dfm.set_index('cgid')
    dfm.to_csv(output[0])
    """

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
  shell:
"""
  module load python
  cd {params.dir}
  python {params.script_dir}/deconvolve.py --atlas_path {params.ref} --plot --residuals {input}  > {output}  2>&1
"""

rule merge_tables:
  input:
    expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
  output:
    join(working_dir, "deconvolution_CSV/total.csv"),
  shell:
    """
    dfm=pd.read_csv(input[0])
    for f in input[1:]:
    df=pd.read_csv(f)
    dfm=pd.merge(dfm,df,on='cgid',how='outer')
    dfm.to_csv(output[0],index=False)
    """

rule run_deconv_merged:
  input:
    join(working_dir, "deconvolution_CSV/total.csv"),
  output:
    join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
    join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),
  params:
    ref=REF_ATLAS,
    dir=join(working_dir,"deconvolution_CSV"),
    script_dir=join(working_dir,"scripts"),
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve.py --atlas_path {params.ref} --plot --residuals {input}
    """
