#####################################################################################################
# Bisulphite sequencing (methyl-seek) analysis workflow
#
# Execution mode: run and dmr
# Last Modified: August, 2022
#
#####################################################################################################

from os.path import join
from snakemake.io import expand, glob_wildcards
from snakemake.utils import R
from os import listdir
import pandas as pd

# Global workflow variables
sample_file= config["samples"]
rawdata_dir= config["rawdata_dir"]
working_dir= config["result_dir"]

# References
hg38_fa= config["hg38_fa"]
phage_fa= config["phage_fa"]
hg38_gtf= config["hg38_gtf"]
hg38_rRNA_intervals= config["hg38_rRNA_intervals"]
hg38_bed_ref= config["hg38_bed_ref"]
hg38_refFlat= config["hg38_refFlat"]
bisulphite_genome_path= config["bisulphite_ref"]
phage_genome_path= config["phage_ref"]
bisulphite_fa= config["bisulphite_fa"]
species= config["species"]
REF_ATLAS=config["REF_ATLAS"]
CpG_MAP_TABLE=config["CpG_MAP_TABLE"]
mode=config["Mode"]

# sample list
df = pd.read_csv(sample_file, header=0, sep='\t')
SAMPLES=list(set(df['samples'].tolist()))
GROUPS=list(set(df['group'].tolist()))


CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX']

RN = ['R1', 'R2']


dmr_dir = working_dir

if mode == "dmr":
  # contrasts list
  df2 = pd.read_csv(working_dir + "/contrasts.txt", header=0, sep='\t')
  GROUPS=list(set(df2['comparisons'].tolist()))
  GRPSAMPLES=list(set(df2['samples'].tolist()))
  dmr_dir = join(working_dir, "dmr", GROUPS[0])
  #mkdir -p dmr_dir



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
        outputs.append(expand(join(dmr_dir, "bsseq/{group}_{chr}_betas_wtFDR_unfiltered.bed"),chr=CHRS,group=GROUPS))
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
        outputs.append(expand(join(dmr_dir, "combP/{group}.pie-cpg.png"),group=GROUPS))
        outputs.append(expand(join(dmr_dir, "combP/{group}.pie-genic.png"),group=GROUPS))

    if mode == "dcv":
        outputs.append(join(working_dir, "multiqc_report.html"))
        outputs.append(expand(join(working_dir, "CpG_CSV/{samples}.csv"),samples=SAMPLES))
        outputs.append(expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES))
        outputs.append(expand(join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),samples=SAMPLES))
        outputs.append(join(working_dir, "deconvolution_CSV/total.csv"))
        outputs.append(join(working_dir, "deconvolution_CSV/total_deconv_output.csv"))
        outputs.append(join(working_dir, "deconvolution_CSV/total_deconv_plot.png"))

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
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.out.txt"),samples=SAMPLES),
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.taxa.txt"),samples=SAMPLES),
      expand(join(working_dir, "kraken","{samples}.trim.kraken_bacteria.krona.html"),samples=SAMPLES),

      #FQscreen output
      expand(join(working_dir, "FQscreen","{samples}_val_1_screen.txt"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen","{samples}_val_1_screen.png"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen","{samples}_val_2_screen.txt"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen","{samples}_val_2_screen.png"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen2","{samples}_val_1_screen.txt"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen2","{samples}_val_1_screen.png"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen2","{samples}_val_2_screen.txt"),samples=SAMPLES),
      expand(join(working_dir, "FQscreen2","{samples}_val_2_screen.png"), samples=SAMPLES),

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
      expand(join(working_dir, "rseqc/{samples}.inner_distance_freq.txt"),samples=SAMPLES),
      expand(join(working_dir, "rseqc/{samples}.strand.info"),samples=SAMPLES),
      expand(join(working_dir, "rseqc/{samples}.Rdist.info"),samples=SAMPLES),
      expand(join(working_dir, "trimGalore/{samples}_insert_sizes.txt"),samples=SAMPLES),
      expand(join(working_dir, "CpG/{samples}/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),samples=SAMPLES),


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
    threads:
      16
    shell:
      """
      module load trimgalore/0.6.7
      module load fastqc/0.11.9
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
    threads: 4
    shell: """
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
        dir=join(working_dir,"kraken"),
        #bacdb="/fdb/kraken/20170202_bacteria",
        bacdb="/fdb/kraken/20210223_standard_kraken2",
        prefix="{samples}",
    threads: 24
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
    threads:
      16
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
    threads:
      16
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
  params:
    rname='pl:bismark_extract',
    bismark_index=join(bisulphite_genome_path,species),
    outdir=join(working_dir, "CpG/{samples}"),
    genome=bisulphite_fa,
  shell:
    """
    mkdir -p {params.outdir}
    module load bismark samtools bowtie
    bismark_methylation_extractor --paired-end --no_overlap --multicore 8 --gzip --report --bedGraph --counts --buffer_size 100G --no_header --cytosine_report --output {params.outdir} --scaffolds --genome_folder {params.bismark_index} {input.bam}
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
    refflat=hg38_refFlat,
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

############### Deconvolution rules begin here
## run following rules in dcv mode only
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
  params:
    rname="merge_tables",
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
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve.py --atlas_path {params.ref} --plot --residuals {input}
    """

############### Differential methylation rules begin here
## run following rules in dmr mode only
#################

rule bsseq_bismark:
  input:
    bizfile=join(working_dir,"contrasts.txt"),
  output:
    bed1=join(dmr_dir, "bsseq/{group}_{chr}_betas_pval.bed"),
    bed2=join(dmr_dir, "bsseq/{group}_{chr}_betas_wtFDR_unfiltered.bed"),
  params:
    rname="bsseq_bismark",
    chr='{chr}',
    dir=directory(join(dmr_dir, "bsseq")),
    cov="2",
    sample_prop="0.25",
    script_dir=join(working_dir, "dmr","scripts"),
  threads:
    4
  shell:
    """
    module load R
    mkdir -p {params.dir}
    Rscript {params.script_dir}/bsseq_lm.R {params.chr} {input.bizfile} {output.bed1} {params.sample_prop} {params.cov} {output.bed2}
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
    bed=join(dmr_dir, "bsseq/{group}_{chr}_betas_wtFDR_unfiltered.bed"),
  output:
    miami=join(dmr_dir, "combP/{group}.miami.png"),
    violin=join(dmr_dir, "combP/{group}.violin.png"),
    pie_cpg=join(dmr_dir, "combP/{group}.pie-cpg.png"),
    pie_genic=join(dmr_dir, "combP/{group}.pie-genic.png"),
  params:
    rname="mvp_plots",
    groups='{group}_{chr}',
    dir=join(dmr_dir, "combP"),
    fdr_cutoff="0.05",
    highlight="20"
  shell:
    """
    mkdir -p {params.dir}
    module load R/4.1
    Rscript {params.script_dir}/mvp_plots.R {input.dir} {params.fdr_cutoff} {params.highlight} {output.miami} {output.violin} {output.pie_cpg} {output.pie_genic}
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
  shell:
    """
    mkdir -p {params.dir}
    ls {params.dir}/{wildcards.group}_*_betas_pval.bed > {output.LIST}
    module load R
    Rscript {params.script_dir}/bsseq_Manhatten.R {output.LIST} {output.MAN}
    """
