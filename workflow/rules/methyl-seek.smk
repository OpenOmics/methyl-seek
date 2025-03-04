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


# Standard data processing and QC related
# rules that always run regardless of the
# selected run mode (i.e. dmr or dcv)
rule raw_fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
      r1 = join(working_dir, "{samples}.R1.fastq.gz"),
      r2 = join(working_dir, "{samples}.R2.fastq.gz"),
    output:
      join(working_dir, "rawQC", "{samples}.R1_fastqc.zip"),
      join(working_dir, "rawQC", "{samples}.R1_fastqc.html"),
      join(working_dir, "rawQC", "{samples}.R2_fastqc.zip"),
      join(working_dir, "rawQC", "{samples}.R2_fastqc.html"),
    params:
      rname  = "raw_fastqc",
      outdir = join(working_dir,"rawQC"), 
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
      fastqc \\
        -o {params.outdir} \\
        -f fastq \\
        --threads {threads} \\
        --extract {input.r1} {input.r2}
      """


rule trimGalore:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome. Adapters are composed of synthetic
    sequences and should be removed prior to alignment. TrimGalore is a wrapper
    to cutadapt and fastqc. It will run both tools under the hood.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Trimmed FastQ file
    """
    input:
      r1 = join(working_dir, "{samples}.R1.fastq.gz"),
      r2 = join(working_dir, "{samples}.R2.fastq.gz"),
    output:
      r1 = join(working_dir, "trimGalore", "{samples}_val_1.fq.gz"),
      r2 = join(working_dir, "trimGalore", "{samples}_val_2.fq.gz")
    params:
      rname     = "trimgalore",
      outdir    = join(working_dir, "trimGalore"),
      tag       = "{samples}",
      fastqcdir = join(working_dir, "postTrimQC"),
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
      trim_galore \\
        --paired \\
        --cores {threads} \\
        --fastqc \\
        --clip_R1 10 \\
        --clip_R2 10 \\
        --three_prime_clip_R1 10 \\
        --three_prime_clip_R2 10 \\
        --length 50 \\
        --gzip \\
        --basename {params.tag} \\
        --output_dir {params.outdir} \\
        --fastqc_args "--outdir {params.fastqcdir}" \\
        {input.r1} \\
        {input.r2}
      """


rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
      file1 = join(working_dir, "trimGalore", "{samples}_val_1.fq.gz"),
      file2 = join(working_dir, "trimGalore", "{samples}_val_2.fq.gz"),
    output:
      out1 = join(working_dir,"FQscreen","{samples}_val_1_screen.txt"),
      out2 = join(working_dir,"FQscreen","{samples}_val_1_screen.png"),
      out3 = join(working_dir,"FQscreen","{samples}_val_2_screen.txt"),
      out4 = join(working_dir,"FQscreen","{samples}_val_2_screen.png"),
      out5 = join(working_dir,"FQscreen2","{samples}_val_1_screen.txt"),
      out6 = join(working_dir,"FQscreen2","{samples}_val_1_screen.png"),
      out7 = join(working_dir,"FQscreen2","{samples}_val_2_screen.txt"),
      out8 = join(working_dir,"FQscreen2","{samples}_val_2_screen.png")
    params:
      rname   ='fqscreen',
      outdir  = join(working_dir,"FQscreen"),
      outdir2 = join(working_dir,"FQscreen2"),
      fastq_screen_config  = "/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf",
      fastq_screen_config2 = "/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen_2.conf",
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
      # First pass of databases
      fastq_screen \\
        --conf {params.fastq_screen_config} \\
        --outdir {params.outdir} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input.file1} {input.file2}
      # Second pass of databases
      fastq_screen \\
        --conf {params.fastq_screen_config2} \\
        --outdir {params.outdir2} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input.file1} {input.file2}
      """


rule kraken_pe:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        fq1 = join(working_dir, "trimGalore", "{samples}_val_1.fq.gz"),
        fq2 = join(working_dir, "trimGalore", "{samples}_val_2.fq.gz"),
    output:
        krakenout  = join(working_dir, "kraken", "{samples}.trim.kraken_bacteria.out.txt"),
        krakentaxa = join(working_dir, "kraken", "{samples}.trim.kraken_bacteria.taxa.txt"),
        kronahtml  = join(working_dir, "kraken", "{samples}.trim.kraken_bacteria.krona.html"),
    params:
        rname  = "kraken",
        outdir = join(working_dir,"kraken"),
        # TODO: Add this as a shared references
        bacdb  = "/fdb/kraken/20210223_standard_kraken2",
        prefix = "{samples}",
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
      cd /lscratch/$SLURM_JOBID;
      # Copy kraken2 db to /lscratch or temp 
      # location to reduce filesystem strain
      cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
      kdb_base=$(basename {params.bacdb})
      kraken2 --db /lscratch/$SLURM_JOBID/${{kdb_base}} \\
        --threads {threads} --report {output.krakentaxa} \\
        --output {output.krakenout} \\
        --gzip-compressed \\
        --paired {input.fq1} {input.fq2}
      # Generate Krona Report
      cut -f2,3 {output.krakenout} | \\
        ktImportTaxonomy - -o {output.kronahtml}
      """


# TODO: This rule should run conditionally,
# it should only run if the user provided a
# their own genomic fasta file. Need to up-
# date the cli to allow for this. If user
# provided a genomic fasta file, then this
# rule should run to build the required
# reference files for Bismark.
rule prep_bisulphite_genome:
    """
    Reference-building step to prepare the reference genome of interest for bisulfite
    alignments. You need to specify a directory containing the genome you want to align
    your reads against (please be aware that the bismark_genome_preparation script
    expects FASTA files in this folder (with either .fa or .fasta extension, single
    or multiple sequence entries per file). Bismark will create two individual folders
    within this directory, one for a C->T converted genome and the other one for the
    G->A converted genome.
    @Input:
        User provided genomic FASTA file (singleton)
    @Output:
        C to T converted genomic FASTA file,
        G to A converted genomic FASTA file
    """
    input:
      bisulphite_fa
    output:
      join(working_dir, "refs", "Bisulfite_Genome", "CT_conversion", "genome_mfa.CT_conversion.fa"),
      join(working_dir, "refs", "Bisulfite_Genome", "GA_conversion", "genome_mfa.GA_conversion.fa"),
    params:
      rname   = "prep_bisulphite_genome",
      outdir  = join(working_dir, "refs"),
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
      cd {params.outdir}
      # Build custom reference genome for Bismark 
      bismark_genome_preparation \\
        --verbose \\
        --parallel {threads} \\
        {params.dir}
      """


rule bismark_align:
    """
    Data-processing step perform bisulfite alignment and methylation calling.
    Bismark requires the user to specify only two things: 
      1. The directory containing the genome of interest. This folder must contain
         the unmodified genome (as .fa or .fasta files) as well as the two bisulfite
         genome subdirectories which were generated in the Bismark Genome prep 
         step (see above).
      2. The FastQ file(s) to be analysed.
    @Input:
        Trimmed FastQ files (scatter),
        Output of bismark_genome_preparation step
    @Output:
        BAM file containing reads aligned to the prepared reference genome,
        Samtools flagstat output file
    """
    input:
      F1 = join(working_dir, "trimGalore", "{samples}_val_1.fq.gz"),
      F2 = join(working_dir, "trimGalore", "{samples}_val_2.fq.gz"),
      bismark_ref = bisulphite_genome_path
    output:
      B1  = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.bam"),
      B2  = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.flagstat"),
      FQ1 = temp(join(working_dir, "bismarkAlign", "{samples}_val_1.fq.gz_unmapped_reads_1.fq.gz")),
      FQ2 = temp(join(working_dir, "bismarkAlign", "{samples}_val_2.fq.gz_unmapped_reads_2.fq.gz")),
      FQ3 = temp(join(working_dir, "bismarkAlign", "{samples}_val_1.fq.gz_ambiguous_reads_1.fq.gz")),
      FQ4 = temp(join(working_dir, "bismarkAlign", "{samples}_val_2.fq.gz_ambiguous_reads_2.fq.gz")),
    params:
      rname   = "bismark_align",
      outdir  = join(working_dir, "bismarkAlign"),
      outbam  = join(working_dir, "bismarkAlign", "{samples}_val_1_bismark_bt2_pe.bam"),
      R1      = join(working_dir, "bismarkAlign", "{samples}_val_1_bismark_bt2_PE_report.txt"),
      R2      = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_PE_report.txt"),
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
      bismark \\
        --multicore {threads} \\
        --temp_dir /lscratch/$SLURM_JOBID/ \\
        --bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 4 --score_min L,-0.6,-0.6 \\
        --output_dir {params.outdir} \\
        --genome {input.bismark_ref} \\
        -1 {input.F1} \\
        -2 {input.F2}
      mv {params.outbam} {output.B1}
      mv {params.R1} {params.R2}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """


rule bismark_dedup:
    """
    Data-processing step to remove alignments to the same position in the genome
    from the Bismark mapping output. Sequences which align to the same genomic
    position but on different strands are scored individually. It is important to
    note that deduplication is not recommended for RRBS, amplicon or other target
    enrichment-type libraries!
    @Input:
        BAM file containing reads aligned to the prepared reference genome (scatter),
    @Output:
        Dedup CRAM file
    """
    input:
      F1 = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.bam"),
    output:
      T1 = temp(join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.temp.bam")),
      B1 = temp(join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.bam")),
      B2 = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.flagstat"),
      C1 = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.cram"),
    params:
      rname  = "bismark_dedup",
      outdir = join(working_dir, "bismarkAlign"),
      prefix = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe"),
      genome = bisulphite_fa,
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
      cd {params.outdir}
      samtools view \\
        -hb {input.F1} \\
      | samtools sort \\
        -n -@ {threads} \\
        -O BAM \\
        -o {output.T1}
      deduplicate_bismark \\
        --paired \\
        --bam \\
        --outfile {params.prefix} \\
        {output.T1}
      samtools view \\
        -T {params.genome} \\
        -@ 16 \\
        -h \\
        -C {output.B1} \\
      > {output.C1}
      samtools flagstat \\
        -@ {threads} \\
        {output.B1} \\
      > {output.B2}
      """


rule bismark_extract:
  """
  Data-processing step extracts methylation calls for every single C analysed. The position 
  of every single C will be written out to a new output file, depending on its context
  (CpG, CHG or CHH), whereby methylated Cs will be labelled as forward reads (+), 
  non-methylated Cs as reverse reads (-), etc.
  @Input:
      BAM file containing reads aligned to the prepared reference genome (scatter),
  @Output:
      Text and BedGraph file containing CpG calls
  """
  input:
    bam = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    cov = join(working_dir, "CpG", "{samples}", "{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
    bed = join(working_dir, "CpG", "{samples}", "{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
  params:
    rname  = 'bismark_extract',
    outdir = join(working_dir, "CpG", "{samples}"),
    genome = bisulphite_fa,
    bismark_index = join(bisulphite_genome_path,species),
  resources:
    # Scaling memory and partion allocation
    # based on the number of retries
    mem       = lambda wc, attempt: allocated("mem", "bismark_extract", cluster) if attempt == 1 else "900g",
    partition = lambda wc, attempt: allocated("partition", "bismark_extract", cluster) if attempt == 1 else "largemem",
    gres      = allocated("gres", "bismark_extract", cluster),
    time      = allocated("time", "bismark_extract", cluster),
  threads:
    int(allocated("threads", "bismark_extract", cluster)),
  shell:
    """
    module load bismark samtools bowtie
    bismark_methylation_extractor \\
      --paired-end \\
      --no_overlap \\
      --multicore {threads} \\
      --gzip \\
      --report \\
      --bedGraph \\
      --counts \\
      --buffer_size '50%' \\
      --no_header \\
      --cytosine_report \\
      --output {params.outdir} \\
      --scaffolds \\
      --genome_folder {params.bismark_index} \\
      {input.bam}
    """


rule multiqc:
  """
  Reporting step to aggregate sample summary statistics and quality-control
  information across all samples. This will be one of the last steps of the 
  pipeline. The inputs listed here are to ensure that this step runs last. 
  During runtime, MultiQC will recurively crawl through the working directory
  and parse files that it supports.
  @Input:
      Lists of files to ensure this step runs last (gather)
  @Output:
      Interactive MulitQC report
  """
  input:
    expand(join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
    expand(join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.bam"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore", "{samples}_val_1.fq.gz"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore", "{samples}_val_2.fq.gz"),samples=SAMPLES),
  output:
    join(working_dir, "bismark_summary_report.txt"),
    join(working_dir, "bismark_summary_report.html"),
    join(working_dir, "multiqc_report.html"),
  params:
    rname="multiqc",
    outdir=working_dir,
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
    mv bismark_summary_report.html {params.outdir}/
    mv bismark_summary_report.txt {params.outdir}/
    cd {params.outdir}
    multiqc --ignore '*/.singularity/*' -f --interactive .
    """
    