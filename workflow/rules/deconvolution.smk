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


# Deconvolution related rules, this rules
# should only run if the pipeline is run
# with the hg38 reference genome
rule get_CpG:
    input:
        join(working_dir, "CpG", "{samples}", "{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
    output:
        join(working_dir, "CpG_CSV", "{samples}.csv"),
    params:
        rname      = "get_CpG",
        cutoff     = 5,
        script_dir = join(working_dir, "workflow", "scripts"),
        dir1       = join(working_dir, "CpG_CSV"),
        dir2       = join(working_dir, "deconvolution_CSV"),
    resources:
        mem       = allocated("mem",       "get_CpG", cluster),
        gres      = allocated("gres",      "get_CpG", cluster),
        time      = allocated("time",      "get_CpG", cluster),
        partition = allocated("partition", "get_CpG", cluster),
    threads:
        int(allocated("threads", "get_CpG", cluster)),
    shell:
        """
        module load R
        mkdir -p {params.dir1}
        mkdir -p {params.dir2}
        Rscript {params.script_dir}/get_methy.R \\
            {input} \\
            {wildcards.samples} \\
            {params.cutoff} \\
            {output}
        """


rule get_overlap_meth:
    input:
        join(working_dir, "CpG_CSV", "{samples}.csv"),
    output:
        join(working_dir, "deconvolution_CSV", "{samples}.csv"),
    params:
        rname     = "get_overlap_meth",
        map_table = CpG_MAP_TABLE,
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
        join(working_dir, "deconvolution_CSV", "{samples}.csv"),
    output:
        join(working_dir, "deconvolution_CSV", "{samples}_deconv.log"),
    params:
        script_dir = join(working_dir, "workflow", "scripts"),
        outdir     = join(working_dir, "deconvolution_CSV"),
        rname      = "run_deconv",
        ref        = REF_ATLAS,
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
        cd {params.outdir}
        python {params.script_dir}/deconvolve2.py \\
            --atlas_path {params.ref} \\
            --residuals {input} \\
        > {output}  2>&1
        """


rule merge_tables:
    input:
        expand(join(working_dir, "deconvolution_CSV", "{samples}.csv"), samples=SAMPLES),
    output:
        join(working_dir, "deconvolution_CSV", "total.csv"),
    params:
        rname = "merge_tables",
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
        join(working_dir, "deconvolution_CSV", "total.csv"),
    output:
        join(working_dir, "deconvolution_CSV", "total_deconv_output.csv"),
    params:
        rname      = "run_deconv_merged",
        ref        = REF_ATLAS,
        outdir     = join(working_dir, "deconvolution_CSV"),
        script_dir = join(working_dir, "workflow", "scripts"),
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
        cd {params.outdir}
        python {params.script_dir}/deconvolve2.py \\
            --atlas_path {params.ref} \\
            --residuals {input}
        """


# cfDNAme rules 
rule sorting_CpG_bedgraph:
    input:
        bed = join(working_dir, "CpG", "{samples}", "{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
    output:
        bed = temp(join(working_dir, "CpG", "{samples}", "{samples}.mapped_autosomal_CpG.bedGraph.tmp")),
    params:
        rname = "format1",
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
        zcat {input.bed} \\
            | tail -n +2 \\
            | bedtools sort -i - \\
        > {output.bed}
        """


rule liftover_bedgraph:
    input:
        bed = join(working_dir,"CpG", "{samples}", "{samples}.mapped_autosomal_CpG.bedGraph.tmp"),
    output:
        graph = temp(join(working_dir,"CpG", "{samples}", "{samples}.liftover.bedGraph.tmp")),
    params:
        rname     = "format2",
        lift_file = HG38TOHG19,
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
        graph = join(working_dir, "CpG", "{samples}", "{samples}.liftover.bedGraph.tmp"),
    output:
        sort  = temp(join(working_dir, "CpG", "{samples}", "{samples}.sorted.bedGraph")),
    params:
        rname   = "format3",
        markers = GOLD_MARKERS,
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
        bedtools sort -i {input.graph} \\
            | bedtools intersect -wo -a {params.markers} -b stdin -sorted \\
            | awk '$6-$5==1 {{print $0}}' \\
            | awk 'NF{{NF-=1}};1' \\
        > {output.sort}
        """


rule aggregate_over_regions:
    input:
        sort = join(working_dir, "CpG", "{samples}", "{samples}.sorted.bedGraph"),
    output:
        tsv  = join(working_dir, "CpG", "{samples}", "{samples}.cfDNAmeInput.bedGraph"),
    params:
        rname      = "format4",
        script_dir = join(working_dir, "workflow", "scripts"),
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
        Rscript {params.script_dir}/aggregate_over_regions.R \\
            {input.sort} \\
            {output.tsv}
        """


rule cfDNAme:
    input:
        bed = join(working_dir, "CpG", "{samples}", "{samples}.cfDNAmeInput.bedGraph"),
    output:
        tsv = join(working_dir, "CpG", "{samples}", "{samples}.cfDNAmeDeconvolution.tsv"),
    params:
        rname      = "cfDNAme",
        script_dir = join(working_dir, "workflow", "scripts"),
        sampleName = "{samples}",
        reference_markers = REF_MARKERS,
        reference_IDs     = REF_IDS,
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
        Rscript {params.script_dir}/tissues_of_origin_v2.R \\
            {input.bed}  \\
            {params.reference_markers} \\
            {output.tsv} \\
            {params.reference_IDs} \\
            {params.sampleName} \\
            TRUE \\
            FALSE \\
            TRUE \\
            colon_5/hsc_5/hsc_2/spleen_1/spleen_2/spleen_3/spleen_4/
        """


rule bamsort:
    input:
        bam = join(working_dir, "bismarkAlign", "{samples}.bismark_bt2_pe.deduplicated.bam"),
    output:
        bam = temp(join(working_dir, "bismarkAlign", "bams", "{samples}.bam")),
        bai = temp(join(working_dir, "bismarkAlign", "bams", "{samples}.bam.bai")),
    params:
        rname     = "bamsort",
        reference = bisulphite_fa,
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
        samtools sort \\
            -@ {threads} \\
            --reference {params.reference} \\
            {input.bam} \\
        > {output.bam}
        samtools index \\
            -@ {threads} \\
            {output.bam}
        """


rule wgbstools:
    input:
        bam = join(working_dir, "bismarkAlign", "bams", "{samples}.bam"),
        bai = join(working_dir, "bismarkAlign", "bams", "{samples}.bam.bai"),
    output:
        pat = join(working_dir, "UXM", "{samples}.pat.gz"),
    params:
        rname      = "wgbstools",
        ref        = species,
        outdir     = join(working_dir, "UXM"),
        script_dir = join(working_dir, "workflow", "scripts"),
        atlas_bed  = config["references"][species].get("UXM_ATLAS_BED", "__no_ref_available"),
        tool_path  = str(os.path.dirname(config["tools"]["WGBS_TOOLS"])),
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
        # Add PATH to wgbstools and UXM
        export PATH=$PATH:{params.tool_path}
        wgbstools bam2pat \\
            --genome {params.ref} \\
            -L {params.atlas_bed} \\
            --out_dir {params.outdir} \\
            -@ {threads} {input.bam}
        """


rule UXM:
    input:
        pat = join(working_dir, "UXM", "{samples}.pat.gz"),
    output:
        pat = join(working_dir, "UXM", "{samples}_deconv.250.csv"),
    params:
        rname = "UXM",
        atlas_ref  = config["references"][species].get("UXM_ATLAS_REF", "__no_ref_available"),
        # UXM tools needs wbstools in the PATH,
        # assumption being made here is that 
        # both tools are in the same directory
        tool_path  = str(os.path.dirname(config["tools"]["UXM_DECONV"])),
        tmpdir    = tmp_dir,
    resources:
        mem       = allocated("mem",       "UXM", cluster),
        gres      = allocated("gres",      "UXM", cluster),
        time      = allocated("time",      "UXM", cluster),
        partition = allocated("partition", "UXM", cluster),
    threads:
        int(allocated("threads", "UXM", cluster)),
    shell:
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        module load samtools bedtools bamtools
        # Add PATH to wgbstools and UXM
        export PATH=$PATH:{params.tool_path}
        uxm deconv \\
            {input.pat} \\
            -o {output.pat} \\
            --atlas {params.atlas_ref} \\
            --tmp_dir ${{tmp}} \\
            --ignore Colon-Fibro Dermal-Fibro Gallbladder Bone-Osteob
        """
