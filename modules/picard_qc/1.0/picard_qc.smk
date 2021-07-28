#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["picard_qc"]`
CFG = op.setup_module(
    name = "picard_qc",
    version = "1.0",
    subdirectories = ["inputs", "metrics", "merged_metrics", "outputs"]
)


# Define rules to be run locally when using a compute cluster
localrules:
    _picard_qc_input_bam,
    _picard_qc_rrna_int,
    _picard_qc_merge_metrics,
    _picard_qc_merged_output,
    _picard_qc_flagstats_output,
    _picard_qc_merged_dispatch,
    _picard_qc_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _picard_qc_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"],
        sample_bai = CFG["inputs"]["sample_bai"]
    output:
        sample_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        sample_bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.sample_bam, output.sample_bam)
        op.absolute_symlink(input.sample_bai, output.sample_bai)


rule _picard_qc_alignment_summary:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary_metrics",
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary.stderr.log"
    params:
        opts = CFG["options"]["alignment_summary"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["alignment_summary"]
    resources:
        mem_mb = CFG["mem_mb"]["alignment_summary"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectAlignmentSummaryMetrics
        {params.opts} 
        I={input.bam} O={output.metrics} R={input.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_insert_size:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/insert_size_metrics",
        histogram = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/insert_size_histogram.pdf"
    log:
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/insert_size.stderr.log",
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/insert_size.stdout.log"
    params:
        opts = CFG["options"]["insert_size"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["insert_size"]
    resources:
        mem_mb = CFG["mem_mb"]["insert_size"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectInsertSizeMetrics
        {params.opts} 
        I={input.bam} O={output.metrics} H={output.histogram}
        > {log.stdout} 2> {log.stderr}
        """)


#capture metrics
rule _picard_qc_hs_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        interval = op.switch_on_column("genome_build", CFG["samples"], CFG["switches"]["capture_intervals"], match_on = "sample")
    output:
        hs = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics",
        intervals = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/interval_hs_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics.stderr.log"
    params:
        opts = CFG["options"]["hs_metrics"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["hs_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["hs_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectHsMetrics {params.opts} 
        I={input.bam} O={output.hs} PER_TARGET_COVERAGE={output.intervals} 
        R={input.fasta} TI={input.interval} BI={input.interval} 
        > {log.stdout} 2> {log.stderr}
        """)

# mRNA metrics
rule _picard_qc_rrna_int:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        rrna_int = reference_files("genomes/{genome_build}/rrna_intervals/rRNA_int_gencode-33.txt")
    output:
        sample_rrna = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rrna_int_list"
    log:
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rrna_int.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        op.as_one_line("""
        samtools view -H {input.bam} | 
        grep -v "@RG" | cat - {input.rrna_int} 
        > {output.sample_rrna}
        2> {log.stderr}
        """)


rule _picard_qc_rnaseq_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
        sample_rrna = rules._picard_qc_rrna_int.output.sample_rrna,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        refFlat = reference_files("genomes/{genome_build}/annotations/refFlat_gencode-33.txt")
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics.stderr.log"
    params:
        opts = CFG["options"]["rnaseq_metrics"]["base"],
        strand = op.switch_on_column("strand", CFG["samples"], CFG["switches"]["rnaseq_metrics"], match_on = "sample")
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["rnaseq_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["rnaseq_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectRnaSeqMetrics {params.opts} 
        I={input.bam} O={output.metrics} 
        R={input.fasta} 
        STRAND_SPECIFICITY={params.strand} 
        REF_FLAT={input.refFlat}
        RIBOSOMAL_INTERVALS={input.sample_rrna}
        CHART_OUTPUT={output.metrics}.pdf 
        > {log.stdout} 2> {log.stderr}
        """)

# genome metrics
rule _picard_qc_wgs_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics.stderr.log"
    params:
        opts = CFG["options"]["wgs_metrics"],
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["wgs_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["wgs_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectWgsMetrics {params.opts} 
        I={input.bam} O={output.metrics} 
        R={input.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)


def _get_sample_metrics(metrics_dir):
    DIR = metrics_dir
    CFG = config["lcr-modules"]["picard_qc"]
    def _get_sample_metrics_custom(wildcards):
        fsample = op.filter_samples(CFG["samples"], seq_type=wildcards.seq_type, genome_build=wildcards.genome_build)
        return expand("{dir}{seq_type}--{genome_build}/{sample_id}/{metrics}", dir = DIR, sample_id = list(fsample["sample_id"]), **wildcards)
    return _get_sample_metrics_custom


rule _picard_qc_merge_metrics:
    input: 
        _get_sample_metrics(metrics_dir = CFG["dirs"]["metrics"])
    output: 
        metrics = CFG["dirs"]["merged_metrics"] + "{seq_type}--{genome_build}/all.{metrics}.txt"
    run:
        samples = input
        with open(output[0], "w") as out:
            for i in range(0,len(samples)):
                s_id = samples[i].split("/")[-2]
                with open(samples[i], "r") as f:
                    data = [l for l in f.readlines() if not (l.startswith('#') or l == '\n')]
                    if i == 0:
                        header = "SAMPLEID\t" + "\t".join(data[0].split("\t")[0:])
                        out.write(header)
                    if wildcards.metrics == "alignment_summary_metrics":
                        line = s_id + "\t" + "\t".join(data[3].split("\t")[0:])
                    else:
                        line = s_id + "\t" + "\t".join(data[1].split("\t")[0:])
                    out.write(line)
                    f.close()
        out.close()


rule _picard_qc_flagstats:
    input:
        bam = rules._picard_qc_input_bam.output.sample_bam
    output:
        flagstats = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/flagstats"
    log:
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/flagstats.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["flagstats"]
    resources:
        mem_mb = CFG["mem_mb"]["flagstats"]
    shell:
        op.as_one_line("""
        samtools flagstat -@ {threads}
        {input.bam} > {output.flagstats} 
        2> {log.stderr}
        """)


rule _picard_qc_merged_output:
    input:
        metrics = rules._picard_qc_merge_metrics.output.metrics
    output:
        metrics = CFG["dirs"]["outputs"] + "merged_metrics/{seq_type}--{genome_build}/all.{metrics}.txt"
    run:
        op.relative_symlink(input.metrics, output.metrics, in_module=True)


rule _picard_qc_flagstats_output:
    input:
        flagstats = rules._picard_qc_flagstats.output.flagstats
    output:
        flagstats = CFG["dirs"]["outputs"] + "flagstats/{seq_type}--{genome_build}/{sample_id}.flagstats"
    run:
        op.relative_symlink(input.flagstats, output.flagstats, in_module=True)


def _get_picard_qc_files(wildcards):
    # base metrics to calculate for all seq_type
    base = ["alignment_summary_metrics", "insert_size_metrics"]
    
    # get seq_type specific metrics
    if wildcards.seq_type == 'mrna':
        m = base + ["rnaseq_metrics"]
    elif wildcards.seq_type == 'genome':
        m = base + ["wgs_metrics"]
    elif wildcards.seq_type == 'capture':
        m = base + ["hs_metrics", "interval_hs_metrics"]
    else:
        m = base

    targets = expand(rules._picard_qc_merged_output.output.metrics, **wildcards, metrics = m)

    return targets

 
rule _picard_qc_merged_dispatch:
    input:
        _get_picard_qc_files
    output:
        touch(CFG["dirs"]["outputs"] + "merged_metrics/{seq_type}--{genome_build}/merged.dispatched")


rule _picard_qc_all:
    input: 
        expand(rules._picard_qc_merged_dispatch.output, zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"]),
        expand(rules._picard_qc_flagstats_output.output.flagstats, zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
