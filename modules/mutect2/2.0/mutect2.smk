#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import inspect

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
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mutect2"]`
CFG = op.setup_module(
    name = "mutect2",
    version = "2.1",
    subdirectories = ["inputs", "mutect2", "filter", "passed", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutect2_input_bam,
    _mutect2_input_chrs,
    _mutect2_get_sm,
    _mutect2_output_vcf,
    _mutect2_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutect2_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)

rule _mutect2_dummy_positions:
    # creates a dummy vcf if users do not specify candidateSmallIndels file
    output:
        touch(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.dummy.vcf")


# Symlink chromosomes used for parallelization
checkpoint _mutect2_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.chrs, output.chrs)


# Retrieves from SM tag from BAM and writes to file
rule _mutect2_get_sm:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    output:
        sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}_sm.txt",
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}_mutect2_get_sm.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        "samtools view -H {input.bam} | grep '^@RG' | "
        r"sed 's/.*SM:\([^\t]*\).*/\1/g'"" | uniq > {output.sm} 2> {log.stderr}"

# This generates a command-line argument for the Mutect2 functions by combining the 
# input intervals file with the intervals arguments specified in the config. 
# The first function loads the wildcard-containing file path and additional args from the config. 
# The second replaces wildcards with those used in the rule. 
def _mutect2_get_interval_cli_arg(
    vcf_in = config["lcr-modules"]["mutect2"]["inputs"]["candidate_positions"], 
    interval_arg_in = config["lcr-modules"]["mutect2"]["options"]["mutect2_interval_rules"]
):
    def _mutect2_get_interval_cli_custom(wildcards, input):
        if vcf_in:
            param = f"-L {input.candidate_positions} {interval_arg_in}"
        else: 
            param = ""
        return param
    return _mutect2_get_interval_cli_custom

# Launces Mutect2 in matched and unmatched mode
rule _mutect2_run_matched_unmatched:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        normal_sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{normal_id}_sm.txt", 
        pon = reference_files("genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.gz"), 
        candidate_positions = CFG["inputs"]["candidate_positions"] if CFG["inputs"]["candidate_positions"] else str(rules._mutect2_dummy_positions.output)
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.tbi"),
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.stats"), 
        f1r2 = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.f1r2.tar.gz")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stderr.log"
    resources:
        **CFG["resources"]["mutect2_run"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        opts = CFG["options"]["mutect2_run"], 
        interval_arg = lambda w: get_capture_space(w.tumour_id, w.normal_id, w.genome_build, w.seq_type, "interval_list", default = _mutect2_get_interval_cli_arg())
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "matched|unmatched"
    shell:
        op.as_one_line("""
        gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" {params.opts} 
        -I {input.tumour_bam} -I {input.normal_bam}
        -R {input.fasta} -normal "$(cat {input.normal_sm})" -O {output.vcf}
        --germline-resource {input.gnomad} 
        -L {wildcards.chrom} -L {params.interval_arg}
        -pon {input.pon} --f1r2-tar-gz {output.f1r2}
        > {log.stdout} 2> {log.stderr}
        """)


# Launches Mutect2 in no normal mode
rule _mutect2_run_no_normal:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        pon = reference_files("genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.gz"), 
        candidate_positions = CFG["inputs"]["candidate_positions"] if CFG["inputs"]["candidate_positions"] else str(rules._mutect2_dummy_positions.output)
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.tbi"),
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.stats"), 
        f1r2 = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.f1r2.tar.gz"
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stderr.log"
    resources:
        **CFG["resources"]["mutect2_run"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2_run"], 
        interval_arg = _mutect2_get_interval_cli_arg 
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "no_normal"
    shell:
        op.as_one_line("""
        gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" 
        {params.opts} -I {input.tumour_bam} -R {input.fasta} 
        -O {output.vcf} --germline-resource {input.gnomad} 
        -L {wildcards.chrom} {params.interval_arg}
        -pon {input.pon} --f1r2-tar-gz {output.f1r2}
        > {log.stdout} 2> {log.stderr}
        """)


def _mutect2_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _mutect2_get_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


# Merge chromosome mutect2 VCFs from the same sample
rule _mutect2_merge_vcfs:
    input:
        vcf = _mutect2_get_chr_vcfs,
        tbi = _mutect2_get_chr_tbis
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["mutect2_merge_vcfs"]
    resources:
        **CFG["resources"]["mutect2_merge_vcfs"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -m {params.mem_mb}M -O z -o {output.vcf} 2>> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)


def _mutect2_get_chr_stats(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    stats = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz.stats",
        chrom = chrs
    )
    return(stats)


# Merge chromosome mutect2 stats for FilterMutectCalls rule
rule _mutect2_merge_stats:
    input:
        stat = _mutect2_get_chr_stats
    output:
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz.stats")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_stats.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_stats.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for i in {input.stat}; do echo -n "-stats $i "; done)
        -O {output.stat} > {log.stdout} 2> {log.stderr}
        """)

# Learn read orientation model

def _mutect2_get_chr_f1r2(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    f1r2 = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.f1r2.tar.gz",
        chrom = chrs
    )
    return(f1r2)

rule _mutect2_learn_orient_model: 
    input: 
        f1r2 = _mutect2_get_chr_f1r2
    output:
        model =  CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/read-orientation-model.tar.gz"
    log: 
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/read-orientation-model.stdout.log", 
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/read-orientation-model.stderr.log"
    resources: 
        **CFG["resources"]["mutect2_f1r2"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads: 
        CFG["threads"]["mutect2_f1r2"]
    shell: 
        op.as_one_line("""
        inputs=$(for input in {input.f1r2}; do printf -- "-I $input "; done);
        gatk LearnReadOrientationModel 
        --java-options "-Xmx{params.mem_mb}m" 
        $inputs -O {output.model}
        > {log.stdout} 2> {log.stderr}
        """)

# Get pileup summaries
rule _mutect2_pileup_summaries: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam", 
        snps = reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz"), 
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        pileup = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pileupSummary.table"
    log: 
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pileupSummary.stdout.log", 
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pileupSummary.stderr.log"
    resources: 
        **CFG["resources"]["mutect2_pileupsummaries"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads: 
        CFG["threads"]["mutect2_pileupsummaries"]
    shell: 
        op.as_one_line("""
        gatk GetPileupSummaries 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.tumour_bam}
            -R {input.fasta} 
            -V {input.snps}
            -L {input.snps}
            -O {output.pileup}
            > {log.stdout} 2> {log.stderr}
        """)

# Calculate contamination  
rule _mutect2_calc_contamination: 
    input: 
        pileup = str(rules._mutect2_pileup_summaries.output.pileup)
    output: 
        segments = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/segments.table", 
        contamination = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/contamination.table"
    log: 
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/contamination.stdout.log", 
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/contamination.stderr.log"
    resources: 
        **CFG["resources"]["mutect2_contamination"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads: 
        CFG["threads"]["mutect2_contamination"]
    shell: 
        op.as_one_line("""
        gatk CalculateContamination 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.pileup}
            -tumor-segmentation {output.segments}
            -O {output.contamination}
            > {log.stdout} 2> {log.stderr}
        """)
    
# Marks variants filtered or PASS annotations
rule _mutect2_filter:
    input:
        vcf = str(rules._mutect2_merge_vcfs.output.vcf),
        tbi = str(rules._mutect2_merge_vcfs.output.tbi),
        stat = str(rules._mutect2_merge_stats.output.stat),
        segments = str(rules._mutect2_calc_contamination.output.segments), 
        contamination = str(rules._mutect2_calc_contamination.output.contamination), 
        model = str(rules._mutect2_learn_orient_model.output.model),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.unfilt.vcf.gz"
    log:
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter.stdout.log",
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter.stderr.log"
    resources:
        **CFG["resources"]["mutect2_filter"]
    params:
        opts = CFG["options"]["mutect2_filter"], 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_filter"]
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls --java-options "-Xmx{params.mem_mb}m" 
            {params.opts} 
            -V {input.vcf} 
            -R {input.fasta}
            --tumor-segmentation {input.segments}
            --contamination-table {input.contamination}
            --ob-priors {input.model}
            -O {output.vcf} 
            > {log.stdout} 2> {log.stderr}
        """)


# Filters for PASS variants
rule _mutect2_filter_passed:
    input:
        vcf = str(rules._mutect2_filter.output.vcf)
    output:
        vcf = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.passed.vcf.gz",
        tbi = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.passed.vcf.gz.tbi"
    params:
        opts = CFG["options"]["mutect2_filter_passed"]
    log:
        stderr = CFG["logs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter_passed.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["mutect2_passed"]
    resources:
        **CFG["resources"]["mutect2_passed"]
    shell:
        op.as_one_line(""" 
        bcftools view {params.opts} -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutect2_output_vcf:
    input:
        vcf = str(rules._mutect2_filter_passed.output.vcf),
        tbi = str(rules._mutect2_filter_passed.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True)
        op.relative_symlink(input.tbi, output.tbi, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutect2_all:
    input:
        expand(
            [
                str(rules._mutect2_output_vcf.output.vcf),
                str(rules._mutect2_output_vcf.output.tbi)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
