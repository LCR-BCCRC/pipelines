#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
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
# `CFG` is a shortcut to `config["lcr-modules"]["liftover"]`
CFG = op.setup_module(
    name = "liftover",
    version = "1.2",
    subdirectories = ["inputs", "seg2bed", "liftover", "bed2seg", "outputs"])

# Define rules to be run locally when using a compute cluster
localrules:
    _liftover_input_seg,
    _liftover_seg_2_bed,
    _run_liftover,
    _liftover_sort,
    _liftover_bed_2_seg,
    _liftover_fill_segments,
    _liftover_output_seg,
    _liftover_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _liftover_input_seg:
    input:
        seg = CFG["inputs"]["sample_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.seg"
    run:
        op.relative_symlink(input.seg, output.seg)


# Convert initial seg file into bed format
rule _liftover_seg_2_bed:
    input:
        seg = str(rules._liftover_input_seg.output.seg)
    output:
        bed = CFG["dirs"]["seg2bed"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.bed",
        header = temp(CFG["dirs"]["seg2bed"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.bed.header")
    log:
        stderr = CFG["logs"]["seg2bed"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.stderr.log"
    params:
        opts = CFG["options"]["seg2bed2seg"],
        chr_colNum = CFG["options"]["chr_colNum"],
        start_colNum = CFG["options"]["start_colNum"],
        end_colNum = CFG["options"]["end_colNum"],
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.seg} 
        --output {output.bed} 
        --chromColnum {params.chr_colNum} 
        --startColnum {params.start_colNum} 
        --endColnum {params.end_colNum}
        2> {log.stderr}
        """)


def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the bed file in hg38 coordinates into hg19 coordinates
rule _run_liftover:
    input:
        native = rules._liftover_seg_2_bed.output.bed,
        chains = get_chain
    output:
        lifted = CFG["dirs"]["liftover"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.bed",
        unmapped = CFG["dirs"]["liftover"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.unmapped.bed"
    log:
        stderr = CFG["logs"]["liftover"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.stderr.log"
    params:
        mismatch = CFG["options"]["min_mismatch"]
    conda:
        CFG["conda_envs"]["liftover-366"]
    wildcard_constraints:
        chain = "hg38ToHg19|hg19ToHg38"
    shell:
        op.as_one_line("""
        liftOver -minMatch={params.mismatch}
        {input.native} {input.chains} 
        {output.lifted} {output.unmapped}
        2> {log.stderr}
        """)

# Sort liftover output
# Here, the perl line will filter out non-standard chromosomes from the output
rule _liftover_sort:
    input:
        lifted = rules._run_liftover.output.lifted
    output:
        lifted_sorted = CFG["dirs"]["liftover"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.sorted.bed"
    log:
        stderr = CFG["logs"]["liftover"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.sorted.stderr.log"
    shell:
        op.as_one_line("""
        sort -k1,1 -k2,2n -V {input.lifted} |
        perl -ne 'print if /^(chr)*[\dX]+\s.+/'
        > {output.lifted_sorted}
        2> {log.stderr}
        """)


# Convert the bed file in hg19 coordinates into seg format
rule _liftover_bed_2_seg:
    input:
        lifted_sorted = str(rules._liftover_sort.output.lifted_sorted),
        headers = str(rules._liftover_seg_2_bed.output.header)
    output:
        seg_lifted = CFG["dirs"]["bed2seg"] + "from--{seq_type}--{genome_build}/raw_segments/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["bed2seg"] + "from--{seq_type}--{genome_build}/raw_segments/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.stderr.log"
    params:
        opts = CFG["options"]["seg2bed2seg"]
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.lifted_sorted}
        --column-header {input.headers}
        --output {output.seg_lifted} 
        2> {log.stderr}
        """)


# Fill in empty segments after lifting them over
rule _liftover_fill_segments:
    input:
        seg_lifted = str(rules._liftover_bed_2_seg.output.seg_lifted)
    output:
        seg_filled = CFG["dirs"]["bed2seg"] + "from--{seq_type}--{genome_build}/filled_segments/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.filled.seg"
    log:
        stdout = CFG["logs"]["bed2seg"] + "from--{seq_type}--{genome_build}/filled_segments/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.filled.stdout.log",
        stderr = CFG["logs"]["bed2seg"] + "from--{seq_type}--{genome_build}/filled_segments/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.filled.stderr.log"
    params:
        script = CFG["options"]["fill_segments"],
        chromArm = op.switch_on_wildcard("chain", CFG["chromArm"])
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python3 {params.script}
        --input {input.seg_lifted}
        --output {output.seg_filled}
        --chromArm {params.chromArm}
        > {log.stdout}
        2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _liftover_output_seg:
    input:
        seg = str(rules._liftover_fill_segments.output.seg_filled)
    output:
        seg = CFG["dirs"]["outputs"] + "from--{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.{tool}.lifted_{chain}.seg"
    run:
        op.relative_symlink(input.seg, output.seg, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _liftover_all:
    input:
        expand(
            [
                str(rules._liftover_output_seg.output.seg)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_sample_id=CFG["runs"]["tumour_sample_id"],
            normal_sample_id=CFG["runs"]["normal_sample_id"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            pair_status=CFG["runs"]["pair_status"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=[CFG["tool"]] * len(CFG["runs"]["tumour_sample_id"]),
            chain=["hg38ToHg19" if "38" in str(x) else "hg19ToHg38" for x in CFG["runs"]["tumour_genome_build"]]
            )
            
            

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
