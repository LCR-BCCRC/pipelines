#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob

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
# `CFG` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG = op.setup_module(
    name = "battenberg",
    version = "1.0",
    subdirectories = ["inputs", "battenberg", "outputs"],
)

#this preserves the variable when using lambda functions
_battenberg_CFG = CFG

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _install_battenberg,
    _battenberg_input_bam,
    _battenberg_output_seg,
    _battenberg_to_igv_seg,
    _battenberg_cleanup,
    _battenberg_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam + ".bai", output.bai)

# Installs the Battenberg R dependencies and associated software (impute2, alleleCounter)
# Currently I think this rule has to be run twice for it to work properly because the conda environment is created here. 
# I am open to suggestions for how to get around this.
rule _install_battenberg:
    output:
        complete = "config/envs/battenberg_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["battenberg"]
    shell:
        """
        R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")' && ##move some of this to config?
        R -q -e 'devtools::install_github("morinlab/battenberg")' &&              ##move some of this to config?
        touch {output.complete}"""


# This rule runs the entire Battenberg pipeline. Eventually we may want to set this rule up to allow re-starting
# of partially completed jobs (e.g. if they run out of RAM and are killed by the cluster, they can automatically retry)
# TODO: this rule needs to be modified to rely on reference_files and allow setup (downloading) of the Battenberg references
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        tumour_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        normal_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam.bai",
        installed = "config/envs/battenberg_dependencies_installed.success"
    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        ac=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"),
        mb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab"),
        mlrg=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab"),
        mlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab"),
        nlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab"),
        nb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab")
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        reference_path = lambda w: _battenberg_CFG["reference_path"][w.genome_build],
        script = CFG["inputs"]["battenberg_script"],
        calc_sex_status = CFG["inputs"]["calc_sex_status"],
        x_chrom = lambda w: _battenberg_CFG["options"]["x_chrom"][w.genome_build],
        y_chrom = lambda w: _battenberg_CFG["options"]["y_chrom"][w.genome_build],
        chr_prefixed = lambda w: _battenberg_CFG["options"]["chr_prefixed_reference"][w.genome_build],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    shell:
        """sex=$({params.calc_sex_status} {input.normal_bam} {params.x_chrom} {params.y_chrom} | tail -1 | awk '{{if( $3/$2 > 0.1) print "male"; else print "female"}}') ;"""
        "Rscript {params.script} -t {wildcards.tumour_id} "
        "-n {wildcards.normal_id} --tb {input.tumour_bam} --nb {input.normal_bam}  "
        "-o {params.out_dir} --sex $sex --reference {params.reference_path} {params.chr_prefixed} --cpu {threads} > {log.stdout} 2> {log.stderr} "


# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        sub = rules._run_battenberg.output.sub,
        cnv2igv = CFG["inputs"]["cnv2igv"]
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    shell:
        op.as_one_line("""
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2> {log.stderr}
        """)


#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based cleaning in the output directory
rule _battenberg_cleanup:
    input:
        rules._battenberg_to_igv_seg.output.seg
    output:
        complete = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
    shell:
        op.as_one_line("""
        d=$(dirname {output});
        rm $d/*impute_input* &&
        rm $d/*alleleFrequencies* &&
        rm $d/*aplotype* &&
        rm $d/*BAFsegmented* && 
        touch {output.complete}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience

rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        sub = rules._run_battenberg.output.sub
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.igv.seg",
        sub = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.txt"
    params: 
        batt_dir = CFG["dirs"]["battenberg"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        png_dir = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}"
    run:
        plots = glob.glob(params.batt_dir + "/*.png")
        for png in plots:
            bn = os.path.basename(png)
            op.relative_symlink(png, params.png_dir + "/" + bn, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)
        op.relative_symlink(input.sub, output.sub, in_module = True)

# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [
                rules._run_battenberg.output.sub,
                rules._battenberg_output_seg.output.seg,
                rules._battenberg_cleanup.output.complete
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
