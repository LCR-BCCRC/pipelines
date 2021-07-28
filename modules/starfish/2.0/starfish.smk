#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob 
import string
import copy
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

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
# `CFG` is a shortcut to `config["lcr-modules"]["starfish"]`
CFG = op.setup_module(
    name = "starfish",
    version = "2.0",
    subdirectories = ["inputs", "starfish", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _starfish_input_vcf,
    _starfish_rename_output, 
    _starfish_output_vcf,
    _starfish_output_venn,
    _starfish_all,

##### GLOBAL VARIABLES #####

# All variant callers from config

CFG["inputs"]["vcf"] = dict(zip(CFG["inputs"]["names"], CFG["inputs"]["paths"]))
callers = list(CFG["inputs"]["vcf"].keys())
callers = [caller.lower() for caller in callers]

assert len(callers) < 6, (
    "Starfish can only handle a maximum of 6 input VCFs."
)

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _starfish_input_vcf:
    input:
        vcf = lambda w: config["lcr-modules"]["starfish"]["inputs"]["vcf"][w.caller] 
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{caller}.vcf.gz", 
        tbi = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{caller}.vcf.gz.tbi"
    run:
        op.absolute_symlink(input.vcf, output.vcf), 
        op.absolute_symlink(input.vcf + ".tbi", output.tbi)


# Run Starfish
rule _starfish_run:
    input:
        vcfs = expand(
            CFG["dirs"]["inputs"] + "vcf/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{caller}.vcf.gz", 
            caller = callers
            ),
        tbis = expand(
            CFG["dirs"]["inputs"] + "vcf/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{caller}.vcf.gz.tbi", 
            caller = callers
            ),
        reference = ancient(reference_files("genomes/{genome_build}/sdf")),
        starfish_script = CFG["inputs"]["starfish_script"]
    output:
        complete = touch(CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.complete"), 
        venn = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/venn.pdf"
    log:
        stdout = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish_run.stdout.log",
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish_run.stderr.log"
    params:
        opts = CFG["options"]["starfish_run"]
    conda:
        CFG["conda_envs"]["starfish"]
    threads:
        CFG["threads"]["starfish_run"]
    resources:
        **CFG["resources"]["starfish_run"]
    shell:
        op.as_one_line("""
        if [[ -e $(dirname {output.complete})/temp ]]; then 
            rm -rf $(dirname {output.complete})/temp; 
        fi
        &&
        {input.starfish_script} --sdf {input.reference} -O $(dirname {output.complete})
        --names {callers} --threads {threads}
        {params.opts} --venn --verbose
        -V {input.vcfs} > {log.stdout} 2> {log.stderr} 
        """)

checkpoint _starfish_rename_output: 
    input: 
        complete = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.complete"
    output: 
        renamed = touch(CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.renamed")
    run: 
        # Get the starfish letter names for the variant callers
        names = list(string.ascii_uppercase[0:(len(callers))])
        # Generate a dictionary of letter names and variant callers
        name_dict = {names[i]: callers[i] for i in range(len(callers))}
        outdir = os.path.dirname(input.complete)
        # Replace each alphabetically-named VCF with the name(s) of the variant caller(s)
        for src in glob.glob(outdir + "/[A-Z]*.vcf.gz*"): 
            bname = os.path.basename(src)
            lhs, rhs = bname.split(".", 1)
            lhs = "_".join(lhs)
            for key in list(name_dict.keys()): 
                lhs = lhs.replace(key, name_dict[key])
            dest = os.path.join(outdir, ".".join([lhs, rhs]))
            os.rename(src, dest) 


def _starfish_get_output(wildcards): 
    CFG = config["lcr-modules"]["starfish"]
    # Get the path to the output directory from the checkpoint _starfish_rename_output
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = os.path.dirname(checkpoints._starfish_rename_output.get(**wildcards).output.renamed)
    # Obtain the vcf_name wildcards using glob_wildcards and expand that into the vcfs object
    vcfs = expand(
        CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz",
        **wildcards,  
        vcf_name = glob_wildcards(os.path.join(checkpoint_output, "{vcf_name, [a-z]\S+}.vcf.gz")).vcf_name # Constrain to omit 2, 2+, 2-, etc. VCF files. 
        )
    return vcfs

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _starfish_output_vcf:
    input:
        vcf = CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz"
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True), 
        op.relative_symlink(input.vcf + ".tbi", output.tbi, in_module=True)


rule _starfish_output_venn: 
    input: 
        venn = str(rules._starfish_run.output.venn)
    output: 
        venn = CFG["dirs"]["outputs"] + "venn/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.venn.pdf"
    run: 
        op.relative_symlink(input.venn, output.venn, in_module=True)

def _starfish_get_output_target(wildcards): 
    CFG = config["lcr-modules"]["starfish"]
    # Get the path to the output directory from the checkpoint _starfish_rename_output
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = os.path.dirname(checkpoints._starfish_rename_output.get(**wildcards).output.renamed)
    # Obtain the vcf_name wildcards using the glob_wildcards function
    vcfs = expand(
        CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz",
        **wildcards,
        vcf_name = glob_wildcards(os.path.join(checkpoint_output, "{vcf_name}.vcf.gz")).vcf_name
        )
    return vcfs


rule _starfish_dispatch: 
    input: 
        str(rules._starfish_output_venn.output.venn), 
        _starfish_get_output_target        
    output: 
        touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _starfish_all:
    input:
        expand(
            str(rules._starfish_dispatch.output),
            zip,  
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
