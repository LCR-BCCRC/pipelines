#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package oncopipe and inspect packages. 
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 


# Setup module and store module-specific configuration in `CFG_SLMS3`
# `CFG_SLMS3` is a shortcut to `config["lcr-modules"]["slms_3"]`
# This is not assigned to `CFG` as in other modules because the sub-modules here will use `CFG`. 
CFG_SLMS3 = op.setup_module(
    name = "slms_3",
    version = "1.0",
    subdirectories = ["inputs", "strelka_gnomad", "lofreq_gnomad", "strelka_lofreq_union", "sage_gnomad", "mutect2_depth_filt", "union", "outputs"],
)

# Ensure each of the submodule versions meets requirements. 

assert float(CFG_SLMS3["module_versions"]["manta"]) >= 2.0, (
    f"The current manta module version is {CFG_SLMS3['module_versions']['manta']}. "
    "SLMS-3 requires manta module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["starfish"]) >= 2.0, (
    f"The current starfish module version is {CFG_SLMS3['module_versions']['starfish']}. "
    "SLMS-3 requires starfish module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["mutect2"]) >= 2.0, (
    f"The current mutect2 module version is {CFG_SLMS3['module_versions']['mutect2']}. "
    "SLMS-3 requires mutect2 module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["strelka"]) >= 1.1, (
    f"The current strelka module version is {CFG_SLMS3['module_versions']['strelka']}. "
    "SLMS-3 requires strelka module version 1.1 or higher. "
)

# Define rules to be run locally when using a compute cluster
localrules:
    _slms_3_input_sage_vcf,
    _slms_3_input_strelka_vcf,
    _slms_3_input_mutect_vcf, 
    _slms_3_input_lofreq_vcf, 
    _slms_3_mutect2_samples_table, 
    _slms_3_output_vcf, 
    _slms_3_all,

##### MODULE CONFIGFILES #####

# Load default module configs and update with values from slms_3 config

for tool in CFG_SLMS3["module_versions"].keys():
    lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"
    configfile: lcr_modules + "/" +  tool + "/" + CFG_SLMS3["module_versions"][tool] + "/config/default.yaml"
    snakemake.utils.update_config(
        config["lcr-modules"][tool], 
        config["lcr-modules"]["slms_3"][tool]
    )
    config["lcr-modules"][tool]["inputs"]["sample_bam"] = CFG_SLMS3["inputs"]["sample_bam"]
    config["lcr-modules"][tool]["inputs"]["sample_bai"] = CFG_SLMS3["inputs"]["sample_bai"]    


##### FIRST PASS VARIANT CALLING MODULE SNAKEFILES #####

# Load Manta first
include: "../../manta/" + CFG_SLMS3["module_versions"]["manta"] + "/manta.smk"

# Update Strelka config to use Manta Candidate Small Indels    
config["lcr-modules"]["strelka"]["inputs"]["candidate_small_indels"] = expand(str(rules._manta_output_vcf.output.vcf), vcf_name = "candidateSmallIndels", allow_missing=True)

# Load Strelka, SAGE, and Lofreq
include: "../../strelka/" + CFG_SLMS3["module_versions"]["strelka"] + "/strelka.smk"
include: "../../sage/" + CFG_SLMS3["module_versions"]["sage"] + "/sage.smk"
include: "../../lofreq/" + CFG_SLMS3["module_versions"]["lofreq"] + "/lofreq.smk"
 
##### FIRST PASS VARIANT CALLING RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

rule _slms_3_input_strelka_vcf: 
    input:
        vcf = str(rules._strelka_output_filtered_vcf.output.vcf), 
        tbi = str(rules._strelka_output_filtered_vcf.output.vcf_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)

rule _slms_3_input_sage_vcf: 
    input:
        vcf = str(rules._sage_output_vcf.output.combined), 
        tbi = str(rules._sage_output_vcf.output.combined_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)

rule _slms_3_input_lofreq_vcf: 
    input:
        vcf = str(rules._lofreq_output_vcf.output.vcf_all), 
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.vcf + ".tbi", output.tbi, in_module = True)

# Annotate Strelka VCF and remove common GnomAD variants

rule _slms_3_annotate_strelka_gnomad:
    input:
        vcf = str(rules._slms_3_input_strelka_vcf.output.vcf),
        tbi = str(rules._slms_3_input_strelka_vcf.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_gnomad.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    threads:
        CFG_SLMS3["threads"]["strelka_gnomad"]
    resources:
        **CFG_SLMS3["resources"]["strelka_gnomad"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        bcftools view -i 'INFO/AF < 0.0001 && INFO/SomaticEVS >= 10 && FMT/DP[1] >= 10' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Annotate LoFreq VCF and remove common GnomAD variants

rule _slms_3_annotate_lofreq_gnomad:
    input: 
        vcf = str(rules._slms_3_input_lofreq_vcf.output.vcf), 
        tbi = str(rules._slms_3_input_lofreq_vcf.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output: 
        vcf = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_gnomad.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["lofreq_gnomad"]
    threads: 
        CFG_SLMS3["threads"]["lofreq_gnomad"]
    shell: 
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        bcftools view -i 'INFO/AF < 0.0001' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Create a union of Strelka and Lofreq as candidates for Mutect2

rule _slms_3_strelka_lofreq_union: 
    input: 
        strelka = str(rules._slms_3_annotate_strelka_gnomad.output.vcf),
        strelka_tbi = str(rules._slms_3_annotate_strelka_gnomad.output.tbi), 
        lofreq = str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
        lofreq_tbi = str(rules._slms_3_annotate_lofreq_gnomad.output.tbi)
    output: 
        vcf = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_lofreq_union_gnomad.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["strelka_lofreq_union"]
    threads: 
        CFG_SLMS3["threads"]["strelka_lofreq_union"]
    shell: 
        op.as_one_line("""
        bcftools merge --force-samples -m all --threads {threads} 
        -Oz -o {output.vcf} {input.strelka} {input.lofreq} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Modify the SAGE VCF file to use `NORMAL` and `TUMOR` as sample IDs
rule _slms_3_annotate_sage_gnomad: 
    input: 
        vcf = str(rules._slms_3_input_sage_vcf.output.vcf), 
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output: 
        vcf = CFG_SLMS3["dirs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sage.renamed.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sage.renamed.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["sage_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_lofreq_union_gnomad.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["sage_gnomad"]
    threads: 
        CFG_SLMS3["threads"]["sage_gnomad"]
    shell: 
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        sed 's/{wildcards.tumour_id}/TUMOR/g' | 
        sed 's/{wildcards.normal_id}/NORMAL/g' | 
        bcftools view -s "NORMAL,TUMOR" -i 'INFO/AF < 0.0001' -Oz -o {output.vcf} 2> {log.stderr} 
        &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


##### SECOND PASS VARIANT CALLING MODULE SNAKEFILES #####

# Upate the Mutect2 config with the _strelka_lofreq_vcf

config["lcr-modules"]["mutect2"]["inputs"]["candidate_positions"] = str(rules._slms_3_strelka_lofreq_union.output.vcf)

include: "../../mutect2/" + CFG_SLMS3["module_versions"]["mutect2"] + "/mutect2.smk"

##### SECOND PASS VARIANT CALLING RULES #####

rule _slms_3_input_mutect_vcf: 
    input:
        vcf = str(rules._mutect2_output_vcf.output.vcf), 
        tbi = str(rules._mutect2_output_vcf.output.tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Depth/VAF filter the Mutect2 output. This must be completed here and not in the mutect2 module so that the TUMOR column can be properly indexed. 
rule _slms_3_mutect2_samples_table: 
    output: 
        table = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "samples.txt"
    shell: 
        """echo "TUMOR" > {output.table}"""

rule _slms_3_mutect2_depth_filt: 
    input: 
        vcf = str(rules._slms_3_input_mutect_vcf.output.vcf), 
        table = str(rules._slms_3_mutect2_samples_table.output.table)
    output: 
        vcf = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_depth_filt.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["mutect2_depth_filt"]
    threads: 
        CFG_SLMS3["threads"]["mutect2_depth_filt"]
    shell: 
        op.as_one_line("""
        tsamp=$(zgrep "##tumor_sample=" {input.vcf} | sed 's|##tumor_sample=||g');
        nsamp=$(zgrep "##normal_sample=" {input.vcf} | sed 's|##normal_sample=||g');
        bcftools view {input.vcf} | 
        sed "s|$tsamp|TUMOR|g" | sed "s|$nsamp|NORMAL|g" |  
        bcftools view  -s "NORMAL,TUMOR" -i 'FMT/DP[@{input.table}] >= 10 && FMT/AD[@{input.table}:1] >= 4 && FMT/AF[@{input.table}:0] >= 0.1' 
        -Oz -o {output.vcf} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Update Starfish config to use outputs generated here

snakemake.utils.update_config(config["lcr-modules"]["starfish"], {
    "dirs": {
        "_parent": CFG_SLMS3["dirs"]["_parent"] + "starfish-" + CFG_SLMS3["module_versions"]["starfish"]
    },
    "inputs": {
        "names": 
            ["mutect", "sage", "strelka", "lofreq"], 
        "paths": 
            [
                str(rules._slms_3_mutect2_depth_filt.output.vcf),
                str(rules._slms_3_annotate_sage_gnomad.output.vcf),
                str(rules._slms_3_annotate_strelka_gnomad.output.vcf), 
                str(rules._slms_3_annotate_lofreq_gnomad.output.vcf)
            ]
    }
})

include: "../../starfish/" + CFG_SLMS3["module_versions"]["starfish"] + "/starfish.smk"

# Create a final union VCF file summarizing which variants were called by which variant callers
# First rename the sample columns in each VCF as TUMOR_{caller} or NORMAL_{caller}

rule _slms_3_rename_samples_all: 
    input: 
        vcf = lambda w: config["lcr-modules"]["starfish"]["inputs"]["vcf"][w.caller]
    output: 
        vcf = temp(CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.tmp.vcf.gz"),
        tbi = temp(CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.tmp.vcf.gz.tbi"), 
        samples = temp(CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}.samples.txt")
    log:
        stderr = CFG_SLMS3["logs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/rename_samples_{caller}.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["rename_all"]
    threads: 
        CFG_SLMS3["threads"]["rename_all"]
    shell: 
        op.as_one_line("""
        printf "TUMOR TUMOR_{wildcards.caller}\\nNORMAL NORMAL_{wildcards.caller}\\n" > {output.samples} 
        && 
        bcftools reheader -s {output.samples} -o {output.vcf} {input.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

rule _slms_3_union_vcf: 
    input: 
        vcf = expand(
            rules._slms_3_rename_samples_all.output.vcf, 
            caller = config["lcr-modules"]["starfish"]["inputs"]["names"], 
            allow_missing = True
        ), 
        tbi = expand(
            rules._slms_3_rename_samples_all.output.tbi, 
            caller = config["lcr-modules"]["starfish"]["inputs"]["names"], 
            allow_missing = True
        )
    output: 
        vcf = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/union.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["union_vcf"]
    threads: 
        CFG_SLMS3["threads"]["union_vcf"]
    shell: 
        op.as_one_line("""
        bcftools merge --threads {threads} -m both -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)
    

# Symlinks the final output files into the module results directory (under '99-outputs/')


rule _slms_3_output_vcf:
    input:
        union_vcf = str(rules._slms_3_union_vcf.output.vcf), 
        union_tbi = str(rules._slms_3_union_vcf.output.tbi),
        dispatched = str(rules._starfish_dispatch.output)
    output:
        isec_vcf = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz", 
        isec_tbi = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz.tbi", 
        union_vcf = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.union.vcf.gz",
        union_tbi = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.union.vcf.gz.tbi"
    params: 
        vcf_in = config["lcr-modules"]["starfish"]["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.3+.vcf.gz", 
        rm_files = CFG_SLMS3["options"]["cleanup_vcfs"], 
        files_to_rm = [
            str(rules._slms_3_annotate_strelka_gnomad.output.vcf), 
            str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
            str(rules._slms_3_strelka_lofreq_union.output.vcf), 
            str(rules._slms_3_mutect2_depth_filt.output.vcf), 
            str(rules._slms_3_annotate_sage_gnomad.output.vcf)
        ] 
    run:
        op.relative_symlink(params.vcf_in, output.isec_vcf, in_module = True)
        op.relative_symlink(params.vcf_in + ".tbi", output.isec_tbi, in_module = True)
        op.relative_symlink(input.union_vcf, output.union_vcf, in_module = True)
        op.relative_symlink(input.union_tbi, output.union_tbi, in_module = True)
        if params.rm_files: 
            for file in params.files_to_rm: 
                if os.path.exists(file): 
                    os.remove(file)




# Generates the target sentinels for each run, which generate the symlinks
rule _slms_3_all:
    input:
        rules._manta_all.input,
        rules._strelka_all.input, 
        rules._lofreq_all.input, 
        rules._sage_all.input,
        rules._mutect2_all.input,
        rules._starfish_all.input, 
        expand(
            [
                str(rules._slms_3_output_vcf.output.isec_vcf),
                str(rules._slms_3_output_vcf.output.union_vcf)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG_SLMS3["runs"]["tumour_seq_type"],
            genome_build=CFG_SLMS3["runs"]["tumour_genome_build"],
            tumour_id=CFG_SLMS3["runs"]["tumour_sample_id"],
            normal_id=CFG_SLMS3["runs"]["normal_sample_id"],
            pair_status=CFG_SLMS3["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
# op.cleanup_module(CFG_SLMS3)
