#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton and Shaghayegh Soudi
# Module Author:    Laura Hilton and Shaghayegh Soudi
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import hashlib
import glob
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
# `CFG` is a shortcut to `config["lcr-modules"]["phylowgs"]`
CFG = op.setup_module(
    name = "phylowgs",
    version = "1.0",
    subdirectories = ["inputs", "maf_to_vcf", "preprocess_battenberg", "preprocess_inputs", "multievolve", "results",  "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _phylowgs_input_maf,
    _phylowgs_input_battenberg,
    _phylowgs_process_output,
    _phylowgs_output_plots,
    _phylowgs_all,
    _phylowgs_priority_ssms


# Generate a de-duplicated table of patient_ids etc.
PATIENTS = CFG["runs"][["tumour_patient_id", "normal_patient_id", "tumour_genome_build", "tumour_seq_type", "tumour_sex"]].drop_duplicates(subset = None, ignore_index = True)

# Obtain the path to the phylowgs conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")
md5hash.update(conda_prefix.encode())
f = open(CFG['conda_envs']['phylowgs'], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
PHYLO = "".join(glob.glob(conda_prefix + "/" + h[:8] + "*/share/phylowgs/"))

##### FUNCTIONS #####

# Input function to get all MAFs per patient
def get_input_mafs(wildcards):
    CFG = config["lcr-modules"]["phylowgs"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    inputs = expand(
        [
            str(rules._phylowgs_input_maf.output.maf)
        ],
        zip,
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        allow_missing = True
    )
    return(inputs)

def get_maf_cli(wildcards):
    CFG = config["lcr-modules"]["phylowgs"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    inputs = expand(
        [
            str(rules._phylowgs_input_maf.output.maf)
        ],
        zip,
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        genome_build = PATIENT["tumour_genome_build"],
        seq_type = PATIENT["tumour_seq_type"]
    )
    cli =  ",".join([str(elem) for elem in inputs])
    return(cli)

# Define the order of sample labels by time point
def order_samples(wildcards):
    CFG = config["lcr-modules"]["phylowgs"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id)
    samples = str(",".join(PATIENT.sort_values(by = ["tumour_time_point"]).tumour_sample_id.tolist()))
    return(samples)

# Expand the input files to create a command-line argument for create_phylowgs_inputs.py
def create_phylowgs_inputs_cli(wildcards):
    CFG = config["lcr-modules"]["phylowgs"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    cnvs = expand(
        "--cnvs {time_point}=" + CFG['dirs']['preprocess_battenberg'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.cnvs.txt",
        zip,
        time_point = PATIENT["tumour_time_point"],
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        seq_type = PATIENT["tumour_seq_type"],
        genome_build = PATIENT["tumour_genome_build"],
        patient_id = PATIENT["tumour_patient_id"]
    )
    vcf_types = expand(
        "--vcf-type {time_point}=mutect_smchet",
        zip,
        time_point = PATIENT["tumour_time_point"]
    )
    vcfs = expand(
        "{time_point}=" + CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf_to.vcf.gz",
        zip,
        time_point = PATIENT["tumour_time_point"],
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        seq_type = PATIENT["tumour_seq_type"],
        genome_build = PATIENT["tumour_genome_build"],
        patient_id = PATIENT["tumour_patient_id"]
    )
    cli = cnvs + vcf_types + vcfs
    cli = " ".join([str(elem) for elem in cli])
    return(cli)

# Input function to pull in input VCF and preprocessed CNV data
def create_phylowgs_inputs(wildcards):
    CFG = config["lcr-modules"]["phylowgs"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    inputs = expand(
        [
            CFG["dirs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.cnvs.txt",
            CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf_to.vcf.gz"
        ],
        zip,
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        allow_missing = True
    )
    return(inputs)

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _phylowgs_input_maf:
    input:
        maf = CFG["inputs"]["maf"],
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf",
    group: "input_maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _phylowgs_input_battenberg:
    input:
        cellularity = CFG["inputs"]["cellularity"],
        subclones = CFG["inputs"]["subclones"]
    output:
        cellularity = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.cellularity_ploidy.txt",
        subclones = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.subclones.txt"
    group: "input_battenberg"
    run:
        op.absolute_symlink(input.cellularity, output.cellularity)
        op.absolute_symlink(input.subclones, output.subclones)


rule _phylowgs_parse_battenberg:
    input:
        cellularity = str(rules._phylowgs_input_battenberg.output.cellularity),
        subclones = str(rules._phylowgs_input_battenberg.output.subclones)
    output:
        txt = CFG["dirs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.cnvs.txt"
    log:
        stderr = CFG["logs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.preprocess_battenberg.stderr.log",
        stdout = CFG["logs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.preprocess_battenberg.stdout.log"
    params:
        script = PHYLO + "parser/parse_cnvs.py"
    conda:
         CFG["conda_envs"]["phylowgs"]
    threads:
        CFG["threads"]["create_inputs"]
    resources:
        **CFG["resources"]["create_inputs"]
    group: "input_battenberg"
    shell:
        op.as_one_line("""
        cellularity=$(tail -n +2 {input.cellularity} | cut -f 1);
        python2 {params.script} -f battenberg-smchet -c $cellularity --cnv-output {output.txt} {input.subclones}
        2> {log.stderr} > {log.stdout}
        """)

# Convert the input maf file to a vcf file
rule _phylowgs_maf_to_vcf:
    input:
        maf = str(rules._phylowgs_input_maf.output.maf),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf_to.vcf")
    conda:
        CFG["conda_envs"]["vcf2maf"]
    group: "input_maf"
    shell:
        op.as_one_line("""
        maf2vcf.pl --input-maf {input.maf} --output-dir $(dirname {output.vcf}) --output-vcf {output.vcf} --ref-fasta {input.fasta}
        """)

rule _phylowgs_priority_ssms:
    input:
        mafs = get_input_mafs
    output:
        ssms = CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/coding_ssms.txt"
    params:
        noncoding = CFG["scripts"]["noncoding"]
    conda:
        CFG["conda_envs"]["coreutils"]
    shell:
        op.as_one_line("""
        grep -hvf {params.noncoding} {input.mafs} | awk '{{FS="\\t"}} {{OFS="_"}} {{print $5, $6}}' | sed 's/chr//g' > {output.ssms}
        """)

rule _phylowgs_bgzip_vcf:
    input:
        vcf = str(rules._phylowgs_maf_to_vcf.output.vcf)
    output:
        vcf = temp(CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf_to.vcf.gz"),
        tbi = CFG['dirs']['maf_to_vcf'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf_to.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools sort {input.vcf} | bcftools view -s "{wildcards.normal_id},{wildcards.tumour_id}" -i 'FMT/DP[0] > 0 && FMT/AD[0:1] > 1' -Oz -o {output.vcf} && tabix -p vcf {output.vcf}
        """)


# Preprocess vcf and battenberg inputs together
rule _phylowgs_create_inputs:
    input:
        create_phylowgs_inputs,
        priority = str(rules._phylowgs_priority_ssms.output.ssms)
    output:
        ssms = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/ssm_data.txt",
        cnvs = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/cnv_data.txt",
        params = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/params.json"
    log:
        stderr = CFG["logs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/create_phylowgs_inputs.stderr.log",
        stdout = CFG["logs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/create_phylowgs_inputs.stdout.log"
    params:
        cli = create_phylowgs_inputs_cli,
        opts = CFG["options"]["create_inputs"]["opts"],
        sex = lambda w: config["lcr-modules"]["phylowgs"]["switches"]["sex"][op.filter_samples(PATIENTS, tumour_patient_id = w.patient_id)["tumour_sex"].values[0]] if op.filter_samples(PATIENTS, tumour_patient_id = w.patient_id)["tumour_sex"].values[0] in config["lcr-modules"]["phylowgs"]["switches"]["sex"].keys() else "auto",
        script = PHYLO + "parser/create_phylowgs_inputs.py"
    conda:
        CFG["conda_envs"]["phylowgs"]
    threads:
        CFG["threads"]["create_inputs"]
    resources:
        **CFG["resources"]["create_inputs"]
    shell:
        op.as_one_line("""
        python2 {params.script}
        --output-cnvs {output.cnvs}
        --output-variants {output.ssms}
        --output-params {output.params}
        --priority-ssms {input.priority}
        --sex {params.sex}
        {params.opts}
        {params.cli}
        2> {log.stderr} > {log.stdout}
        """)

# Run multievolve to sample trees and reconstruct phylogeny
rule _phylowgs_multievolve:
    input:
        **rules._phylowgs_create_inputs.output
    output:
        trees = CFG["dirs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/trees.zip"
    log:
        stderr = CFG["logs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/multievolve.stderr.log",
        stdout = CFG["logs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/multievolve.stdout.log"
    params:
        script = PHYLO + "multievolve.py",
        opts = CFG["options"]["multievolve"]
    conda:
        CFG["conda_envs"]["phylowgs"]
    threads:
        CFG["threads"]["multievolve"]
    resources:
        **CFG["resources"]["multievolve"]
    shell:
        op.as_one_line("""
        python2 {params.script}
        {params.opts}
        -n {threads}
        -O $(dirname {output.trees})
        --ssms {input.ssms}
        --cnvs {input.cnvs}
        2> {log.stderr} > {log.stdout}
        """)

# Write the results
rule _phylowgs_write_results:
    input:
        trees = str(rules._phylowgs_multievolve.output.trees)
    output:
        muts = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.muts.json",
        summ = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.summ.json",
        mutass = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.mutass.zip"
    params:
        script = PHYLO + "write_results.py",
        opts = CFG["options"]["write_results"]
    conda:
        CFG["conda_envs"]["phylowgs"]
    threads:
        CFG["threads"]["write_results"]
    resources:
        **CFG["resources"]["write_results"]
    shell:
        op.as_one_line("""
        python2 {params.script}
        {params.opts}
        {wildcards.patient_id}
        {input.trees}
        {output.summ}.gz
        {output.muts}.gz
        {output.mutass} &&
        gunzip -f $(dirname {output.mutass})/*.gz &&
        rm -rf $(dirname {input.trees})/chain*
        """)



# Symlinks the final output files to the witness directory in preparation for HTTP browsing
rule _phylowgs_process_output:
    input:
        mafs = get_input_mafs,
        **rules._phylowgs_create_inputs.output,
        **rules._phylowgs_write_results.output,
    output:
        tree_summary = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/results/tree_summary.tsv",
        maf = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/results/merged_ssm_cluster_assignments.maf",
        cnvs = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/results/merged_cnvs_cluster_assignments.tsv",
        CCF = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/results/CCF.tsv",
        plots = directory(CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/results/plots/")
    params:
        sample_order = order_samples,
        maf_list = get_maf_cli,
        drivers = CFG['inputs']['drivers'],
        script = CFG["scripts"]["process_outputs"]
    conda:
        CFG["conda_envs"]["phylowgs_results"]
    script:
        "{params.script}"

rule _phylowgs_output_plots:
    input:
        plots = str(rules._phylowgs_process_output.output.plots)
    output:
        plots = directory(CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plots/{patient_id}")
    run:
        op.relative_symlink(input.plots, output.plots)



# Generates the target sentinels for each run, which generate the symlinks
rule _phylowgs_all:
    input:
        expand(
            [
                str(rules._phylowgs_output_plots.output.plots),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=PATIENTS["tumour_seq_type"],
            genome_build=PATIENTS["tumour_genome_build"],
            patient_id=PATIENTS["tumour_patient_id"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
