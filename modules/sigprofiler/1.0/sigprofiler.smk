#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import os.path

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
# `CFG` is a shortcut to `config["lcr-modules"]["sigprofiler"]`
CFG = op.setup_module(
    name = "sigprofiler",
    version = "1.0",
    subdirectories = ["inputs", "estimate", "extract", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sigprofiler_input_maf,
    _sigprofiler_output_tsv,
    _sigprofiler_all,


##### FUNCTIONS #####

def get_dir(wildcards):
    if "ID" in wildcards.type:
        topdir = 'ID'
    elif "SBS" in wildcards.type:
        topdir = 'SBS'
    elif "DBS" in wildcards.type:
        topdir = 'DBS'

    cc = config["lcr-modules"]["sigprofiler"]
    mat = cc["dirs"]["inputs"]+f"matrices/{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.sample_set}/output/{topdir}/{wildcards.sample_set}.{wildcards.type}.all"
    ret = {'script' : cc["inputs"]["extractor"], 'mat' : mat}
    return(ret)

max_sigs = {
    "SBS6"    : 10,
    "SBS18"   : 10,
    "SBS24"   : 20,
    "SBS96"   : 20,
    "SBS288"  : 25,
    "SBS384"  : 25,
    "SBS4608" : 30,
    "SBS1536" : 30,
    "SBS6144" : 35,
    "ID28"    : 5,
    "ID83"    : 5,
    "ID96"    : 10,
    "ID332"   : 10,
    "ID415"   : 15,
    "ID8628"  : 20,
    "DBS78"   : 15,
    "DBS150"  : 15,
    "DBS186"  : 15,
    "DBS1248" : 20,
    "DBS2400" : 20,
    "DBS2976" : 20
}

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sigprofiler_input_maf:
    input:
        maf = CFG["inputs"]["maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "matrices/{seq_type}--{genome_build}/{sample_set}/{sample_set}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)

# Generates sample by k-mer context matrices from MAF
rule _sigprofiler_run_generator:
    input:
        mg = reference_files("genomes/{genome_build}/sigprofiler_genomes/{genome_build}.installed"),
        script = CFG["inputs"]["generator"],
        maf = str(rules._sigprofiler_input_maf.output.maf)
    output:
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/SBS/{{sample_set}}.{type}.all", type = ['SBS6','SBS18','SBS24','SBS96','SBS288','SBS384','SBS4608','SBS1536','SBS6144']),
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/ID/{{sample_set}}.{type}.all", type = ['ID28','ID83','ID96','ID332','ID415','ID8628']),
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/DBS/{{sample_set}}.{type}.all", type = ['DBS78','DBS150','DBS186','DBS1248','DBS2400','DBS2976'])
    log:
        stdout = CFG["logs"]["inputs"] + "{seq_type}--{genome_build}/{sample_set}/generator.stdout.log",
        stderr = CFG["logs"]["inputs"] + "{seq_type}--{genome_build}/{sample_set}/generator.stderr.log"
    params:
        ref = lambda w: config['lcr-modules']['sigprofiler']["sigpro_genomes"][w.genome_build]
    conda: 
        CFG["conda_envs"]["sigprofiler"]
    threads: 
        CFG["threads"]["generator"]
    resources:
        mem_mb = CFG["mem_mb"]["generator"]
    shell:
        "python {input.script} {wildcards.sample_set} {params.ref} {input.maf} > {log.stdout} 2> {log.stderr}"

# Performs NMF on a wide range of ranks to obtain the optimal signature rank
rule _sigprofiler_run_estimate:
    input:
        unpack(get_dir)
    output:
        stat = CFG["dirs"]["estimate"]+"{seq_type}--{genome_build}/{sample_set}/{type}/All_solutions_stat.csv"
    log:
        stdout = CFG["logs"]["estimate"] + "{seq_type}--{genome_build}/{sample_set}/estimate.{type}.stdout.log",
        stderr = CFG["logs"]["estimate"] + "{seq_type}--{genome_build}/{sample_set}/estimate.{type}.stderr.log"
    params:
        opts = CFG["options"]["estimate"],
        ref = lambda w: config['lcr-modules']['sigprofiler']["sigpro_genomes"][w.genome_build],
        context_type = '96,DINUC,ID',
        exome = lambda w: {'genome': 'False', 'capture': 'True'}[w.seq_type],
        min_sig = 1,
        max_sig = lambda w: max_sigs[w.type],
        outpath = CFG["dirs"]["estimate"]+"{seq_type}--{genome_build}/{sample_set}"
    conda: CFG["conda_envs"]["sigprofiler"]
    threads: CFG["threads"]["estimate"]
    resources:
        mem_mb = CFG["mem_mb"]["estimate"]
    shell:
        op.as_one_line("""
        python {input.script} {params.opts} {input.mat} {params.outpath} 
        {params.ref} {params.context_type} {params.exome} 
        {params.min_sig} {params.max_sig} {threads} > {log.stdout} 2> {log.stderr}
        """)

# Runs NMF around the optimal rank from _sigprofiler_run_estimate
# using higher NMF replications to get final signatures
rule _sigprofiler_run_extract:
    input:
        unpack(get_dir),
        stat = str(rules._sigprofiler_run_estimate.output.stat)
    output:
        decomp = CFG["dirs"]["extract"]+"{seq_type}--{genome_build}/{sample_set}/{type}/Suggested_Solution/COSMIC_{type}_Decomposed_Solution/De_Novo_map_to_COSMIC_{type}.csv"
    log:
        stdout = CFG["logs"]["extract"] + "{seq_type}--{genome_build}/{sample_set}/extract.{type}.stdout.log",
        stderr = CFG["logs"]["extract"] + "{seq_type}--{genome_build}/{sample_set}/extract.{type}.stderr.log"
    params:
        opts = CFG["options"]["extract"],
        ref = lambda w: config['lcr-modules']['sigprofiler']["sigpro_genomes"][w.genome_build],
        context_type = '96,DINUC,ID',
        exome = lambda w: {'genome': 'False', 'capture': 'True'}[w.seq_type],
        outpath = CFG["dirs"]["extract"]+"{seq_type}--{genome_build}/{sample_set}"
    conda: CFG["conda_envs"]["sigprofiler"]
    threads: CFG["threads"]["extract"]
    resources:
        mem_mb = CFG["mem_mb"]["extract"]
    shell:
        op.as_one_line("""
        python {input.script} {params.opts} {input.mat} {params.outpath} 
        {params.ref} {params.context_type} {params.exome}
        $(awk 'BEGIN {{OFS=FS=","}} $1 ~ "*" {{S=substr($1,1,length($1)-1); if (S-2<1) {{print 1}} else {{print S-2}}}}' {input.stat})
        $(awk 'BEGIN {{OFS=FS=","}} $1 ~ "*" {{S=substr($1,1,length($1)-1); print S+2}}' {input.stat})
        {threads} > {log.stdout} 2> {log.stderr}
        """)

# Generate contents file
rule _sigprofiler_make_contents:
    input:
        decomp = str(rules._sigprofiler_run_extract.output.decomp)
    output:
        contents = CFG["dirs"]["outputs"] + "contents/{seq_type}--{genome_build}/{sample_set}.{genome_build}.{type}.contents"
    run:
        import pandas as pd

        # Import
        meta = pd.read_table(config['lcr-modules']['sigprofiler']["inputs"]["samples_metadata"])
        tbl = pd.read_table(config['lcr-modules']['sigprofiler']["inputs"]["sample_set_table"])
        samples = tbl.loc[tbl[wildcards.sample_set]==1, "sample_id"].tolist()
        sample_subset = meta[meta["sample_id"].isin(samples)]

        # Subset
        contents = sample_subset[["sample_id", "seq_type", "genome_build"]]
        contents = contents[contents["seq_type"] == wildcards.seq_type]

        # Mutate
        contents = contents.assign(normal_id = '.', pair_status = '.')
        contents = contents.assign(current_genome_build = wildcards.genome_build)
        contents = contents.assign(filename = '.')
        contents = contents.assign(is_lifted = ['native' if x == wildcards.genome_build else 'lifted' for x in contents['genome_build']])

        # Rename, reorder
        contents = contents.rename({'sample_id':'tumour_id'}, axis=1)
        contents = contents[["tumour_id","normal_id","pair_status","is_lifted","seq_type","genome_build","current_genome_build","filename"]]

        # Write contents file
        contents.to_csv(output.contents, sep = '\t', index = False)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sigprofiler_output_tsv:
    input:
        decomp = str(rules._sigprofiler_run_extract.output.decomp),
    output:
        decomp = CFG["dirs"]["outputs"] + "cosmic_sigs/{seq_type}--{genome_build}/{sample_set}.COSMIC_{type}.De_Novo_map_to_COSMIC_{type}.csv" 
    run:
        op.relative_symlink(input.decomp, output.decomp)

# Generates the target sentinels for each run, which generate the symlinks
rule _sigprofiler_all:
    input:
        expand(
            expand(
                [
                    CFG["dirs"]["outputs"] + "cosmic_sigs/{seq_type}--{genome_build}/{sample_set}.COSMIC_{{type}}.De_Novo_map_to_COSMIC_{{type}}.csv",
                    CFG["dirs"]["outputs"] + "contents/{seq_type}--{genome_build}/{sample_set}.{genome_build}.{{type}}.contents"
                ],
                zip,
                seq_type = CFG["mafs"]["seq_type"],
                genome_build = CFG["mafs"]["projected_build"],
                sample_set = CFG["mafs"]["sample_set"]),
            type = CFG["inputs"]["type"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
