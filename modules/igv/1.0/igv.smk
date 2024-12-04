#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd
from datetime import datetime
import os

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
# `CFG` is a shortcut to `config["lcr-modules"]["igv"]`
CFG = op.setup_module(
    name = "igv",
    version = "1.0",
    subdirectories = ["inputs", "batch_scripts", "igv", "snapshots", "outputs"],
)

# Rename genome_build values in sample metadata to correlate with MAF values
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["options"]["genome_map"]["grch37"]), "grch37", inplace=True)
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["options"]["genome_map"]["hg38"]), "hg38", inplace=True)
CFG["runs"]["normal_genome_build"].mask(CFG["runs"]["normal_genome_build"].isin(CFG["options"]["genome_map"]["grch37"]), "grch37", inplace=True)
CFG["runs"]["normal_genome_build"].mask(CFG["runs"]["normal_genome_build"].isin(CFG["options"]["genome_map"]["hg38"]), "hg38", inplace=True)

# Setup variables 
SUFFIX = ".pad" + str(CFG["options"]["generate_batch_script"]["padding"])
if "launch_date" in CFG:
    LAUNCH_DATE = CFG["launch_date"]
else:
    LAUNCH_DATE = datetime.today().strftime('%Y-%m-%d')

# Define rules to be run locally when using a compute cluster
localrules:
    _igv_symlink_bams,
    _igv_symlink_maf,
    _igv_reduce_maf_cols,
    _igv_merge_regions,
    _igv_format_regions,
    _igv_liftover_regions,
    _igv_merge_lifted_regions,
    _igv_filter_maf,
    _igv_create_batch_script_per_variant,
    _igv_batches_to_merge,
    _igv_download_igv,
    _igv_run,
    _igv_track_failed,
    _igv_quality_control,
    _igv_symlink_snapshot,
    _igv_check_snapshots,
    _igv_check_samples,
    _igv_setup_estimates,
    _igv_estimate_snapshots,
    _igv_check_estimates


##### FUNCTIONS #####


def get_bam(wildcards):
    CFG = config["lcr-modules"]["igv"]
    metadata = CFG["samples"]
    genome_build = metadata[(metadata.sample_id == wildcards.sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"]

    return expand(
        CFG["inputs"]["bam_path"],
        seq_type = wildcards.seq_type,
        sample_id = wildcards.sample_id,
        genome_build = genome_build
    )

def get_bai(wildcards):
    CFG = config["lcr-modules"]["igv"]
    metadata = CFG["samples"]
    genome_build = metadata[(metadata.sample_id == wildcards.sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"]

    return expand(
        CFG["inputs"]["bai_path"],
        seq_type = wildcards.seq_type,
        sample_id = wildcards.sample_id,
        genome_build = genome_build
    )

##### RULES #####



# Symlinks the input files into the module results directory (under '00-inputs/')
rule _igv_symlink_bams:
    input:
        bam = get_bam,
        bai = get_bai
    output:
        bam = CFG["dirs"]["inputs"] + "bams/{seq_type}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bams/{seq_type}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

rule _igv_symlink_maf:
    input:
        maf = CFG["inputs"]["maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Reduce MAF columns to prevent parsing errors in Pandas and reduce file size
rule _igv_reduce_maf_cols:
    input:
        maf = str(rules._igv_symlink_maf.output.maf)
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf.temp")
    shell:
        op.as_one_line("""
        cut -f 1,5,6,7,9,10,11,13,16,17 {input.maf} > {output.maf}
        """)

rule _igv_merge_regions:
    input:
        input_regions = lambda w: config["lcr-modules"]["igv"]["regions"][w.tool_type][w.tool_build]
    output:
        merged_regions = CFG["dirs"]["inputs"] + "regions/{tool_type}_merged.{tool_build}.tsv"
    log:
        stdout = CFG["logs"]["inputs"] + "regions/merge_{tool_type}_regions.{tool_build}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "regions/merge_{tool_type}_regions.{tool_build}.stderr.log"
    run:
        merged_df = pd.DataFrame()
        for result in input.input_regions:
            try:
                df = pd.read_table(result, comment="#", sep="\t")
                merged_df = pd.concat([merged_df, df])
            except:
                with open(log.stdout, "a") as header:
                    header.write(f"Error reading or merging file {result}\n")
        merged_df.to_csv(output.merged_regions, sep="\t", na_rep="NA", index=False)

# Convert input regions file into BED format
rule _igv_format_regions:
    input:
        regions = str(rules._igv_merge_regions.output.merged_regions)
    output:
        regions = CFG["dirs"]["inputs"] + "regions/{tool_type}_formatted.{tool_build}.tsv"
    params:
        regions_format = lambda w: w.tool_type,
        regions_build = lambda w: w.tool_build
    log:
        stdout = CFG["logs"]["inputs"] + "regions/format_regions_{tool_type}.{tool_build}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "regions/format_regions_{tool_type}.{tool_build}.stderr.log"
    script:
        config["lcr-modules"]["igv"]["scripts"]["format_regions"]

def _igv_get_chain(wildcards):
    if "38" in str({wildcards.tool_build}):
        return reference_files("genomes/{tool_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{tool_build}/chains/grch37/hg19ToHg38.over.chain")

rule _igv_liftover_regions:
    input:
        regions = str(rules._igv_format_regions.output.regions),
        igv_chain = _igv_get_chain
    output:
        regions = CFG["dirs"]["inputs"] + "regions/{tool_type}.{tool_build}To{genome_build}.lifted.txt"
    params:
        liftover_script = CFG["scripts"]["region_liftover_script"],
        liftover_minmatch = CFG["options"]["liftover_regions"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    resources:
        **CFG["resources"]["_igv_liftover_regions"]
    log:
        stderr = CFG["logs"]["inputs"] + "regions/liftover_regions_{tool_type}.{tool_build}To{genome_build}.stderr.log"
    shell:
        op.as_one_line("""
        echo "Running {rule} for {wildcards.tool_type} {wildcards.tool_build}" > {log.stderr} ;
        if [ {wildcards.tool_build} == {wildcards.genome_build} ] ; then
            echo "{wildcards.tool_type} build {wildcards.tool_build} is already target {wildcards.genome_build}. Nothing to be done, copying {input.regions} to {output.regions}..." >> {log.stderr} ;
            cat {input.regions} > {output.regions} ;
        else 
            echo "Lifting over {wildcards.tool_type} {wildcards.tool_build} regions to {wildcards.genome_build}..." >> {log.stderr} ;
            bash {params.liftover_script} 
            BED 
            {input.regions} 
            {output.regions} 
            {input.igv_chain} 
            YES 
            {params.liftover_minmatch}
            2>> {log.stderr} ;
        fi
        """)

def _get_lifted_regions(wildcards):
    CFG = config["lcr-modules"]["igv"]
    return expand(
        str(rules._igv_liftover_regions.output.regions),
        tool_type = list(CFG["regions"]),
        tool_build = ["grch37","hg38"],
        allow_missing = True
    )

rule _igv_merge_lifted_regions:
    input:
        regions = _get_lifted_regions
    output:
        regions = CFG["dirs"]["inputs"] + "regions/regions.{genome_build}.txt"
    run:
        merged_df = pd.DataFrame()
        for region in input.regions:
            try:
                df = pd.read_table(region, comment = "#", sep = "\t", header=None, usecols=[0,1,2])
                df.drop(df[df[0] == "chrom"].index, inplace = True)
                merged_df = pd.concat([merged_df, df])
            except:
                print(f"Lifted regions file is empty: {region}")
        merged_df = merged_df.drop_duplicates()
        merged_df.to_csv(output.regions, sep="\t", na_rep="NA", index=False, header=False)

# Filter MAF to lines containing positions of interest
rule _igv_filter_maf:
    input:
        maf = str(rules._igv_reduce_maf_cols.output.maf),
        regions = str(rules._igv_merge_lifted_regions.output.regions)
    output:
        maf = CFG["dirs"]["inputs"] + "maf/filtered_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf"
    params:
        regions_format = "bed",
        metadata = CFG["runs"]
    log:
        stdout = CFG["logs"]["inputs"] + "filter_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "filter_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.stderr.log"
    wildcard_constraints:
        seq_type = "[a-zA-Z]+"
    script:
        config["lcr-modules"]["igv"]["scripts"]["filter_script"]

def _get_maf(wildcards):
    CFG = config["lcr-modules"]["igv"]

    if wildcards.sample_id in list(CFG["runs"]["tumour_sample_id"]):
        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        return expand(
            str(rules._igv_filter_maf.output.maf),
            seq_type = wildcards.seq_type,
            genome_build = wildcards.genome_build,
            tumour_id = wildcards.sample_id,
            normal_sample_id = normal_sample_id,
            pair_status = pair_status
        )

    if wildcards.sample_id in list(CFG["runs"]["normal_sample_id"]):
        these_samples = op.filter_samples(CFG["runs"], normal_sample_id = wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build, pair_status = "matched")

        assert sorted(these_samples["tumour_genome_build"].drop_duplicates()) == sorted(these_samples["normal_genome_build"].drop_duplicates()), f"Different genome builds between normal ID and tumour ID for {wildcards.sample_id}"

        return expand(
            str(rules._igv_filter_maf.output.maf),
            zip,
            seq_type = these_samples["tumour_seq_type"],
            genome_build = these_samples["tumour_genome_build"],
            tumour_id = these_samples["tumour_sample_id"],
            normal_sample_id = these_samples["normal_sample_id"],
            pair_status = these_samples["pair_status"]
        )

checkpoint _igv_create_batch_script_per_variant:
    input:
        bam_file = str(rules._igv_symlink_bams.output.bam),
        bai_file = str(rules._igv_symlink_bams.output.bai),
        filter_maf = _get_maf
    output:
        finished = CFG["dirs"]["batch_scripts"] + "completed/{seq_type}--{genome_build}/{sample_id}.finished",
        variant_batch = expand(CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{{seq_type}}--{{genome_build}}/{preset}/{{sample_id}}" + SUFFIX + ".batch", preset = CFG["options"]["igv_presets"], allow_missing=True)
    params:
        tissue_status = lambda w: "normal" if w.sample_id in list(config["lcr-modules"]["igv"]["runs"]["normal_sample_id"]) else "tumour" if w.sample_id in list(config["lcr-modules"]["igv"]["runs"]["tumour_sample_id"]) else "unknown",
        batch_dir = CFG["dirs"]["batch_scripts"],
        snapshot_dir = CFG["dirs"]["snapshots"],
        genome_build = "{genome_build}",
        seq_type = "{seq_type}",
        batch_options = CFG["options"]["generate_batch_script"],
        suffix = SUFFIX,
        igv_presets = CFG["options"]["igv_presets"]
    log:
        stdout = CFG["logs"]["batch_scripts"] + "{seq_type}--{genome_build}/{sample_id}" + SUFFIX + ".stdout.log",
        stderr = CFG["logs"]["batch_scripts"] + "{seq_type}--{genome_build}/{sample_id}" + SUFFIX + ".stderr.log"
    script:
        config["lcr-modules"]["igv"]["scripts"]["batch_script_per_variant"]
 
rule _igv_batches_to_merge:
    input:
        batch_script = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{sample_id}" + SUFFIX + ".batch"
    output:
        dispatched_batch_script = CFG["dirs"]["batch_scripts"] + "dispatched_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{sample_id}" + SUFFIX + ".batch"
    params:
        merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{sample_id}" + SUFFIX + ".batch",
        igv_options = lambda w: config["lcr-modules"]["igv"]["options"]["generate_batch_script"]["igv_options"][w.preset]
    threads: (workflow.cores / 10)
    run:
        batch_script_path = os.path.abspath(input.batch_script)
        output_file = os.path.abspath(params.merged_batch)

        batch_script = open(batch_script_path, "r")

        with open(output_file, "r") as f:
            merged_lines = len(f.readlines())

        with open(output_file, "a") as handle:
            for line in batch_script:
                if merged_lines > 0:
                    if line.startswith(("load", "maxPanelHeight", "genome","setSleepInterval", "collapse")):
                        continue
                    if line.startswith(tuple(params.igv_options)) and not line.startswith("sort"):
                        continue
                handle.write(line)

        batch_script.close()

        output_touch = open(output.dispatched_batch_script, "w")
        output_touch.close()

def _evaluate_batches(wildcards):
    CFG = config["lcr-modules"]["igv"]
    checkpoint_output = checkpoints._igv_create_batch_script_per_variant.get(**wildcards).output.variant_batch

    if wildcards.sample_id in list(CFG["runs"]["tumour_sample_id"]):
        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf),
            seq_type = wildcards.seq_type,
            genome_build = wildcards.genome_build,
            tumour_id = wildcards.sample_id,
            normal_sample_id = normal_sample_id,
            pair_status = pair_status
        )

        maf_table = pd.read_table(maf[0], comment = "#", sep="\t")

    if wildcards.sample_id in list(CFG["runs"]["normal_sample_id"]):
        these_samples = op.filter_samples(CFG["runs"], normal_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build, pair_status="matched") # only get MAFs from matched tumour + normal combos

        mafs = expand(
                str(rules._igv_filter_maf.output.maf),
                zip,
                seq_type = these_samples["tumour_seq_type"],
                genome_build = these_samples["tumour_genome_build"],
                tumour_id = these_samples["tumour_sample_id"],
                normal_sample_id = these_samples["normal_sample_id"],
                pair_status = these_samples["pair_status"]
            )

        maf_table = pd.concat([pd.read_table(m, comment="#", sep="\t") for m in mafs])

    return expand(
            expand(
                str(rules._igv_batches_to_merge.output.dispatched_batch_script),
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                allow_missing = True
            ),
            preset = wildcards.preset,
            allow_missing = True
        )

IGV_VERSION = CFG["options"]["igv_download_path"]
IGV_VERSION = IGV_VERSION.split("/")[-1].replace(".zip","")

rule _igv_download_igv:
    output:
        igv_zip = CFG["dirs"]["igv"] + IGV_VERSION + ".zip",
        igv_installed = CFG["dirs"]["igv"] + IGV_VERSION + ".installed"
    params:
        igv = CFG["options"]["igv_download_path"]
    conda:
        CFG["conda_envs"]["wget"]
    log:
        stdout = CFG["logs"]["igv"] + "download/igv_" + IGV_VERSION + "_download.stdout.log",
        stderr = CFG["logs"]["igv"] + "download/igv_" + IGV_VERSION + "_download.stderr.log"
    shell:
        op.as_one_line("""
        wget -O {output.igv_zip} {params.igv} &&
        unzip {output.igv_zip} -d $(dirname {output.igv_zip}) > {log.stdout} 2> {log.stderr} &&
        touch {output.igv_installed}
        """)

checkpoint _igv_run:
    input:
        igv = ancient(str(rules._igv_download_igv.output.igv_installed)),
        finished_batches = str(rules._igv_create_batch_script_per_variant.output.finished),
        batch_script = _evaluate_batches
    output:
        complete = CFG["dirs"]["snapshots"] + "completed/{seq_type}--{genome_build}/{preset}/{sample_id}.completed"
    params:
        merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{sample_id}" + SUFFIX + ".batch",
        igv = CFG["dirs"]["igv"] + IGV_VERSION + "/igv.sh",
        sleep_time = CFG["options"]["generate_batch_script"]["sleep_timer"],
        server_number = "-n " + CFG["options"]["xvfb_parameters"]["server_number"] if CFG["options"]["xvfb_parameters"]["server_number"] is not None else "--auto-servernum",
        server_args = CFG["options"]["xvfb_parameters"]["server_args"]
    resources:
        **CFG["resources"]["_igv_run"]
    threads: (workflow.cores)
    log:
        stdout = CFG["logs"]["igv"] + "igv_run/{seq_type}--{genome_build}/{preset}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["igv"] + "igv_run/{seq_type}--{genome_build}/{preset}/{sample_id}.stderr.log"
    shell:
        op.as_one_line("""
        lines=$(wc -l < {params.merged_batch}) ;
        if [ $lines -gt 0 ] ;
        then
        if ! grep -q -e "exit" {params.merged_batch} ;
        then
        echo 'exit' >> {params.merged_batch} ;
        fi ; 
        maxtime=$(($(wc -l < {params.merged_batch}) * 10 + 15 + {params.sleep_time})) ;
        timeout $maxtime xvfb-run -s "-screen 0 1920x1080x24" {params.server_number} {params.server_args} {params.igv} -b {params.merged_batch} > {log.stdout} 2> {log.stderr} && touch {output.complete} ;
        exit=$? ;
        if [ $exit -ne 0 ] ;
        then 
        if grep -q -e "No such process" {log.stderr} && grep -q -e "Executing Command: exit" {log.stdout} ;
        then 
        echo "All IGV batch script commands have completed succesfully, but an Xvfb-run kill error has occurred." >> {log.stdout} && touch {output.complete} ;
        else 
        false ;
        fi ;
        fi ;
        else 
        echo 'Skipping sample {wildcards.sample_id} because it either has no variants to snapshot or all variants have already been snapshot.' >> {log.stdout} ;
        touch {output.complete} ;
        fi
        """)

rule _igv_track_failed:
    output:
        failed_summary = CFG["dirs"]["outputs"] + "snapshot_summaries/quality_control/qc_summary_" + LAUNCH_DATE + ".txt"
    run:
        header = "\t".join(["sample_id","seq_type","genome_build","gene","chromosome","position","preset","status","snapshot_path"])
        with open(output.failed_summary, "w") as handle:
            handle.write(header + "\n")

rule _igv_quality_control:
    input:
        igv = str(rules._igv_run.output.complete),
        failed_summary = ancient(str(rules._igv_track_failed.output.failed_summary)),
        snapshot = CFG["dirs"]["snapshots"] + "{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".png"
    output:
        snapshot_qc = CFG["dirs"]["snapshots"] + "qc/{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".qc"
    params:
        batch_script = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{sample_id}" + SUFFIX + ".batch",
        merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{sample_id}" + SUFFIX + ".batch",
        igv = CFG["dirs"]["igv"] + "IGV_Linux_2.7.2/igv.sh",
        server_number = "-n " + CFG["options"]["xvfb_parameters"]["server_number"] if CFG["options"]["xvfb_parameters"]["server_number"] is not None else "--auto-servernum",
        server_args = CFG["options"]["xvfb_parameters"]["server_args"],
        batch_temp = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{sample_id}" +  SUFFIX + ".batch.temp",
        thresholds = CFG["options"]["quality_control"]
    resources:
        **CFG["resources"]["_igv_quality_control"]
    threads: (workflow.cores)
    log:
        stdout = CFG["logs"]["snapshots"] + "quality_control/{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".stdout.log",
        stderr = CFG["logs"]["snapshots"] + "quality_control/{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".stderr.log"
    script:
        config["lcr-modules"]["igv"]["scripts"]["quality_control"]

rule _igv_symlink_snapshot:
    input:
        snapshot = CFG["dirs"]["snapshots"] + "{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".png",
        snapshot_qc = str(rules._igv_quality_control.output.snapshot_qc)
    output:
        snapshot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tissue_status}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{ref_allele}_{alt_allele}--{sample_id}" + SUFFIX + ".png"
    resources:
        **CFG["resources"]["_igv_symlink"]
    run:
        op.relative_symlink(input.snapshot, output.snapshot)

def _symlink_snapshot(wildcards):
    CFG = config["lcr-modules"]["igv"]
    checkpoint_outputs = checkpoints._igv_run.get(**wildcards).output.complete

    if wildcards.sample_id in list(CFG["runs"]["tumour_sample_id"]):
        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf),
            seq_type = wildcards.seq_type,
            genome_build = wildcards.genome_build,
            tumour_id = wildcards.sample_id,
            normal_sample_id = normal_sample_id,
            pair_status = pair_status
        )

        maf_table = pd.read_table(maf[0], comment = "#", sep="\t")

        tumour_snaps = expand(
            expand(
                str(rules._igv_symlink_snapshot.output.snapshot),
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                ref_allele = maf_table["Reference_Allele"],
                alt_allele = maf_table["Tumor_Seq_Allele2"],
                sample_id = maf_table["Tumor_Sample_Barcode"],
                allow_missing = True
            ),
            genome_build = wildcards.genome_build,
            seq_type = wildcards.seq_type,
            preset = wildcards.preset,
            tissue_status = "tumour"
        )

        return tumour_snaps

    if wildcards.sample_id in list(CFG["runs"]["normal_sample_id"]):
        these_samples = op.filter_samples(CFG["runs"], normal_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build, pair_status="matched")

        mafs = expand(
                str(rules._igv_filter_maf.output.maf),
                zip,
                seq_type = these_samples["tumour_seq_type"],
                genome_build = these_samples["tumour_genome_build"],
                tumour_id = these_samples["tumour_sample_id"],
                normal_sample_id = these_samples["normal_sample_id"],
                pair_status = these_samples["pair_status"]
            )

        maf_table = pd.concat([pd.read_table(m, comment="#", sep="\t") for m in mafs])

        normal_snaps = expand(
            expand(
                str(rules._igv_symlink_snapshot.output.snapshot),
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                ref_allele = maf_table["Reference_Allele"],
                alt_allele = maf_table["Reference_Allele"],
                sample_id = maf_table["Matched_Norm_Sample_Barcode"],
                allow_missing = True
            ),
            genome_build = wildcards.genome_build,
            seq_type = wildcards.seq_type,
            preset = wildcards.preset,
            tissue_status = "normal"
        )

        return normal_snaps

def _quality_control(wildcards):
    CFG = config["lcr-modules"]["igv"]
    checkpoint_outputs = checkpoints._igv_run.get(**wildcards).output.complete

    if wildcards.sample_id in list(CFG["runs"]["tumour_sample_id"]):
        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf),
            seq_type = wildcards.seq_type,
            genome_build = wildcards.genome_build,
            tumour_id = wildcards.sample_id,
            normal_sample_id = normal_sample_id,
            pair_status = pair_status
        )

        maf_table = pd.read_table(maf[0], comment = "#", sep="\t")

        tumour_snaps = expand(
            expand(
                str(rules._igv_quality_control.output.snapshot_qc),
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                ref_allele = maf_table["Reference_Allele"],
                alt_allele = maf_table["Tumor_Seq_Allele2"],
                sample_id = maf_table["Tumor_Sample_Barcode"],
                allow_missing = True
            ),
            genome_build = wildcards.genome_build,
            seq_type = wildcards.seq_type,
            preset = wildcards.preset,
            tissue_status = "tumour"
        )

        return tumour_snaps

    if wildcards.sample_id in list(CFG["runs"]["normal_sample_id"]):
        these_samples = op.filter_samples(CFG["runs"], normal_sample_id=wildcards.sample_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build, pair_status="matched") 

        mafs = expand(
                str(rules._igv_filter_maf.output.maf),
                zip,
                seq_type = these_samples["tumour_seq_type"],
                genome_build = these_samples["tumour_genome_build"],
                tumour_id = these_samples["tumour_sample_id"],
                normal_sample_id = these_samples["normal_sample_id"],
                pair_status = these_samples["pair_status"]
            )

        maf_table = pd.concat([pd.read_table(m, comment="#", sep="\t") for m in mafs])

        normal_snaps = expand(
            expand(
                str(rules._igv_quality_control.output.snapshot_qc), 
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                ref_allele = maf_table["Reference_Allele"],
                alt_allele = maf_table["Reference_Allele"],
                sample_id = maf_table["Matched_Norm_Sample_Barcode"],
                allow_missing = True
            ),
            genome_build = wildcards.genome_build,
            seq_type = wildcards.seq_type,
            preset = wildcards.preset,
            tissue_status = "normal"
        )

        return normal_snaps

rule _igv_check_snapshots:
    input:
        snapshots = _symlink_snapshot,
        quality_control = _quality_control
    output:
        complete = CFG["dirs"]["outputs"] + "completed/check_snapshots/{preset}/{seq_type}--{genome_build}/{sample_id}.checked"
    shell:
        "touch {output.complete}"

def _get_finished_samples(wildcards):
    if wildcards.pair_status == "matched":
        tumour_complete = ancient(
            expand(
                str(rules._igv_check_snapshots.output.complete), 
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.tumour_id,
                preset = wildcards.preset
            )
        )

        normal_complete = ancient(
            expand(
                str(rules._igv_check_snapshots.output.complete),
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.normal_sample_id,
                preset = wildcards.preset
            )
        )

        return (tumour_complete + normal_complete)

    else:
        return ancient(
            expand(
                str(rules._igv_check_snapshots.output.complete),
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.tumour_id,
                preset = wildcards.preset
            )
        )

rule _igv_check_samples:
    input:
        igv_completed = _get_finished_samples,
    output:
        checked = CFG["dirs"]["outputs"] + "completed/check_samples/{preset}/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.completed"
    shell:
        "touch {output.checked}"

##### Rules below will only run if CFG["estimate_only"] set to True

rule _igv_setup_estimates:
    output:
        snapshot_summary = CFG["dirs"]["outputs"] + "snapshot_summaries/estimates/snapshot_summary_" + LAUNCH_DATE + ".txt",
        snapshot_estimate = CFG["dirs"]["outputs"] + "snapshot_summaries/estimates/snapshot_estimate_" + LAUNCH_DATE + ".txt"
    run:
        header = "\t".join(["sample_id","seq_type","genome_build","gene","chromosome","position", "igv_preset"])
        with open(output.snapshot_summary,"w") as handle:
            handle.write(header + "\n")
        with open(output.snapshot_estimate, "w") as handle:
            handle.write(header + "\n")

rule _igv_estimate_snapshots:
    input:
        maf = _get_maf,
        snapshot_summary = str(rules._igv_setup_estimates.output.snapshot_summary),
        snapshot_estimate = str(rules._igv_setup_estimates.output.snapshot_estimate)
    output:
        complete = temp(CFG["dirs"]["batch_scripts"] + "estimate_batch_scripts/{seq_type}--{genome_build}/{preset}/{sample_id}.temp")
    threads: (workflow.cores)
    run:
        CFG = config["lcr-modules"]["igv"]
        snapshot_summary = open(input.snapshot_summary, "a")
        snapshot_estimate = open(input.snapshot_estimate, "a")

        maf_table = pd.concat([pd.read_table(file, sep="\t", comment="#") for file in input.maf])

        seq_type = wildcards.seq_type
        genome_build = wildcards.genome_build
        sample_id = wildcards.sample_id
        preset = wildcards.preset

        for index, row in maf_table.iterrows():
            gene = row["Hugo_Symbol"]
            chromosome = row["chr_std"]
            start_position = str(row["Start_Position"])

            dispatch_path = CFG["dirs"]["batch_scripts"] + f"dispatched_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{sample_id}" + SUFFIX + ".batch"
        
            outline = "\t".join([sample_id, seq_type, genome_build, gene, chromosome, start_position, preset])
            snapshot_summary.write(outline + "\n")

            if not os.path.exists(dispatch_path):
                snapshot_estimate.write(outline + "\n")

        finished = open(output.complete, "w")
        finished.close()
        snapshot_summary.close()
        snapshot_estimate.close()

def _check_estimates(wildcards):
    if wildcards.pair_status == "matched":
        tumour_complete = ancient(
            expand(
                str(rules._igv_estimate_snapshots.output.complete), 
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.tumour_id,
                preset = wildcards.preset
            )
        )

        normal_complete = ancient(
            expand(
                str(rules._igv_estimate_snapshots.output.complete),
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.normal_sample_id,
                preset = wildcards.preset
            )
        )

        return (tumour_complete + normal_complete)

    else:
        return ancient(
            expand(
                str(rules._igv_estimate_snapshots.output.complete),
                seq_type = wildcards.seq_type,
                genome_build = wildcards.genome_build,
                sample_id = wildcards.tumour_id,
                preset = wildcards.preset
            )
        )

rule _igv_check_estimates:
    input:
        estimates_completed = _check_estimates
    output:
        sample_estimated = temp(CFG["dirs"]["outputs"] + "snapshot_summaries/temp/{seq_type}--{genome_build}/{preset}/{tumour_id}--{normal_sample_id}--{pair_status}.temp")
    shell:
        "touch {output.sample_estimated}"

if CFG["estimate_only"] is False:
    rule _igv_all:
        input:
            expand(
                expand(
                    [
                        str(rules._igv_check_samples.output.checked)
                    ], 
                    zip, 
                    tumour_id=CFG["runs"]["tumour_sample_id"], 
                    normal_sample_id=CFG["runs"]["normal_sample_id"],
                    seq_type=CFG["runs"]["tumour_seq_type"], 
                    genome_build=CFG["runs"]["tumour_genome_build"],
                    pair_status = CFG["runs"]["pair_status"],
                    allow_missing=True
                ),
                preset=CFG["options"]["igv_presets"]
            )


if CFG["estimate_only"] is True:
    rule _igv_all:
        input:
            expand(
                expand(
                    str(rules._igv_check_estimates.output.sample_estimated),
                    zip,
                    tumour_id=CFG["runs"]["tumour_sample_id"],
                    normal_sample_id=CFG["runs"]["normal_sample_id"],
                    seq_type = CFG["runs"]["tumour_seq_type"],
                    genome_build = CFG["runs"]["tumour_genome_build"],
                    pair_status=CFG["runs"]["pair_status"],
                    allow_missing=True
                ),
                preset=CFG["options"]["igv_presets"]
            )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
