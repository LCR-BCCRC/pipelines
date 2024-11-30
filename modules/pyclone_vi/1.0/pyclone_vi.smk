#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Andrew Roth
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import hashlib
import glob

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["pyclone_vi"]`
CFG = op.setup_module(
    name = "pyclone_vi",
    version = "1.0",
    subdirectories = ["inputs", "build_inputs", "fit", "og-pyclone", "phyclone", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _pyclone_vi_write_results,
    _pyclone_vi_all

# Install GAMBLR

# Obtain the path to the GAMBLR conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open("config/envs/GAMBLR.yaml", 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
GAMBLR = glob.glob(conda_prefix + "/" + h[:8] + "*")
for file in GAMBLR:
    if os.path.isdir(file):
        GAMBLR = file

rule _pyclone_vi_install_GAMBLR:
    params:
        branch = ", ref = \"" + CFG["options"]["build_input"]['gamblr_branch'] + "\"" if CFG["options"]["build_input"]['gamblr_branch'] != "" else "",
        config_url = CFG["options"]["build_input"]["gamblr_config_url"]
    output:
        installed = directory(GAMBLR + "/lib/R/library/GAMBLR"),
        config = "gamblr.yaml"
    conda:
        CFG['conda_envs']['gamblr']
    shell:
        op.as_one_line("""
        wget -qO {output.config} {params.config_url} &&
        R -q -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR"{params.branch})'
        """)


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _pyclone_vi_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.maf"
    group: "input_and_build"
    run:
        op.relative_symlink(input.maf, output.maf)

rule _pyclone_vi_input_battenberg:
    input:
        subclones = CFG["inputs"]["sample_subclones"],
        cellularity = CFG["inputs"]["sample_cellularity"],
        sex = CFG["inputs"]["sample_sex"]
    output:
        subclones = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.battenberg.subclones.txt",
        cellularity = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.battenberg.cellularity_ploidy.txt",
        sex = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.battenberg.inferred_sex.txt"
    group: "input_and_build"
    wildcard_constraints:
        seq_type = "genome"
    run:
        op.absolute_symlink(input.subclones, output.subclones)
        op.absolute_symlink(input.cellularity, output.cellularity)
        op.absolute_symlink(input.sex, output.sex)


# Prepare Pyclone inputs


rule _pyclone_vi_subset_maf:
    input:
        maf = str(rules._pyclone_vi_input_maf.output.maf),
        GAMBLR = ancient(rules._pyclone_vi_install_GAMBLR.output.installed)
    output:
        maf = CFG["dirs"]["build_inputs"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.subset.maf"
    params:
        script = CFG["scripts"]["subset_maf"]
    conda:
        CFG["conda_envs"]["gamblr"]
    script:
        "{params.script}"

subset_maf = CFG["options"]["build_input"]["subset_maf"]

rule  _pyclone_vi_build_input:
    input:
        maf = str(rules._pyclone_vi_subset_maf.output.maf) if subset_maf else str(rules._pyclone_vi_input_maf.output.maf),
        cnv = expand(str(rules._pyclone_vi_input_battenberg.output.subclones), seq_type = "genome", allow_missing = True),
        cellularity = expand(str(rules._pyclone_vi_input_battenberg.output.cellularity), seq_type = "genome", allow_missing = True),
        sex = expand(str(rules._pyclone_vi_input_battenberg.output.sex), seq_type = "genome", allow_missing = True)
    output:
        inputs = CFG["dirs"]["build_inputs"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.inputs.tsv"
    params:
        script = CFG["scripts"]["build_input"],
        version = lambda w: {"genome": "pyclone-vi", "capture": "pyclone"}[w.seq_type]
    log:
        CFG["logs"]["build_inputs"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}/build_input.log"
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["build_input"]
    resources:
        **CFG["resources"]["build_input"]
    group: "input_and_build"
    shell:
        op.as_one_line("""
        cellularity=$(tail -n +2 {input.cellularity} | cut -f 1);
        if [[ $(tail -n+2 {input.sex} | cut -f4) == "female" ]]; then sex="F"; else sex="M"; fi;
        echo "Prepping PyClone-vi inputs with cellularity $cellularity and sex $sex. ";
        python {params.script} -c battenberg -s $sex -t $cellularity -ic {input.cnv} -is {input.maf} -o {output.inputs} -id {wildcards.tumour_id} -p {params.version}
        > {log} 2>&1
        """)

# Merge all built inputs into a single tsv file per patient
def get_built_inputs(wildcards):
    CFG = config["lcr-modules"]["pyclone_vi"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    inputs = expand(
        [
            str(rules._pyclone_vi_build_input.output.inputs)
        ],
        zip,
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"],
        allow_missing = True
    )
    return(inputs)

def get_cellularity(wildcards):
    CFG = config["lcr-modules"]["pyclone_vi"]
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id).sort_values(by = ["tumour_time_point"])
    inputs = expand(
        rules._pyclone_vi_input_battenberg.input.cellularity,
        zip,
        seq_type = ["genome"]*len(PATIENT["tumour_genome_build"]),
        genome_build = PATIENT["tumour_genome_build"],
        tumour_id = PATIENT["tumour_sample_id"],
        normal_id = PATIENT["normal_sample_id"],
        pair_status = PATIENT["pair_status"]
    )
    cell_list = []
    for file in inputs:
        cellularity = pd.read_csv(file, sep = "\t")
        cellularity = cellularity["cellularity"].tolist()
        cell_list = cell_list + cellularity
    cell_list = [str(i) for i in cell_list]
    cellularities = " ".join(cell_list)
    sample_ids = " ".join(PATIENT.tumour_sample_id.tolist())
    return dict(cellularity = cellularities, sample_ids = sample_ids)

rule _pyclone_vi_merge_input:
    input:
        get_built_inputs
    output:
        merged = CFG["dirs"]["build_inputs"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.merged_inputs.tsv"
    group: "input_and_build"
    run:
        df_list = []
        for tsv in input:
            df = pd.read_csv(tsv, sep = "\t", header = 0, index_col = None)
            df_list.append(df)
        df_merged = pd.concat(df_list, ignore_index = True, axis = 0)
        df_merged.sort_values("mutation_id", inplace=True)
        df_merged.rename(columns = {"var_counts": "alt_counts"}, inplace=True)
        df_merged.to_csv(output.merged, sep = "\t", index=False)

# PyClone for capture data
rule _pyclone_run_analysis_pipeline:
    input:
        tsv = get_built_inputs
    output:
       loci = CFG["dirs"]["og-pyclone"] + "{seq_type}--{genome_build}/{patient_id}/tables/loci.tsv"
    params:
        workdir = CFG["dirs"]["og-pyclone"] + "{seq_type}--{genome_build}/{patient_id}",
        cellularity = lambda w: get_cellularity(w)["cellularity"],
        sample_ids = lambda w: get_cellularity(w)["sample_ids"]
    conda: CFG["conda_envs"]["pyclone"]
    wildcard_constraints: seq_type = "capture"
    shell:
        op.as_one_line("""
        PyClone run_analysis_pipeline
            --in_files {input.tsv}
            --working_dir {params.workdir}
            --tumour_contents {params.cellularity}
            --samples {params.sample_ids}
            --plot_file_format pdf
        """)


# Fit pyclone-vi
rule _pyclone_vi_fit:
    input:
        tsv = str(rules._pyclone_vi_merge_input.output.merged)
    output:
        trace = CFG["dirs"]["fit"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.h5"
    log:
        CFG["logs"]["fit"] + "{seq_type}--{genome_build}/{patient_id}/fit.log"
    params:
        **CFG["options"]["fit"]
    threads:
        CFG["threads"]["fit"]
    resources:
        **CFG["resources"]["fit"]
    conda:
        CFG["conda_envs"]["pyclone-vi"]
    group: "fit_and_write"
    wildcard_constraints: seq_type = "genome"
    shell:
        op.as_one_line("""
        pyclone-vi fit
        -i {input.tsv}
        -o {output.trace}
        -c {params.num_clusters}
        -d binomial
        -r {params.num_restarts}
        {params.opts}
        > {log} 2>&1
        """)

# Fit pyclone-vi
rule _pyclone_vi_write_results:
    input:
        trace = str(rules._pyclone_vi_fit.output.trace)
    output:
        results = CFG["dirs"]["fit"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.pyclone_vi.tsv"
    log:
        CFG["logs"]["fit"] + "{seq_type}--{genome_build}/{patient_id}/write.log"
    conda:
        CFG["conda_envs"]["pyclone-vi"]
    group: "fit_and_write"
    wildcard_constraints: seq_type = "genome"
    shell:
        op.as_one_line("""
        pyclone-vi write-results-file -i {input.trace} -o {output.results} > {log} 2>&1
        """)

# Run PhyClone to get phylogenetic tree
rule _pyclone_vi_run_phyclone:
    input:
        merged = str(rules._pyclone_vi_merge_input.output.merged),
        pyclone = lambda w: str(rules._pyclone_vi_write_results.output.results) if w.seq_type == "genome" else str(rules._pyclone_run_analysis_pipeline.output.loci)
    output:
        trace = CFG["dirs"]["phyclone"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.phyclone.pkl.gz"
    params:
        **CFG["options"]["phyclone"]
    threads:
        CFG["threads"]["phyclone"]
    resources:
        **CFG["resources"]["phyclone"]
    conda:
        CFG["conda_envs"]["phyclone"]
    shell:
        op.as_one_line("""
        phyclone run
        -i {input.merged}
        -c {input.pyclone}
        -o {output.trace}
        -b {params.burnin}
        -d {params.density}
        -n {params.num_iters}
        """)

rule _pyclone_vi_phyclone_consensus:
    input:
        merged = str(rules._pyclone_vi_merge_input.output.merged),
        pyclone = str(rules._pyclone_vi_run_phyclone.output.trace)
    output:
        tree = CFG["dirs"]["phyclone"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.phyclone.tree.nwk",
        clusters = CFG["dirs"]["phyclone"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.phyclone.clusters.tsv"
    params:
        **CFG["options"]["phyclone"]
    threads:
        CFG["threads"]["phyclone"]
    resources:
        **CFG["resources"]["phyclone"]
    conda:
        CFG["conda_envs"]["phyclone"]
    shell:
        op.as_one_line("""
        phyclone consensus
        -i {input.pyclone}
        -t {output.tree}
        -o {output.clusters}
        """)

rule _pyclone_vi_compute_tree_stats:
    input:
        trace = str(rules._pyclone_vi_run_phyclone.output.trace)
    output:
        stats = CFG["dirs"]["phyclone"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.phyclone.stats.tsv"
    params:
        **CFG["options"]["phyclone"],
        script = CFG["scripts"]["compute_stats"]
    threads:
        CFG["threads"]["phyclone"]
    resources:
        **CFG["resources"]["phyclone"]
    conda:
        CFG["conda_envs"]["phyclone"]
    shell:
        op.as_one_line("""
        python {params.script}
        -i {input.trace}
        -o {output.stats}
        -p {wildcards.patient_id}
        -b {params.burnin}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _pyclone_vi_output_tsv:
    input:
        pyclone = lambda w: str(rules._pyclone_vi_write_results.output.results) if w.seq_type == "genome" else str(rules._pyclone_run_analysis_pipeline.output.loci),
        phyclone = str(rules._pyclone_vi_compute_tree_stats.output.stats),
        tree = str(rules._pyclone_vi_phyclone_consensus.output.tree),
        clusters = str(rules._pyclone_vi_phyclone_consensus.output.clusters),
    output:
        pyclone = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{patient_id}.pyclone.results.tsv",
        phyclone = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{patient_id}.phyclone.stats.tsv",
        tree = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{patient_id}.phyclone.tree.nwk",
        clusters = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{patient_id}.phyclone.clusters.tsv"
    run:
        op.relative_symlink(input.pyclone, output.pyclone, in_module = TRUE)
        op.relative_symlink(input.phyclone, output.phyclone)
        op.relative_symlink(input.tree, output.tree)
        op.relative_symlink(input.clusters, output.clusters)


# Generates the target sentinels for each run, which generate the symlinks

# Generate a de-duplicated table of patient_ids etc.
PATIENTS_GENOMES = op.filter_samples(CFG["runs"][["tumour_patient_id", "normal_patient_id", "tumour_genome_build", "tumour_seq_type"]], tumour_seq_type = "genome")\
    .drop_duplicates(subset = None, ignore_index = True)
PATIENTS_CAPTURE = op.filter_samples(CFG["runs"][["tumour_patient_id", "normal_patient_id", "tumour_genome_build", "tumour_seq_type"]], tumour_seq_type = "capture")\
    .drop_duplicates(subset = None, ignore_index = True)
if isinstance(PATIENTS_GENOMES, pd.DataFrame) and isinstance(PATIENTS_CAPTURE, pd.DataFrame):
    PATIENTS = pd.concat([PATIENTS_GENOMES, PATIENTS_CAPTURE])

rule _pyclone_vi_all:
    input:
        expand(
            rules._pyclone_vi_input_maf.output.maf,
            zip,
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            patient_id = CFG["runs"]["tumour_patient_id"]
        ),
        expand(
            [
                str(rules._pyclone_vi_output_tsv.output.phyclone),
                str(rules._pyclone_vi_output_tsv.output.pyclone),
                str(rules._pyclone_vi_output_tsv.output.tree),
                str(rules._pyclone_vi_output_tsv.output.clusters)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=PATIENTS["tumour_seq_type"],
            genome_build=PATIENTS["tumour_genome_build"],
            patient_id=PATIENTS["tumour_patient_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
