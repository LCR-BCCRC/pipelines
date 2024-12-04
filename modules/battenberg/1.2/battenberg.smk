#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob
import pandas as pd

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
    version = "1.2",
    subdirectories = ["inputs", "infer_sex", "battenberg", "convert_coordinates", "fill_regions", "normalize", "outputs"],
)

#set variable for prepending to PATH based on config
BATTENBERG_SCRIPT_PATH = CFG['inputs']['src_dir']
#this is used in place of the shell.prefix() because that was not working consistently. This is not ideal. 

#this preserves the variable when using lambda functions
#_battenberg_CFG = CFG

# Define rules to be run locally when using a compute cluster
localrules:
    _battenberg_all

BATTENBERG_VERSION_MAP = {
    "hg19": "grch37",
    "grch37": "grch37",
    "hs37d5": "grch37",
    "hg38": "hg38",
    "grch38": "hg38",
    "grch38-legacy": "hg38"

}

##### RULES #####

# Downloads the reference files into the module results directory (under '00-inputs/') from https://www.bcgsc.ca/downloads/morinlab/reference/ . 
rule _battenberg_get_reference:
    output:
        battenberg_impute =  directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_impute_v3"),
        impute_info = CFG["dirs"]["inputs"] + "reference/{genome_build}/impute_info.txt",
        probloci =  CFG["dirs"]["inputs"] + "reference/{genome_build}/probloci.txt.gz",
        battenberg_wgs_replic_correction = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_replic_correction_1000g_v3"),
        battenberg_gc_correction = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_gc_correction_1000g_v3"),  
        genomesloci = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_1000genomesloci2012_v3")
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/reference",
        alt_build = lambda w: BATTENBERG_VERSION_MAP[w.genome_build],
        folder = CFG["dirs"]["inputs"] + "reference/{genome_build}",
        build = "{genome_build}",
        battenberg_path = CFG['inputs']['src_dir']
    resources:
        **CFG["resources"]["reference"]
    threads:
        CFG["threads"]["reference"]
    conda: CFG["conda"]["wget"]
    shell:
        op.as_one_line("""
        wget -qO-  {params.url}/battenberg_impute_{params.alt_build}.tar.gz  |
        tar -xvz > {output.battenberg_impute} -C {params.folder}
        &&
        wget -qO- {params.url}/battenberg_{params.alt_build}_gc_correction.tar.gz  |
        tar -xvz > {output.battenberg_gc_correction} -C {params.folder}
        &&
        wget -qO- {params.url}/battenberg_1000genomesloci_{params.alt_build}.tar.gz | 
        tar -xvz > {output.genomesloci} -C {params.folder}
        &&
        wget -O {output.impute_info} {params.url}/impute_info_{params.alt_build}.txt
        &&
        python {params.battenberg_path}/reference_correction.py {params.build} $(dirname $(readlink -f {output.impute_info}))
        &&
        wget -qO-  {params.url}/battenberg_{params.alt_build}_replic_correction.tar.gz |
        tar -xvz > {output.battenberg_wgs_replic_correction} -C {params.folder}
        &&
        wget -O {output.probloci} {params.url}/probloci_{params.alt_build}.txt.gz
        
        """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    group: "setup_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam + ".bai", output.bai)
        op.absolute_symlink(input.bam + ".bai", output.crai)

# this process is very fast on bam files and painfully slow on cram files. 
# The result of calc_sex_status.sh is stored in a file to avoid having to rerun it unnecessarily
rule _infer_patient_sex:
    input: 
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex"
    resources:
        **CFG["resources"]["infer_sex"]
    log:
        stderr = CFG["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    group: "setup_run"
    threads: 8
    shell:
        op.as_one_line(""" 
        PATH={BATTENBERG_SCRIPT_PATH}:$PATH;
        echo "running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr} ;
        calc_sex_status.sh {input.normal_bam} {input.fasta} {wildcards.normal_id} > {output.sex_result} 2>> {log.stderr} &&
        echo "DONE running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" >> {log.stderr} 
        """)


# This rule runs the entire Battenberg pipeline. Eventually we may want to set this rule up to allow re-starting
# of partially completed jobs (e.g. if they run out of RAM and are killed by the cluster, they can automatically retry)
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        impute_info = str(rules._battenberg_get_reference.output.impute_info)

    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        ac=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"),
        mb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab"),
        mlrg=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab"),
        mlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab"),
        nlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab"),
        nb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab"),
        cp=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cellularity_ploidy.txt"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        script = CFG["inputs"]["battenberg_script"],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        ref = CFG["dirs"]["inputs"] + "reference/{genome_build}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    shell:
       op.as_one_line("""
        if [[ $(head -c 4 {params.fasta}) == ">chr" ]]; then chr_prefixed='true'; else chr_prefixed='false'; fi;
        echo "$chr_prefixed"
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stdout};
        sex=$(cut -f 4 {input.sex_result}| tail -n 1); 
        echo "setting sex as $sex";
        Rscript --vanilla {params.script} -t {wildcards.tumour_id} 
        -n {wildcards.normal_id} --tb $(readlink -f {input.tumour_bam}) --nb $(readlink -f {input.normal_bam}) -f {input.fasta} --reference $(readlink -f {params.ref})
        -o {params.out_dir} --chr_prefixed_genome $chr_prefixed --sex $sex --cpu {threads} >> {log.stdout} 2>> {log.stderr} &&  
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" >> {log.stdout}; 
        """)


# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        sub = rules._run_battenberg.output.sub,
        cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2>> {log.stderr}
        """)


# Fill subclones.txt with empty regions for compatibility with downstream tools
rule _battenberg_fill_subclones:
    input:
        sub = str(rules._run_battenberg.output.sub)
    output:
        sub = CFG["dirs"]["fill_regions"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.filled.txt"
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_fill_subclones.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/",
        script = "fill_segments.sh",
        arm_file = lambda w: "src/chromArm.hg38.bed" if "38" in str({w.genome_build}) else "src/chromArm.grch37.bed",
        blacklist_file = lambda w: "src/blacklisted.hg38.bed" if "38" in str({w.genome_build}) else "src/blacklisted.grch37.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.path}{params.script}
        {params.path}{params.arm_file}
        {input.sub}
        {params.path}{params.blacklist_file}
        {output.sub}
        {wildcards.tumour_id}
        subclones
        2>> {log.stderr}
        """)


#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based cleaning in the output directory
rule _battenberg_cleanup:
    input:
        rules._battenberg_to_igv_seg.output.seg
    output:
        complete = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
    group: "battenberg_post_process"
    shell:
        op.as_one_line("""
        d=$(dirname {output});
        rm -f $d/*impute_input* &&
        rm -f $d/*alleleFrequencies* &&
        rm -f $d/*aplotype* &&
        rm -f $d/*BAFsegmented* && 
        touch {output.complete}
        """)


def _battenberg_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

# Convert the coordinates of seg file to a different genome build
rule _battenberg_convert_coordinates:
    input:
        battenberg_native = str(rules._battenberg_to_igv_seg.output.seg),
        battenberg_chain = _battenberg_get_chain
    output:
        battenberg_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.liftover_script}
        SEG
        {input.battenberg_native}
        {output.battenberg_lifted}
        {input.battenberg_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)

# ensure to request the correct files for each projection and drop wildcards that won't be used downstream
def _battenberg_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["battenberg"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    
    prefixed_projections = CFG["options"]["prefixed_projections"]
    non_prefixed_projections = CFG["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        hg38_projection = str(rules._battenberg_to_igv_seg.output.seg).replace("{genome_build}", this_genome_build[0])
        grch37_projection = str(rules._battenberg_convert_coordinates.output.battenberg_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the hg19 (prefixed) separately
        if "38" in str(this_genome_build[0]):
            grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
        else:
            grch37_projection = grch37_projection.replace("{chain}", "hg19ToHg38")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        grch37_projection = str(rules._battenberg_to_igv_seg.output.seg).replace("{genome_build}", this_genome_build[0])
        hg38_projection = str(rules._battenberg_convert_coordinates.output.battenberg_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the grch38 (non-prefixed) separately
        if "38" in str(this_genome_build[0]):
            hg38_projection = hg38_projection.replace("{chain}", "hg38ToHg19")
        else:
            hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")

    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Fill segments of both native and filled file
rule _battenberg_fill_segments:
    input:
        unpack(_battenberg_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/"
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)


def _battenberg_determine_projection(wildcards):
    CFG = config["lcr-modules"]["battenberg"]
    if any(substring in wildcards.projection for substring in ["hg19", "grch37", "hs37d5"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"
    elif any(substring in wildcards.projection for substring in ["hg38", "grch38"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    return (this_file)


# Normalize chr prefix of the output file
rule _battenberg_normalize_projection:
    input:
        filled = _battenberg_determine_projection,
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_battenberg"]
    threads: 1
    group: "battenberg_post_process"
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False)
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _battenberg_output_projection:
    input:
        projection = str(rules._battenberg_normalize_projection.output.projection)
    output:
        projection = CFG["output"]["seg"]["projection"]
    threads: 1
    group: "battenberg_post_process"
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience

rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        sub = rules._battenberg_fill_subclones.output.sub,
        cp = rules._run_battenberg.output.cp
    output:
        seg = CFG["output"]["seg"]["original"],
        sub = CFG["output"]["txt"]["subclones"],
        cp = CFG["output"]["txt"]["cell_ploidy"]
    params: 
        batt_dir = CFG["dirs"]["battenberg"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        png_dir = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}"
    group: "battenberg_post_process"
    run:
        plots = glob.glob(params.batt_dir + "/*.png")
        for png in plots:
            bn = os.path.basename(png)
            op.relative_symlink(png, params.png_dir + "/" + bn,in_module=True)
        op.relative_symlink(input.seg, output.seg,in_module=True)
        op.relative_symlink(input.sub, output.sub,in_module=True)
        op.relative_symlink(input.cp, output.cp,in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [ 
                rules._battenberg_output_seg.output.sub,
                rules._battenberg_output_seg.output.seg,
                rules._battenberg_cleanup.output.complete
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]),
        expand(
            expand(
            [
                str(rules._battenberg_output_projection.output.projection)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            pair_status=CFG["runs"]["pair_status"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=["battenberg"] * len(CFG["runs"]["tumour_sample_id"]),
            allow_missing=True),
            projection=CFG["output"]["requested_projections"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
