#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
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
# `CFG` is a shortcut to `config["lcr-modules"]["jabba"]`
CFG = op.setup_module(
    name = "jabba",
    version = "1.0",
    subdirectories = ["inputs", "fragcounter", "dryclean", "jabba",  "outputs"],
)

os.environ["CPLEX_DIR"] = str(CFG["CPLEX_DIR"])


# Define rules to be run locally when using a compute cluster
localrules:
    _jabba_install_jabba,
    _jabba_input_bam,
    _jabba_link_graph_rds,
    _jabba_all


##### RULES #####

### Adjust coverage for GC and mappability ###

# Install JaBba suite from github
rule _jabba_install_fragcounter:
    output:
        complete = CFG["dirs"]["fragcounter"] + "fragcounter.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'if (!"fragCounter" %in% rownames(installed.packages())) {{remotes::install_github("morinlab/fragCounter", upgrade = TRUE, force = TRUE)}}'
                -e 'library(fragCounter)'
            &&
        touch {output.complete}
        """)

rule _jabba_install_dryclean:
    input:
        installed1 = str(rules._jabba_install_fragcounter.output.complete)
    output:
        complete = CFG["dirs"]["dryclean"] + "dryclean.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)' 
                -e 'if (!"dryclean" %in% rownames(installed.packages())) {{remotes::install_github("morinlab/dryclean", upgrade = FALSE)}}'
                -e 'library(dryclean)'
            &&
        touch {output.complete}
        """)

rule _jabba_install_jabba:
    input:
        installed1 = str(rules._jabba_install_dryclean.output.complete),
        installed2 = str(rules._jabba_install_fragcounter.output.complete)
    output:
        complete = CFG["dirs"]["jabba"] + "jabba.installed"
    conda: CFG["conda_envs"]["jabba"]
    params:
        cplex_dir = CFG["CPLEX_DIR"]
    shell:
        op.as_one_line(""" 
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'Sys.setenv(CPLEX_DIR = "{params.cplex_dir}")'
                -e 'if ("copynumber" %in% installed.packages()==FALSE) {{BiocManager::install("copynumber")}}'
                -e 'if (is.null(packageDescription("copynumber")$GithubUsername)) {{remotes::install_github("ShixiangWang/copynumber", dependencies = TRUE)}}'
                -e 'if (!"JaBbA" %in% rownames(installed.packages())) {{remotes::install_github("mskilab/JaBbA", upgrade = FALSE)}}'
                -e 'library(JaBbA)'
            &&
        touch {output.complete}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _jabba_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"]+ "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)


rule _jabba_input_gridss_vcf:
    input:
        junc = CFG["inputs"]["gridss_junc"]
    output:
        junc = CFG["dirs"]["inputs"] + "junc/gridss/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf"
    run:
        op.relative_symlink(input.junc, output.junc)


rule _jabba_input_manta_vcf:
    input:
        junc = CFG["inputs"]["manta_junc"]
    output:
        junc = CFG["dirs"]["inputs"] + "junc/manta/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf"
    run:
        op.relative_symlink(input.junc, output.junc)


rule _jabba_merge_svs:
    input:
        installed = str(rules._jabba_install_jabba.output.complete),
        manta = str(rules._jabba_input_manta_vcf.output.junc),
        gridss = str(rules._jabba_input_gridss_vcf.output.junc)
    output:
        junc = CFG["dirs"]["inputs"] + "junc/merged/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.rds"
    log:
        stdout = CFG["logs"]["inputs"] + "run/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merge_svs.stdout.log",
        stderr = CFG["logs"]["inputs"] + "run/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merge_svs.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["merge_svs"]
    params:
        merge_svs = CFG["inputs"]["merge_svs"]
    resources:
        mem_mb = CFG["mem_mb"]["merge_svs"]
    shell:
        op.as_one_line("""
        Rscript {params.merge_svs} --manta {input.manta} --gridss {input.gridss} --genome {wildcards.genome_build} --rds {output} > {log.stdout} 2> {log.stderr}
        """)
#Rscript {input.merge_svs} --manta {input.manta} --gridss {input.gridss} --genome {wildcards.genome_build} --rds {output}
# Rscript {input.merge_svs} --manta {input.manta} --gridss {input.gridss} --genome {wildcards.genome_build} | sort -k1,1V -k2,2n -k4,4V -k5,5n > {output}

# Runs fragcounter on individual samples
rule _jabba_run_fragcounter:
    input:
        installed = str(rules._jabba_install_fragcounter.output.complete),
        bam = str(rules._jabba_input_bam.output.bam),
        gc = reference_files("genomes/{genome_build}/annotations/jabba/gc1000.rds"),
        map = reference_files("genomes/{genome_build}/annotations/jabba/map1000.rds"),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/cov.rds"
    log:
        stdout = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/fc.stdout.log",
        stderr = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/fc.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["fragcounter"]
    resources:
        mem_mb = CFG["mem_mb"]["fragcounter"]
    shell:
        op.as_one_line(""" 
        FRAG=$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))'); $FRAG/frag -b {input.bam} -r {input.ref} -w 1000 -d `dirname {input.gc}` -o `dirname {output.rds}` > {log.stdout} 2> {log.stderr}
        """)


### Run dryclean on tumours to remove background and germline variation ###


rule _jabba_run_dryclean_tumour:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{tumour_id}/cov.rds",
        pon = reference_files("genomes/{genome_build}/jabba/pon/detergent.rds"),
        germline = reference_files("genomes/{genome_build}/jabba/pon/germline.markers.rds")
    output:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}/drycleaned.cov.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["dryclean"]
    resources:
        mem_mb = CFG["mem_mb"]["dryclean"]
    wildcard_constraints:
        tumour_id = "|".join(CFG["runs"]["tumour_sample_id"])
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, mc.cores = {threads}, germline.file = "{input.germline}", germline.filter = TRUE)'
                -e 'saveRDS(decomp, "{output.rds}")' > {log.stdout} 2> {log.stderr}
        """)


#def _get_junc_file(wildcards):
#    CFG = config["lcr-modules"]["jabba"]
#    filename, ext = os.path.splitext(CFG["inputs"]["sample_junc"])
#    out = CFG["dirs"]["inputs"] + "junc/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}" + ext
#    return(out)


rule _jabba_run_jabba:
    input:
        installed = str(rules._jabba_install_jabba.output.complete),
        rds = str(rules._jabba_run_dryclean_tumour.output.rds),
        junc = str(rules._jabba_merge_svs.output.junc),
        script = CFG["inputs"]["run_jabba_main"]
    output:
        rds = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/jabba.simple.gg.rds"
    log:
        stdout = CFG["logs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stderr.log"
    params:
        CPLEX = CFG["CPLEX_DIR"]
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["jabba"]
    resources:
        mem_mb = CFG["mem_mb"]["jabba"]
    shell:
        op.as_one_line(""" 
        JABBA_PATH=$(Rscript -e 'cat(paste0(installed.packages()["JaBbA", "LibPath"], "/JaBbA/extdata/"))'); $JABBA_PATH/jba `readlink -e {input.junc}` `readlink -e {input.rds}` --field foreground --cfield tier --outdir `dirname {output.rds}` --cores {threads} --slack 100 > {log.stdout} 2> {log.stderr}
        """)


rule _jabba_link_graph_rds:
    input:
        rds = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/jabba.simple.gg.rds"
    output:
        rds = CFG["dirs"]["outputs"] + "rds/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.jabba.simple.gg.rds"
    run:
        op.relative_symlink(input.rds, output.rds)


# Generates the target sentinels for each run, which generate the symlinks
rule _jabba_all:
    input:
        expand(
            [
                str(rules._jabba_link_graph_rds.output.rds),
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
