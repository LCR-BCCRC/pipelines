#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper
# Module Author:    Jasper
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
# `CFG` is a shortcut to `config["lcr-modules"]["controlfreec"]`
CFG = op.setup_module(
    name = "controlfreec",
    version = "1.2",
    subdirectories = ["inputs", "mpileup", "run", "outputs"]
)



# Define rules to be run locally when using a compute cluster
localrules:
    _controlfreec_input_bam,
    _controlfreec_config,
    _controlfreec_output,
    _controlfreec_all


##### RULES #####

#### generate references ####
rule _controlfreec_get_map_refs:
    output:
        tar = temp(CFG["dirs"]["inputs"] + "references/{genome_build}/freec/out100m2_{genome_build}.tar.gz"),
        gem = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/out100m2_{genome_build}.gem"
    params:
        provider = "ensembl",
        url = lambda w: {"grch37": "http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz",
                        "hg19": "http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz",
                        "grch38": "http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip",
                        "hg38": "http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip"}[w.genome_build],
        command1 = lambda w: {"grch37": "tar -xvf ",
                            "hg19": "tar -xvf ",
                            "grch38": "unzip ",
                            "hg38": "unzip "}[w.genome_build],
        command2 = lambda w: {"grch37": " --wildcards --no-anchored 'out100m2*gem' && mv out100m2_hg19.gem ",
                            "hg19": " --wildcards --no-anchored 'out100m2*gem' && mv out100m2_hg19.gem ",
                            "grch38": " -d ",
                            "hg38": " -d "}[w.genome_build],
        outdir = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/"
    shell:
        "wget -O {output.tar} {params.url} "
        "&& {params.command1} {output.tar} {params.command2} {params.outdir}"

# no chr for grch37 and grch38
# chr for hg19 and hg38
# Symlink chromosomes used (i.e. chr1-22,X,Y)
checkpoint _controlfreec_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "references/{genome_build}/main_chromosomes_withY.txt"
    run:
        op.absolute_symlink(input.chrs, output.chrs)


def _controlfreec_get_chr_fastas(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    chrs = checkpoints._controlfreec_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chromosome = file.read().rstrip("\n").split("\n")
    fastas = expand(
        CFG["dirs"]["inputs"] +  "references/{{genome_build}}/freec/chr/{chromosome}.fa",
        chromosome = chromosome
    )
    return(fastas)


#generates file with chromomsome lengths from genome.fa.fai
rule _controlfreec_generate_chrLen:
    input:
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/{genome_build}.len"
    shell:
        op.as_one_line("""
            grep -P '^chr[0-9,X,Y]+\t|^[0-9,X,Y]' {input.fai} | awk '{{print $1"\t"$2}}' > {output.chrLen}
        """)

#generates a per-chromosome fasta file from genome.fa
rule _controlfreec_generate_chrFasta:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        fasta = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/{chromosome}.fa"
    conda:
        CFG["conda_envs"]["controlfreec"]
    shell:
        "samtools faidx {input.fasta} {wildcards.chromosome} > {output.fasta} "

#checks that all chroms are accounted for
rule _controlfreec_check_chrFiles:
    input:
        _controlfreec_get_chr_fastas
    output:
        touch(CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/.all_done")


rule _controlfreec_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/dbsnp.common_all-151.bed"
    shell:
        op.as_one_line(""" gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed} """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _controlfreec_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


#### set-up mpileups for BAF calling ####
def _controlfreec_get_chr_mpileups(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    chrs = checkpoints._controlfreec_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["mpileup"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}.{chrom}.minipileup.pileup.gz", 
        chrom = chrs
    )
    return(mpileups)


rule _controlfreec_mpileup_per_chrom:
    input:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        fastaFile = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bed = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/dbsnp.common_all-151.bed"
    output: # creates a temporary file for mpileup
        pileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.minipileup.pileup.gz")
    conda:
        CFG["conda_envs"]["controlfreec"]
    resources: 
        **CFG["resources"]["mpileup"]
    group: "controlfreec"
    log:
        stderr = CFG["logs"]["inputs"] + "mpileup/{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup.stderr.log",
    shell:
        "samtools mpileup -l {input.bed} -r {wildcards.chrom} -Q 20 -f {input.fastaFile} {input.bam} | gzip -c > {output.pileup} 2> {log.stderr}"


rule _controlfreec_concatenate_pileups:
    input: 
        _controlfreec_get_chr_mpileups
    output: 
        mpileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam_minipileup.pileup.gz")
    resources: 
        **CFG["resources"]["cat"]
    group: "controlfreec"
    shell: 
        "cat {input} > {output.mpileup} "


#### Run control-FREEC ####

# set-up controlfreec
rule _controlfreec_config:
    input:
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
        fastaFile = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        reference = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/out100m2_{genome_build}.gem",
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/{genome_build}.len",
        done = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/.all_done"
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/config_WGS.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    params:
        config = CFG["options"]["configFile"],
        dbSNP = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        shiftInQuality = CFG["options"]["shiftInQuality"],
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/",
        window = CFG["options"]["window"],
        ploidy = CFG["options"]["ploidy"],
        breakPointValue = CFG["options"]["breakPointThreshold"],
        breakPointType = CFG["options"]["breakPointType"],
        coefVar = CFG["options"]["coefficientOfVariation"],
        numCon = CFG["options"]["contamination"],
        booCon = CFG["options"]["contaminationAdjustment"],
        bedGraphOutput = CFG["options"]["BedGraphOutput"],
        degree = CFG["options"]["degree"],
        forceGC = CFG["options"]["forceGCcontentNormalization"],
        chrFiles = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/",
        intercept = CFG["options"]["intercept"],
        minCNAlength = CFG["options"]["minCNAlength"],
        minimalCoveragePerPosition = CFG["options"]["minimalCoveragePerPosition"],
        minimalQualityPerPosition = CFG["options"]["minQualityPerPosition"],
        minMapPerWindow = CFG["options"]["minMappabilityPerWindow"],
        minimumSubclonePresence = CFG["options"]["minimalSubclonePresence"],
        naBoo = CFG["options"]["printNA"],
        noisyData = CFG["options"]["noisyData"],
        step = CFG["options"]["step"],
        telocentromeric = CFG["options"]["telocentromeric"],
        threads = CFG["threads"]["controlfreec_run"],
        uniqBoo = CFG["options"]["uniqueMatch"]
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.tumour_bam}|g\" {params.config} | "
        "sed \"s|CONTROLFILE|{input.normal_bam}|g\" | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|fastaPath|{input.fastaFile}|g\" | "
        "sed \"s|DBsnpFile|{params.dbSNP}|g\" | "
        "sed \"s|phredQuality|{params.shiftInQuality}|g\" | "
        "sed \"s|windowSize|{params.window}|g\" | "
        "sed \"s|ploidyInput|{params.ploidy}|g\" | "
        "sed \"s|chrLenFileInput|{input.chrLen}|g\" | "
        "sed \"s|chrFilesPath|{params.chrFiles}|g\" | "
        "sed \"s|sambambaPath|$sambambaPathName|g\" | "
        "sed \"s|bedtoolsPath|$bedtoolsPathName|g\" | "
        "sed \"s|samtoolsPath|$samtoolsPathName|g\" | "
        "sed \"s|breakPointValue|{params.breakPointValue}|g\" | "
        "sed \"s|breakPointTypeBe|{params.breakPointType}|g\" | "
        "sed \"s|bedgraphBoo|{params.bedGraphOutput}|g\" | "
        "sed \"s|coefVar|{params.coefVar}|g\" | "
        "sed \"s|numCon|{params.numCon}|g\" | "
        "sed \"s|booCon|{params.booCon}|g\" | "
        "sed \"s|forceGCvalue|{params.forceGC}|g\" | "
        "sed \"s|interPoly|{params.intercept}|g\" | "
        "sed \"s|numDegree|{params.degree}|g\" | "
        "sed \"s|minCNAvalue|{params.minCNAlength}|g\" | "
        "sed \"s|minMapPerWindow|{params.minMapPerWindow}|g\" | "
        "sed \"s|minimumSubclonePresenceValue|{params.minimumSubclonePresence}|g\" | "
        "sed \"s|minCovPerPos|{params.minimalCoveragePerPosition}|g\" | "
        "sed \"s|minQualPerPos|{params.minimalQualityPerPosition}|g\" | "
        "sed \"s|booNoise|{params.noisyData}|g\" | "
        "sed \"s|stepValue|{params.step}|g\" | "
        "sed \"s|teloValue|{params.telocentromeric}|g\" | "
        "sed \"s|uniqBoo|{params.uniqBoo}|g\" | "
        "sed \"s|naBoo|{params.naBoo}|g\" | "
        "sed \"s|numThreads|{params.threads}|g\" | "
        "sed \"s|referenceFile|{input.reference}|g\" > {output}"


rule _controlfreec_run:
    input:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/config_WGS.txt",
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
    output:
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
        ratio = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        CNV = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
        BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
    conda: CFG["conda_envs"]["controlfreec"]
    threads: CFG["threads"]["controlfreec_run"]
    resources: **CFG["resources"]["controlfreec_run"]
    log:
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run.stderr.log"
    shell:
        "freec -conf {input.config} > {log.stdout} 2> {log.stderr} "


rule _controlfreec_calc_sig:
    input:
        CNVs = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
    output:
        txt = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"]
    threads: CFG["threads"]["calc_sig"]
    resources: **CFG["resources"]["calc_sig"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:         
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/calc_sig.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/calc_sig.stderr.log"
    shell:
        "cat {params.calc_sig} | R --slave --args {input.CNVs} {input.ratios} > {log.stdout} 2> {log.stderr}"


rule _controlfreec_plot:
    input:
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt",
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_info.txt"
    output:
        plot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.png",
        log2plot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.log2.png",
        bafplot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt.png"
    params:
        plot = CFG["software"]["FREEC_graph"]
    threads: CFG["threads"]["plot"]
    resources: **CFG["resources"]["plot"]
    conda: CFG["conda_envs"]["controlfreec"]
    log: 
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot.stderr.log"
    shell:
        "cat {params.plot} | R --slave --args `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} {input.BAF} > {log.stdout} 2> {log.stderr} "


rule _controlfreec_freec2bed:
    input:
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_info.txt"
    output:
        bed = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bed"
    params:
        freec2bed = CFG["software"]["freec2bed"]
    threads: CFG["threads"]["freec2bed"]
    resources: **CFG["resources"]["freec2bed"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/freec2bed.stderr.log"
    shell:
        "ploidy=$(grep Output_Ploidy {input.info} | cut -f 2); "
        "perl {params.freec2bed} -f {input.ratios} -p $ploidy > {output.bed} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _controlfreec_output:
    input:
        plot = str(rules._controlfreec_plot.output.plot),
        log2plot = str(rules._controlfreec_plot.output.log2plot),
        CNV = str(rules._controlfreec_calc_sig.output.txt),
        bed = str(rules._controlfreec_freec2bed.output.bed),
        BAF = str(rules._controlfreec_run.output.BAF),
        BAFgraph = str(rules._controlfreec_plot.output.bafplot),
        ratio = str(rules._controlfreec_run.output.ratio)
    output:
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plots/{tumour_id}--{normal_id}--{pair_status}.ratio.png",
        log2plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/log2plots/{tumour_id}--{normal_id}--{pair_status}.ratio.log2.png",
        CNV = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/CNV/{tumour_id}--{normal_id}--{pair_status}.CNVs.p.value.txt",
        bed = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/bed/{tumour_id}--{normal_id}--{pair_status}.CNVs.bed",
        BAF = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/BAF/{tumour_id}--{normal_id}--{pair_status}.BAF.txt",
        BAFgraph = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/BAFplot/{tumour_id}--{normal_id}--{pair_status}.BAF.png",
        ratio = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/ratio/{tumour_id}--{normal_id}--{pair_status}.ratio.txt"
    run:
        op.relative_symlink(input.plot, output.plot, in_module=True)
        op.relative_symlink(input.log2plot, output.log2plot, in_module=True)
        op.relative_symlink(input.CNV, output.CNV, in_module=True)
        op.relative_symlink(input.bed, output.bed, in_module=True)
        op.relative_symlink(input.BAF, output.BAF, in_module=True)
        op.relative_symlink(input.BAFgraph, output.BAFgraph, in_module=True)
        op.relative_symlink(input.ratio, output.ratio, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _controlfreec_all:
    input:
        expand(
            [
                str(rules._controlfreec_output.output.plot),
                str(rules._controlfreec_output.output.log2plot),
                str(rules._controlfreec_output.output.CNV),
                str(rules._controlfreec_output.output.bed),
                str(rules._controlfreec_output.output.BAF),
                str(rules._controlfreec_output.output.BAFgraph),
                str(rules._controlfreec_output.output.ratio)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            pair_status=CFG["runs"]["pair_status"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"])



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
