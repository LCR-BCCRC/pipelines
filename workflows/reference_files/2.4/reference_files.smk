##### HEADER #####


include: "reference_files_header.smk"


##### SRC/ETC DIRECTORY #####


for k,v in config['pon'].items():
    config['pon'][k] = os.path.join(workflow.basedir, v)

for k,v in config['jabba'].items():
    config['jabba'][k] = os.path.join(workflow.basedir, v)


##### SEQUENCE AND INDICES #####


rule get_genome_fasta_download:
    input: 
        fasta = rules.download_genome_fasta.output.fasta
    output: 
        fasta = "genomes/{genome_build}/genome_fasta/genome.fa"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.fasta} {output.fasta}"


rule index_genome_fasta:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        fai = "genomes/{genome_build}/genome_fasta/genome.fa.fai"
    log: 
        "genomes/{genome_build}/genome_fasta/genome.fa.fai.log"
    conda: CONDA_ENVS["samtools"]
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"


rule create_bwa_index:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        prefix = touch("genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa")
    log: 
        "genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa.log"
    conda: CONDA_ENVS["bwa"]
    resources:
        mem_mb = 20000
    shell:
        "bwa index -p {output.prefix} {input.fasta} > {log} 2>&1"


rule create_gatk_dict:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        fai = rules.index_genome_fasta.output.fai
    output:
        dict = "genomes/{genome_build}/genome_fasta/genome.dict"
    log:
        "genomes/{genome_build}/gatk_fasta/genome.dict.log"
    conda: CONDA_ENVS["gatk"]
    resources:
        mem_mb = 20000
    shell:
        op.as_one_line(""" 
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict} > {log} 2>&1
        """)


rule create_star_index:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output: 
        index = directory("genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}")
    log: 
        "genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}.log"
    conda: CONDA_ENVS["star"]
    threads: 12
    resources:
        mem_mb = 42000
    shell:
        op.as_one_line("""
        mkdir -p {output.index}
            &&
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index}
        --genomeFastaFiles {input.fasta} --sjdbOverhang {wildcards.star_overhang}
        --sjdbGTFfile {input.gtf} --outTmpDir {output.index}/_STARtmp
        --outFileNamePrefix {output.index}/ > {log} 2>&1
        """)


rule get_liftover_chains:
    input:
        chains = rules.download_liftover_chains.output.chains
    output:
        chains = "genomes/{genome_build}/chains/{version}/{chain_version}.over.chain"
    shell:
        "ln -srf {input.chains} {output.chains}"

rule get_sdf_refs: 
    input: 
        sdf = ancient(rules.download_sdf.output.sdf)
    output: 
        sdf = directory("genomes/{genome_build}/sdf")
    wildcard_constraints: 
        genome_build = "hg38|hg19|grch37|hs37d5"
    shell: 
        "ln -srfT {input.sdf} {output.sdf}"



##### METADATA #####


rule store_genome_build_info:
    output: 
        version = "genomes/{genome_build}/version.txt",
        provider = "genomes/{genome_build}/provider.txt"
    params:
        version = lambda w: config["genome_builds"][w.genome_build]["version"],
        provider = lambda w: config["genome_builds"][w.genome_build]["provider"]
    shell: 
        op.as_one_line("""
        echo "{params.version}" > {output.version}
            &&
        echo "{params.provider}" > {output.provider}
        """)


rule get_main_chromosomes_download:
    input: 
        txt = get_download_file(rules.download_main_chromosomes.output.txt),
        chrx = get_download_file(rules.download_chromosome_x.output.txt),
        fai = rules.index_genome_fasta.output.fai
    output: 
        txt = "genomes/{genome_build}/genome_fasta/main_chromosomes.txt",
        bed = "genomes/{genome_build}/genome_fasta/main_chromosomes.bed",
        chrx = "genomes/{genome_build}/genome_fasta/chromosome_x.txt",
        patterns = temp("genomes/{genome_build}/genome_fasta/main_chromosomes.patterns.txt")
    conda: CONDA_ENVS["coreutils"]
    shell: 
        op.as_one_line("""
        sed 's/^/^/' {input.txt} > {output.patterns}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        cut -f1 > {output.txt}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        awk 'BEGIN {{FS=OFS="\t"}} {{print $1,  0, $2}}' > {output.bed}
            &&
        ln -srf {input.chrx} {output.chrx}
        """)


rule get_main_chromosomes_withY_download:
    input: 
        txt = get_download_file(rules.download_main_chromosomes_withY.output.txt),
        fai = rules.index_genome_fasta.output.fai
    output: 
        txt = "genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt",
        bed = "genomes/{genome_build}/genome_fasta/main_chromosomes_withY.bed",
        patterns = temp("genomes/{genome_build}/genome_fasta/main_chromosomes.patterns.txt")
    conda: CONDA_ENVS["coreutils"]
    shell: 
        op.as_one_line("""
        sed 's/^/^/' {input.txt} > {output.patterns}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        cut -f1 > {output.txt}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        awk 'BEGIN {{FS=OFS="\t"}} {{print $1,  0, $2}}' > {output.bed}
        """)

##### ANNOTATIONS #####


rule get_gencode_download: 
    input:
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output:
        gtf = "genomes/{genome_build}/annotations/gencode_annotation-{gencode_release}.gtf"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.gtf} {output.gtf}"


rule calc_gc_content:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta
    output:
        wig = "genomes/{genome_build}/annotations/gc_wiggle.window_{gc_window_size}.wig.gz"
    log:
        "genomes/{genome_build}/annotations/gc_wiggle.window_{gc_window_size}.wig.gz.log"
    conda: CONDA_ENVS["sequenza-utils"]
    shell:
        op.as_one_line("""
        sequenza-utils gc_wiggle --fasta {input.fasta} -w {wildcards.gc_window_size} -o -
            |
        gzip -c > {output.wig}
        """)

rule get_jabba_gc_rds:
    input:
        bed = get_download_file(rules.download_ucsc_gc.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/gc1000.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line("""
        Rscript
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'names(mcols(gr))[1] <- "score"; gr[gr$score == "nan",]$score <- 0'
            -e 'gr$score <- as.numeric(gr$score)'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)

rule get_jabba_map_rds:
    input:
        bed = get_download_file(rules.download_ucsc_map.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/map1000.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line("""
        Rscript
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'names(mcols(gr))[1] <- "score"; gr[gr$score == "nan",]$score <- 0'
            -e 'gr$score <- as.numeric(gr$score)'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)

rule get_par_rds:
    input:
        bed = get_download_file(rules.download_par_bed.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/PAR_{genome_build}.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line(""" 
        Rscript 
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)


##### JaBbA PON generation ######


rule install_fragcounter:
    output:
        complete = "downloads/jabba_prereqs/fragcounter.installed"
    conda: CONDA_ENVS["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'if (!"fragCounter" %in% rownames(installed.packages())) {{remotes::install_github("morinlab/fragCounter", upgrade = TRUE, force = TRUE)}}'
                -e 'library(fragCounter)'
            &&
        touch {output.complete}
        """)


rule install_dryclean:
    input:
        installed = "downloads/jabba_prereqs/fragcounter.installed"
    output:
        complete = "downloads/jabba_prereqs/dryclean.installed"
    conda: CONDA_ENVS["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'if (!"gUtils" %in% rownames(installed.packages())) {{remotes::install_github("mskilab/gUtils", upgrade = FALSE)}}'
                -e 'if (!"dryclean" %in% rownames(installed.packages())) {{remotes::install_github("morinlab/dryclean", upgrade = FALSE)}}'
                -e 'library(dryclean)'
            &&
        touch {output.complete}
        """)

# Use normal samples specified in config['pon'] table
NORMALS = {k : op.load_samples(v) for k,v in config['pon'].items()}
NORMALS = {k : dict(zip(v.sample_id, v.data_path)) for k,v in NORMALS.items()}


rule jabba_pon_symlink_normal_bams:
    output: 
        bam = "genomes/{genome_build}/jabba/pon/00-normal_bams/{sample_id}.bam",
        bai = "genomes/{genome_build}/jabba/pon/00-normal_bams/{sample_id}.bam.bai"
    params:
        bam = lambda wc: NORMALS[wc.genome_build][wc.sample_id]
    shell: 
        "ln -sf {params.bam} {output.bam}; "
        "ln -sf {params.bam}.bai {output.bai}"

# Run fragcounter on normals to get GC/mappability corrected coverage values
rule jabba_pon_run_fragcounter:
    input:
        installed = str(rules.install_fragcounter.output.complete),
        bam = str(rules.jabba_pon_symlink_normal_bams.output.bam),
        gc = str(rules.get_jabba_gc_rds.output.rds),
        map = str(rules.get_jabba_map_rds.output.rds),
        ref = str(rules.download_genome_fasta.output.fasta) 
    output:
        rds = "genomes/{genome_build}/jabba/pon/01-fragcounter/run/{sample_id}/cov.rds"
    conda: CONDA_ENVS["jabba"]
    threads: 1
    resources:
        mem_mb = 4096
    shell:
        op.as_one_line("""
        FRAG=$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))'); $FRAG/frag -b {input.bam} -r {input.ref} -w 1000 -d `dirname {input.gc}` -o `dirname {output.rds}`
        """)


rule jabba_pon_symlink_fragcounter:
    input:
        rds = str(rules.jabba_pon_run_fragcounter.output.rds)
    output:
        rds = "genomes/{genome_build}/jabba/pon/01-fragcounter/cov/{sample_id}.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Use dryclean to create panel-of-normal (PON) from fragcounter coverages
rule jabba_pon_make_pon:
    input:
        installed = str(rules.install_dryclean.output.complete),
        par = str(rules.get_par_rds.output.rds),
        rds = lambda wc: expand("genomes/{{genome_build}}/jabba/pon/01-fragcounter/cov/{sample_id}.cov.rds",
              sample_id = NORMALS[wc.genome_build].keys()),
        make_pon = config['jabba']['make_pon']
    output:
        tbl = "genomes/{genome_build}/jabba/pon/normal_table.rds",
        pon = "genomes/{genome_build}/jabba/pon/detergent.rds"
    params:
        choose_samples = 'cluster'
    conda: CONDA_ENVS["jabba"]
    threads: 10
    resources:
        mem_mb = 30000
    shell:
        op.as_one_line("""
        Rscript {input.make_pon} {threads} `dirname {input.rds[0]}` {output.tbl} {output.pon} {input.par} {wildcards.genome_build} {params.choose_samples}
        """)

# Improve coverage signal fidelity by running each normal against PON
# Serves to identify and separate background variations (noise) from foreground variation (signal)
rule jabba_pon_run_dryclean_normal:
    input:
        installed = str(rules.install_dryclean.output.complete),
        rds = str(rules.jabba_pon_run_fragcounter.output.rds),
        pon = str(rules.jabba_pon_make_pon.output.pon)
    output:
        rds = "genomes/{genome_build}/jabba/pon/02-dryclean/run/{sample_id}/drycleaned.cov.rds"
    conda: CONDA_ENVS["jabba"]
    threads: 12
    resources:
       mem_mb = 30000
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, mc.cores = {threads}, germline.filter = FALSE)'
                -e 'saveRDS(decomp, "{output.rds}")'
        """)


rule jabba_pon_link_dryclean_normal_rds:
    input:
        rds = "genomes/{genome_build}/jabba/pon/02-dryclean/run/{sample_id}/drycleaned.cov.rds"
    output:
        rds = "genomes/{genome_build}/jabba/pon/02-dryclean/cov/{sample_id}.drycleaned.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Use dryclean normals to identify germline regions which will be used
# on tumour data to find germline regions.
rule jabba_pon_make_germline_filter:
    input:
        installed = str(rules.install_dryclean.output.complete),
        tbl = str(rules.jabba_pon_make_pon.output.tbl),
        rds = lambda wc: expand("genomes/{{genome_build}}/jabba/pon/02-dryclean/cov/{sample_id}.drycleaned.cov.rds", sample_id = NORMALS[wc.genome_build].keys()),
        make_germline = config["jabba"]["make_germline"]
    output:
        germline = "genomes/{genome_build}/jabba/pon/germline.markers.rds"
    conda: CONDA_ENVS["jabba"]
    threads: 25
    resources:
        mem_mb = 30000
    shell:
        op.as_one_line("""
        Rscript {input.make_germline} {input.tbl} `dirname {input.rds[0]}` {output.germline} 0.5 0.1 {threads}
        """)


##### VARIATION #####


rule get_dbsnp_download: 
    input:
        vcf = get_download_file(rules.download_dbsnp_vcf.output.vcf)
    output:
        vcf = "genomes/{genome_build}/variation/dbsnp.common_all-{dbsnp_build}.vcf.gz"
    conda: CONDA_ENVS["samtools"]
    shell:
        op.as_one_line("""
        bgzip -c {input.vcf} > {output.vcf}
            &&
        tabix {output.vcf}
        """)

##### PICARD METRICS

rule create_rRNA_interval:
    input:
        gtf = rules.get_gencode_download.output.gtf
    output: 
        rrna_int = "genomes/{genome_build}/rrna_intervals/rRNA_int_gencode-{gencode_release}.txt"
    log: 
        "genomes/{genome_build}/rrna_intervals/rRNA_int_gencode-{gencode_release}.log"
    conda: CONDA_ENVS["picard"]
    shell:
        op.as_one_line("""
        grep 'gene_type "rRNA"' {input.gtf} |
        awk '$3 == "transcript"' |
        cut -f1,4,5,7,9 |
        perl -lane '
            /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
            print join "\t", (@F[0,1,2,3], $1)
        ' | 
        sort -k1V -k2n -k3n >> {output.rrna_int}
        &&
        chmod a-w {output.rrna_int}
        """)


rule create_refFlat:
    input:
        gtf = rules.get_gencode_download.output.gtf
    output:
        txt = "genomes/{genome_build}/annotations/refFlat_gencode-{gencode_release}.txt"
    log: "genomes/{genome_build}/annotations/gtfToGenePred-{gencode_release}.log"
    conda: CONDA_ENVS["ucsc-gtftogenepred"]
    threads: 4
    resources:
        mem_mb = 6000
    shell:
        op.as_one_line("""
        gtfToGenePred -genePredExt -geneNameAsName2 
        {input.gtf} {output.txt}.tmp 
        2> {log} &&
        paste <(cut -f 12 {output.txt}.tmp) <(cut -f 1-10 {output.txt}.tmp) > {output.txt}
        """)

##### ENCODE #####

rule get_blacklist_download: 
    input: 
        bed = get_download_file(rules.download_blacklist.output.bed)
    output: 
        bed = "genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.bed} {output.bed}"

##### REPEATMASKER #####

rule get_repeatmasker_download:
    input: 
        bed = get_download_file(rules.download_repeatmasker.output.bed)
    output: 
        bed = "genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed"
    conda: CONDA_ENVS["coreutils"]
    shell: 
        "ln -srf {input.bed} {output.bed}"

# salmon index
rule _download_salmon_script:
    output: 
        script = "downloads/scripts/salmon/generateDecoyTranscriptome.sh"
    log: 
        "downloads/scripts/salmon/log"
    params: 
        url = "https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh"
    shell:
        op.as_one_line("""
        curl -L {params.url} > {output.script} 2> {log}
            &&
        chmod a-w {output.script}
        """)


rule _create_transcriptome_fasta:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file("downloads/gencode-33/gencode.annotation.{version}.gtf")
    output: 
        fasta = "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/transcriptome.fa"
    log: 
        "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/transcriptome.log"
    conda: CONDA_ENVS["gffread"]
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        op.as_one_line("""
        gffread -F
        -w {output.fasta}
        -g {input.fasta}
        {input.gtf}
        > {log} 2>&1
            &&
        chmod a-w {output.fasta}
        """)

rule create_salmon_index:
    input:
        fasta = rules._create_transcriptome_fasta.output.fasta,
        gtf = get_download_file("downloads/gencode-33/gencode.annotation.{version}.gtf")
    output: 
        index = directory("genomes/{genome_build}/salmon_index/salmon-{salmon_version}/index")
    log: 
        "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/log"
    conda: CONDA_ENVS["salmon"]
    threads: 8
    resources:
        mem_mb = 12000
    shell:
        op.as_one_line("""
        mkdir {output.index}
            &&
        chmod u+w {output.index}
            &&
        salmon index
        --threads {threads}
        -t {input.fasta}
        -i {output.index}
        > {log} 2>&1
        """)

rule get_af_only_gnomad_vcf:
    input:
        vcf = get_download_file(rules.download_af_only_gnomad_vcf.output.vcf)
    output:
        vcf = "genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"
    conda: CONDA_ENVS["samtools"]
    shell:
        op.as_one_line(""" 
        bgzip -c {input.vcf} > {output.vcf}
            &&
        tabix {output.vcf}
        """)

rule get_mutect2_pon:
    input:
        vcf = get_download_file(rules.download_mutect2_pon.output.vcf)
    output:
        vcf = "genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.gz"
    conda: CONDA_ENVS["samtools"]
    shell:
        op.as_one_line(""" 
        bgzip -c {input.vcf} > {output.vcf}
            &&
        tabix {output.vcf}
        """)

rule get_mutect2_small_exac:
    input:
        vcf = get_download_file(rules.download_mutect2_small_exac.output.vcf)
    output:
        vcf = "genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz"
    conda: CONDA_ENVS["samtools"]
    shell:
        op.as_one_line(""" 
        bgzip -c {input.vcf} > {output.vcf}
            &&
        tabix {output.vcf}
        """)
