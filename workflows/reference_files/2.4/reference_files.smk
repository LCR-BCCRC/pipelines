##### HEADER #####


include: "reference_files_header.smk"


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


##### SigProfiler #####

rule download_sigprofiler_genome:
    output:
        complete = "downloads/sigprofiler_prereqs/{sigprofiler_build}.installed"
    conda: CONDA_ENVS["sigprofiler"]
    shell:
        op.as_one_line("""
        python -c 'from SigProfilerMatrixGenerator import install as genInstall;
        genInstall.install("{wildcards.sigprofiler_build}", rsync = False, bash = True)'
            &&
        touch {output.complete}
        """)

def get_sigprofiler_genome(wildcards):
    sigprofiler_build = ''
    if wildcards.genome_build in ['grch37','hg19','hs37d5']:
        sigprofiler_build = "GRCh37"
    elif wildcards.genome_build in ['grch38','grch38-legacy','hg38','hg38-panea']:
        sigprofiler_build = "GRCh38"
    return("downloads/sigprofiler_prereqs/" + sigprofiler_build + ".installed")

rule install_sigprofiler_genome:
    input:
        get_sigprofiler_genome
    output:
        complete = "genomes/{genome_build}/sigprofiler_genomes/{genome_build}.installed"
    run:
        op.relative_symlink(input, output.complete)

