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


<<<<<<< HEAD
=======
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


>>>>>>> master
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

<<<<<<< HEAD


##### PICARD METRICS
rule create_seq_dict:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        seq_dict = "genomes/{genome_build}/genome_fasta/genome.dict"
    log: 
        "genomes/{genome_build}/genome_fasta/genome_dict.log"
    conda: CONDA_ENVS["picard"]
    shell:
        op.as_one_line("""
        picard CreateSequenceDictionary
        R={input.fasta}
        O={output.seq_dict}
        2> {log}
        &&
        chmod a-w {output.seq_dict}
        """)


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
=======
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
>>>>>>> master
        """)
