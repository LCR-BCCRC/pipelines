##### MODULES #####


import os
import re
import tarfile
import urllib.request
from collections import defaultdict

import yaml
import snakemake as smk
from snakemake.logging import logger

import oncopipe as op


##### CONFIG #####
localrules: download_genome_fasta,
            download_main_chromosomes, 
            download_main_chromosomes_withY, 
            download_gencode_annotation,
            download_blacklist, 
            hardlink_download, 
            update_contig_names,
            get_genome_fasta_download, 
            index_genome_fasta,
            get_main_chromosomes_download,
            get_main_chromosomes_withY_download, 
            get_gencode_download, 
            create_star_index, 
            get_gencode_download,
            download_af_only_gnomad_vcf,
            download_liftover_chains


# Check for genome builds
assert "genome_builds" in config and len(config["genome_builds"]) > 0, (
    "No `genome_builds` in snakemake configuration."
)

# Switch between case for version names
VERSION_UPPER = {
    "grch37": "GRCh37",
    "grch38": "GRCh38",
}

GENOME_VERSION_GROUPS = {}
for genome_build in VERSION_UPPER.keys():
    GENOME_VERSION_GROUPS[genome_build] = []

# Check genome build versions, providers, and genome_fasta
possible_versions = list(VERSION_UPPER.keys())
possible_providers = ["ensembl", "ucsc", "gencode", "ncbi"]
for build_name, build_info in config["genome_builds"].items():
    assert "version" in build_info and build_info["version"] in possible_versions, (
        f"`version` not set for `{build_name}` or `version` not among {possible_versions}."
    )
    GENOME_VERSION_GROUPS[build_info["version"]].append(build_name)
    assert "provider" in build_info and build_info["provider"] in possible_providers, (
        f"`provider` not set for `{build_name}` or `provider` not among {possible_providers}."
    )
    assert "genome_fasta_url" in build_info, f"`genome_fasta_url` not set for `{build_name}`."
    url_code = urllib.request.urlopen(build_info["genome_fasta_url"]).getcode()
    assert url_code == 200, (
        f"Pinging `genome_fasta_url` for {build_name} returned HTTP code {url_code} "
        f"(rather than 200): \n{build_info['genome_fasta_url']}"
    )    

# Check parent genome, provider for the capture space
for build_name, build_info in config["capture_space"].items():
    assert "provider" in build_info and build_info["provider"] in possible_providers, (
        f"`provider` not set for `{build_name}` or `provider` not among {possible_providers}."
        )
    assert "genome" in build_info and build_info["genome"] in possible_versions,(
        f"`genome` not set for `{build_name}` or `genome` not among {possible_versions}." )
    assert "capture_bed_url" in build_info
    url_code = urllib.request.urlopen(build_info["capture_bed_url"]).getcode()
    assert url_code == 200, (
         f"Pinging `genome_fasta_url` for {build_name} returned HTTP code {url_code} "
         f"(rather than 200): \n{build_info['genome_fasta_url']}"
        )
    if "default" in build_info:
        assert build_info["default"].lower() in ["true", "false"], (
            f"true/false required for for \'default\' field"
            )


##### TOOLS #####


CONDA_ENVS = { pkg: pkg_info["conda_env"] for pkg, pkg_info in config["tools"].items() }

TOOL_VERSIONS = { pkg: pkg_info["version"] for pkg, pkg_info in config["tools"].items() }


##### WILDCARD CONSTRAINTS


wildcard_constraints:
    genome_build = "|".join(config["genome_builds"].keys()),
    version = "|".join(VERSION_UPPER.keys()),
    bwa_version = TOOL_VERSIONS["bwa"],
    star_version = TOOL_VERSIONS["star"],
    gencode_release = "|".join(config["wildcard_values"]["gencode_release"]),
    dbsnp_build = "|".join(config["wildcard_values"]["dbsnp_build"]),
    rm_version = "|".join(config["wildcard_values"]["rm_version"])


##### CHROMOSOME MAPPINGS #####


# Define method for download ChromosomeMappings repository
def download_chrom_mappings(chrom_mappings_dir):
    chrom_mappings_url = "https://github.com/BrunoGrandePhD/ChromosomeMappings/archive/master.tar.gz"
    tar_file_path = "ChromosomeMappings.tar.gz"
    os.makedirs("./", exist_ok=True)
    urllib.request.urlretrieve(chrom_mappings_url, tar_file_path)
    zip_file = tarfile.open(tar_file_path, "r:gz")
    zip_file.extractall("./")
    zip_file.close()
    os.remove(tar_file_path)
    os.rename("ChromosomeMappings-master", chrom_mappings_dir)

# Download ChromosomeMappings if not already done
CHROM_MAPPINGS_DIR = "chrom_mappings"
if not os.path.exists(CHROM_MAPPINGS_DIR):
    download_chrom_mappings(CHROM_MAPPINGS_DIR)

# Find chromosome mappings for human genome
CHROM_MAPPINGS_FILES = os.listdir(CHROM_MAPPINGS_DIR)
CHROM_MAPPINGS_FILES = [ f for f in CHROM_MAPPINGS_FILES if f.startswith("GRCh") ]

TO_PROVIDERS = set()
CHROM_MAPPINGS = defaultdict(lambda: defaultdict(set))
for chrom_map_file in CHROM_MAPPINGS_FILES:
    chrom_map_file_split = re.split("[_2.]", chrom_map_file)
    version, from_provider, to_provider, _ = [ x.lower() for x in chrom_map_file_split ]
    CHROM_MAPPINGS[version][to_provider].add(from_provider)
    # Include the provider itself since you can always just use a file as is
    CHROM_MAPPINGS[version][to_provider].add(to_provider)
    TO_PROVIDERS.add(to_provider)


##### DOWNLOAD #####


rule download_genome_fasta:
    output: 
        fasta = "downloads/genome_fasta/{genome_build}.fa"
    log: 
        "downloads/genome_fasta/{genome_build}.fa.log"
    params: 
        url = lambda w: config["genome_builds"][w.genome_build]["genome_fasta_url"]
    shell:
        "curl -L {params.url} > {output.fasta} 2> {log}"


rule download_main_chromosomes:
    input:
        mapping = lambda w: f"{CHROM_MAPPINGS_DIR}/{VERSION_UPPER[w.version]}_ensembl2ucsc.txt"
    output:
        txt = "downloads/main_chromosomes/main_chromosomes.{version}.txt"
    params:
        provider = "ensembl"
    shell:
        op.as_one_line("""
        egrep -w "^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X)" {input.mapping}
            |
        cut -f1 > {output.txt}
        """)

rule download_main_chromosomes_withY:
    input:
        mapping = lambda w: f"{CHROM_MAPPINGS_DIR}/{VERSION_UPPER[w.version]}_ensembl2ucsc.txt"
    output:
        txt = "downloads/main_chromosomes/main_chromosomes_withY.{version}.txt"
    params:
        provider = "ensembl"
    shell:
        op.as_one_line("""
        egrep -w "^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)" {input.mapping}
            |
        cut -f1 > {output.txt}
        """)

rule download_chromosome_x:
    output:
        txt = "downloads/main_chromosomes/chromosome_x.{version}.txt"
    params:
        provider = "ensembl"
    shell:
        "echo 'X' > {output.txt}"


rule download_gencode_annotation:
    output:
        gtf = "downloads/gencode-{gencode_release}/gencode.annotation.{version}.gtf"
    params:
        provider = "ucsc"
    run:
        url_parts = [
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human",
            f"release_{wildcards.gencode_release}"
        ]
        release_fmt = f"v{wildcards.gencode_release}"
        if wildcards.version == "grch37":
            url_parts.append("GRCh37_mapping")
            release_fmt += "lift37"
        url_parts.append(f"gencode.{release_fmt}.annotation.gtf.gz")
        url = "/".join(url_parts)
        urllib.request.urlretrieve(url, output.gtf + ".gz")
        shell("gunzip {output.gtf}.gz")

rule download_blacklist:
    output: 
        bed = "downloads/encode_blacklist/blacklist.encode.{version}.bed"
    params:
        file = lambda w: {"grch37": "ENCFF001TDO", "grch38": "ENCFF356LFX"}[w.version], 
        provider = "ucsc"
    wildcard_constraints: 
        version = "grch37|grch38"
    shell: 
        op.as_one_line("""
        wget -qO- https://www.encodeproject.org/files/{params.file}/@@download/{params.file}.bed.gz |
        gzip -dc > {output.bed}
        """)

rule download_repeatmasker: 
    output: 
        bed = "downloads/repeatmasker/repeatmasker.{version}.bed"
    params: 
        provider = "ucsc", 
        version = lambda w: {"grch37": "hg19", "grch38": "hg38"}[w.version],
    conda: CONDA_ENVS["bedops"]
    shell: 
        op.as_one_line("""
        wget -qO- http://www.repeatmasker.org/genomes/{params.version}/RepeatMasker-rm405-db20140131/{params.version}.fa.out.gz | 
        gzip -dc | rmsk2bed > {output.bed}
        """)


rule download_dbsnp_vcf:
    output:
        vcf = "downloads/dbsnp-{dbsnp_build}/dbsnp.common_all.{version}.vcf"
    log:
        "downloads/dbsnp-{dbsnp_build}/dbsnp.common_all.{version}.vcf.log"
    params: 
        provider = "ensembl",
        url = lambda w: {"grch37": "GRCh37p13", "grch38": "GRCh38p7"}[w.version]
    conda: CONDA_ENVS["coreutils"]
    shell:
        op.as_one_line("""
        curl -s https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b{wildcards.dbsnp_build}_{params.url}/VCF/00-common_all.vcf.gz 2> {log}
            | 
        gzip -dc 2>> {log}
            | 
        awk 'BEGIN {{FS=OFS="\t"}} $0 !~ /^#/ {{$8="."}} $0 !~ /^##INFO/' > {output.vcf} 2>> {log}
        """)


rule download_af_only_gnomad_vcf:
    output:
        vcf = "downloads/gnomad/af-only-gnomad.{version}.vcf"
    log:
        "downloads/gnomad/af-only-gnomad.{version}.vcf.log"
    params:
        provider = lambda w: {"grch37": "ensembl", "grch38": "ucsc"}[w.version],
        file = lambda w: {
            "grch37": "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf", 
            "grch38": "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
        }[w.version]
    conda: CONDA_ENVS["gsutil"]
    shell:
        op.as_one_line("""
        if [[ {params.file} == *".gz" ]]; then 
            gsutil cp {params.file} - 2> {log}
                |
            gzip -dc 
                |
            awk '{{FS=OFS="\t"}} {{ 
                if ($0 ~ /^##INFO/ && !($0 ~ /ID=AC,/ || $0 ~ /ID=AF,/)) {{ next }} 
                else {{ print }}
            }}'
                > {output.vcf} 2>> {log}; 
        else
            gsutil cp {params.file} - 2> {log}
                |
            awk '{{FS=OFS="\t"}} {{ 
                if ($0 ~ /^##INFO/ && !($0 ~ /ID=AC,/ || $0 ~ /ID=AF,/)) {{ next }} 
                else {{ print }}
            }}'
                > {output.vcf} 2>> {log}; 
        fi
        """)

rule download_mutect2_pon:
    output:
        vcf = "downloads/mutect2/mutect2_pon.{version}.vcf"
    log:
        "downloads/mutect2/mutect2_pon.{version}.vcf.log"
    params:
        provider = lambda w: {"grch37": "ensembl", "grch38": "ucsc"}[w.version],
        file = lambda w: {
            "grch37": "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf", 
            "grch38": "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
        }[w.version]
    conda: CONDA_ENVS["gsutil"]
    shell:
        op.as_one_line("""
        if [[ {params.file} == *".gz" ]]; then 
            gsutil cp {params.file} - 2> {log}
            |
            gzip -dc > {output.vcf}; 
        else
            gsutil cp {params.file} {output.vcf} 2> {log}; 
        fi
        """)

rule download_mutect2_small_exac:
    output:
        vcf = "downloads/mutect2/mutect2_small_exac.{version}.vcf"
    log:
        "downloads/mutect2/mutect2_small_exac.{version}.vcf.log"
    params:
        provider = lambda w: {"grch37": "ensembl", "grch38": "ucsc"}[w.version],
        file = lambda w: {
            "grch37": "gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf", 
            "grch38": "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"
        }[w.version]
    conda: CONDA_ENVS["gsutil"]
    shell:
        op.as_one_line("""
        if [[ {params.file} == *".gz" ]]; then 
            gsutil cp {params.file} - 2> {log}
            |
            gzip -dc > {output.vcf}; 
        else
            gsutil cp {params.file} {output.vcf} 2> {log}; 
        fi
        """)

rule download_liftover_chains:
    output:
        chains = "downloads/chains/{genome_build}/{chain_version}.{version}.over.chain"
    wildcard_constraints:
        chain_version = "hg19ToHg38|hg38ToHg19"
    params:
        provider = lambda w: {"grch37": "ensembl", "grch38": "ucsc"}[w.version],
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19"
    shell:
        op.as_one_line("""
        wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/{params.build}/liftOver/{wildcards.chain_version}.over.chain.gz |
        gzip -dc > {output.chains}
        """)

rule download_sdf: 
    output: 
        sdf = directory("downloads/sdf/{genome_build}/sdf")
    params: 
        build = lambda w: {
            "grch37": "1000g_v37_phase2.sdf", 
            "hs37d5": "1000g_v37_phase2.sdf",
            "hg19": "hg19.sdf", 
            "hg38": "GRCh38.sdf",
            "hg19-reddy": "hg19.sdf"
        }[w.genome_build]
    shell: 
        op.as_one_line("""
        wget -qO {output.sdf}.zip https://s3.amazonaws.com/rtg-datasets/references/{params.build}.zip && 
        unzip -d $(dirname {output.sdf})/{params.build} {output.sdf}.zip &&
        mv $(dirname {output.sdf})/{params.build}/* {output.sdf}
        """)


##### FUNCTIONS #####


def get_matching_download_rules(file):
    ignored_rules = ["download_genome_fasta", "download_sdf"]
    rule_names = [ r for r in dir(rules) if r.startswith("download_")]
    rule_names = [ r for r in rule_names if r not in ignored_rules ]
    rule_list = [ getattr(rules, name) for name in rule_names ]
    matching_rules = []
    for rule_name, r in zip(rule_names, rule_list):
        # Ensure provider is specified
        assert "provider" in r.params._names, (
                f"The `{rule_name}` download rule doesn't have a `provider` param."
            )
        # At least one output file should be returned
        num_matches = []
        for output_file in r.output:
            if rule_name != "download_capspace_bed":  # Workaround since we don't really care about the provider for the capture space
                assert "{version}" in output_file, (
                    f"The `{rule_name}` download rule doesn't have a `{{version}}` "
                    f"wildcard in the output file ('{output_file}')."
                )
            matches = smk.io.glob_wildcards(output_file, [file])
            num_matches.append(len(matches[0]))
        if any(num > 0 for num in num_matches):
            matching_rules.append(r)
    return matching_rules


def hardlink_same_provider(wildcards):

    genome_build = wildcards.genome_build
    suffix = wildcards.suffix
    output_file = f"genomes/{genome_build}/{suffix}"
    raw_download_file = f"downloads/{suffix}"

    version = config["genome_builds"][genome_build]["version"]
    to_provider = config["genome_builds"][genome_build]["provider"]
    from_provider_options = CHROM_MAPPINGS[version][to_provider]
    
    dependencies = []
    matching_rules = get_matching_download_rules(raw_download_file)

    for r in matching_rules:
        # The provider must be among the ones we can convert from
        from_provider = r.params.provider
        # Handle provider as functions
        if callable(from_provider):
            wildcards.version = version
            from_provider = from_provider(wildcards)
        if from_provider not in from_provider_options:
            logger.warning(
                f"The {r.rule} rule can generate the {output_file} file, but the chromosomes "
                f"from the associated provider ({from_provider}) cannot be converted to "
                f"the destination provider ({to_provider}). Make sure the providers in the "
                "download rules are correct and that no chromosome mappings are missing."
            )
            continue
        dependencies.append(r)

    if len(dependencies) == 0:
        msg = f"Could not find rule to generate {output_file}."
        raise Exception(msg)

    if len(dependencies) > 1:
        msg = f"Found conflicting rules to generate {output_file}."
        raise Exception(msg)

    raw_download_root, raw_download_ext = os.path.splitext(raw_download_file)
    target_download_file = f"{raw_download_root}.{to_provider}{raw_download_ext}"
    return target_download_file


def get_cvbio_params(field):

    def get_cvbio_params_custom(wildcards, input, output):

        dependencies = []
        matching_rules = get_matching_download_rules(input.before)
        assert len(matching_rules) == 1, "Not just one matching rule."

        dependency = matching_rules[0]
        version = VERSION_UPPER[wildcards.version]
        from_provider = dependency.params.provider
        if callable(from_provider):
            from_provider = from_provider(wildcards)
        to_provider = wildcards.to_provider
        file_ext = wildcards.ext

        cvbio_params = {
            "mapping": f"{CHROM_MAPPINGS_DIR}/{version}_{from_provider}2{to_provider}.txt",
            "from_provider": from_provider,
            "comment": config["cvbio_config"][file_ext]["comment"],
            "columns": config["cvbio_config"][file_ext]["columns"],
            "skip": config["cvbio_config"][file_ext]["skip"],
            "delimiter": config["cvbio_config"][file_ext]["delimiter"]
        }
        return cvbio_params[field]

    return get_cvbio_params_custom


def get_download_file(file):
    file_root, file_ext = os.path.splitext(file)
    assert file_ext.lstrip(".") in config["cvbio_config"], (
        f"`{file_ext}` files are not yet configured for cvbio. See `config['cvbio_config']` variable."
    )
    def get_download_file_custom(wildcards):
        genome_build = wildcards.genome_build
        version = config["genome_builds"][genome_build]["version"]
        provider = config["genome_builds"][genome_build]["provider"]
        download_file = file.replace("{version}", version).replace("{provider}", provider)
        download_file = "genomes/{genome_build}/" + download_file
        return download_file
    return get_download_file_custom

### CAPTURE SPACE ###
rule download_capspace_bed:
    output:
        capture_bed = "downloads/capture_space/{capture_space}.{genome_build}.bed"
    log:
        "downloads/capture_space/{capture_space}.{genome_build}.bed.log"
    params:
        url = lambda w: config["capture_space"][w.capture_space]["capture_bed_url"],
        provider = lambda w: config["capture_space"][w.capture_space]["provider"]
    shell:
        "curl -L {params.url} > {output.capture_bed} 2> {log}"

rule add_remove_chr_prefix_bed:
    input:
        capture_bed = rules.download_capspace_bed.output.capture_bed
    output:
        converted_bed = "downloads/capture_space/{capture_space}.{genome_build}.{chr_status}.bed"
    run:
        # Converts the specified BED file and adds/removes chr prefixes
        if wildcards.chr_status == "chr":  # i.e. we need to add a chr prefix
            add_chr = True
        else:
            add_chr = False

        # Process the BED file
        with open(input.capture_bed) as f, open(output.converted_bed, "w") as o:
            i = 0
            for line in f:
                i += 1
                # Make sure that this BED entry is chr-prefixed or not
                if add_chr and line.startswith("chr"):
                    # We were asked to add a chr prefix, but one already exists
                    raise AttributeError("I was asked to add a \'chr\' prefix to \'%s\', but that BED is already \'chr\' prefixed on line %s" % (wildcards.input.capture_bed, line))
                if not add_chr and not line.startswith("chr"):
                    # We were asked to remove a chr prefix, but there isn't one
                    raise AttributeError("I was asked to remove a \'chr\' prefix \'%s\', but it isn't chr prefixed on line %s" % (wildcards.input.capture_bed, line))

                if add_chr:
                    line = "chr" + line
                else:
                    line = line.replace("chr", "")

                o.write(line)


##### SHARED #####


rule hardlink_download:
    input: hardlink_same_provider
    output: "genomes/{genome_build}/downloads/{suffix}"
    shell: "ln -f {input} {output}"


rule update_contig_names:
    input:
        before = "downloads/{parent_dir}/{prefix}.{version}.{ext}"
    output:
        after = "downloads/{parent_dir}/{prefix}.{version}.{to_provider}.{ext}"
    log:
        "downloads/{parent_dir}/{prefix}.{version}.{to_provider}.{ext}.log"
    params:
        mapping = get_cvbio_params("mapping"),
        from_provider = get_cvbio_params("from_provider"),
        comment = get_cvbio_params("comment"),
        columns = get_cvbio_params("columns"),
        skip = get_cvbio_params("skip"),
        delimiter = get_cvbio_params("delimiter")
    wildcard_constraints:
        ext = "|".join(config["cvbio_config"].keys()),
        to_provider = "|".join(TO_PROVIDERS)
    conda: CONDA_ENVS["cvbio"]
    shell:
        op.as_one_line("""
        if [[ '{wildcards.to_provider}' == '{params[from_provider]}' ]]; then
            ln -f {input.before} {output.after};
        else
            cvbio UpdateContigNames --in {input.before} --out {output.after}
            --mapping {params.mapping} --comment-chars '{params.comment}'
            --columns {params.columns} --skip-missing {params.skip}
            --delimiter '{params.delimiter}' > {log} 2>&1; 
        fi
        """)
