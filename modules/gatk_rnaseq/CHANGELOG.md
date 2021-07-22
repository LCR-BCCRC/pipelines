# Changelog

All notable changes to the `gatk_rnaseq` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-04-23

This release was authored by Jasper Wong.

The gatk_rnaseq module follows GATK's variant-calling for RNA-seq pipeline, which involves:
1) SplitNCigarReads - Splitting reads into exon segments (getting rid of N's but keeping grouping information) and hard-clipping overhangs (that cross into intron regions). By default, this also reassigns mapping quality for good alignments (to match DNA alignments)
2) BaseRecalibrator (and applyBQSR) - creates a table that detects and corrects for patterns of systematic errors in the base quality scores. Applies this recalibration step on the BAMs.
3) Variant-calling - I split up the variant calling into one per chrom (24 threads), and run HaplotypeCaller. Default is phred-scaled confidence threshold of 20 for calling variants. Filters out db-snp regions.
4) Variant-filtration - Several hard-filters applied here. Flag any variants that are FS > 30, QD < 2, DP < 5.
5) Select only passed variants.

