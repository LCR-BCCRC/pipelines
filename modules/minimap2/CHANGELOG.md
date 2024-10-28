# Changelog

All notable changes to the `minimap2` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-03-14

This release was authored by Haya Shaalan.

<!-- TODO: Explain each important module design decision below. -->

- This module can take paired and unpaired fastq files.
- Performs short and long read alignment based on {seq_type}. Parameters are switched and can be configured through the config.
- Uses the utils module and writes the outputs to 99-outputs.
- Final output is a bam file with naming format: bam/{seq_type}--{genome_build}/{sample_id}.bam.
