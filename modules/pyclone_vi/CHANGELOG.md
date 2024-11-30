# Changelog

All notable changes to the `pyclone_vi` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-02-22

This release was authored by Laura Hilton.
See details about [PyClone](https://github.com/Roth-Lab/pyclone) and [PyClone-VI](https://github.com/Roth-Lab/pyclone-vi).
The PhyClone GitHub repository is currently private but the tool is available on Conda. Several scripts were adapted from those written by Andrew Roth.

This module runs PyClone (either the original meant for capture data) or PhyClone-VI (newer, meant for genome data), followed by PhyClone to refine clusters and assign a phylogeny. Runs on multi-time point patients, and requires a time_point column in the input metadata. The original PyClone does include some plots in its outputs. PhyClone and PyClone-VI do not generate any plots currently but this will be added as a feature in the near future.

PyClone-VI seems to offer superior performance using only coding and aSHM mutations, rather than all mutations in the genome. This is togglable in the config with the `subset_maf` option.
