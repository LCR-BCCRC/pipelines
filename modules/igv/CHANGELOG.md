# Changelog

All notable changes to the `igv` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-01-10

This release was authored by Manuela Cruz.

- This module requires four file types:
    * Regions file containing desired regions to be snapshot in BED format, MAF format, mutation_id format in which mutations are in "{chromosome}:{start_position}:{end_position}" format or OncodriveCLUSTL / HotMAPS results files. HotMAPS and OncodriveCLUSTL results are preformatted using scripts that are executed in the last steps of their lcr-modules.
    * BAM and 
    * BAI files for each sample to generate IGV screenshots
    * MAF files for each sample to determine which variants are included in desired regions to be snapshot and extract the corresponding Chromosome, Start_Position and Hugo_Symbol values

- BAM and BAI file locations are sourced based on the sample_id and the corresponding genome_build and seq_type metadata values set in CFG["samples"]

- As multiple genome build values other than `grch37` and `hg38` can exist in the metadata genome_build column, add the genome build values to their respective `grch37` or `hg38` lists in the config at ["options"]["genome_map"] to ensure the correct MAF paths are resolved for each sample.

- The regions files must be entered in the config under the correct key denoting the genome build of the regions file. This is required to perform proper liftover if necessary, which allows the workflow to correctly filter sample MAFs that are in builds opposite what is provided in the regions file.

- IGV batch scripts are created for each individual variant of interest. This is a checkpoint rule as it depends on the file contents of the filtered MAF files.

Creation of IGV snapshots:

- The individual IGV batch scripts are appended to a single large "merged" batch script for each sample_id and an empty "dispatched" file is created for each indiviual batch script. Importantly, this rule only runs if the dispatched file has not been created in order to prevent variants that have already undergone IGV snapshots to be appended and rerun again.

- IGV is then run on each sample's merged batch script. This is also a checkpoint rule as the specific snapshots that will be created depend on the contents of the merged batch scripts.

- IGV snapshot presets can be defined in the config at ["options"]["generate_batch_script"]["igv_options"] parameters. Two presets have already been defined, `default` and `pairs`. To specify which presets you want to create images for in the analysis, define them in list format in the config at ["options"]["igv_presets"]. More information on preset options here: https://github.com/igvteam/igv/wiki/Batch-commands

Blank or truncated snapshots:

- Truncated or blank snapshots can occur during the IGV run. The quality control rule checks for blank snapshots based on the kurtosis and skewness values of the snapshot, and checks for truncated snapshots based on the height and width. Multiple attempts are performed in order to resolve the affected snapshots, and if they remain unresolved the quality control rule will fail. Note that the thresholds for blank and truncated snapshots have been determined based on IGV image dimensions 1920x1080x24, and modifying image dimensions may result in blank or truncated snapshots being missed during quality control.

- Increasing the milliseconds set in CFG["options"]["generate_batch_script"]["sleep_timer] typically reduces the amount of snapshots with issues but increases workflow run time.

- A summary file of failed or suspicious snapshots is created for each run in `99-outputs/snapshot_summaries/quality_control/`. Failed snapshots match specific height thresholds known to belong to blank or truncated snapshots. Suspicious snapshots are snapshots that reached the limit of re-run attempts during quality control. These may be blank, truncated, or they may be regular snapshots that display only a few reads and went through the quality control process due to their short heights.

Estimating snapshots:

- To estimate the number of IGV snapshots that will be created, the config parameter ["estimate_only"] can be set to True. Summary files are created in `99-outputs/snapshot_summaries/estimates/` based on the individual batch scripts that have been created and do not have pre-existing "dispatch" files. 


