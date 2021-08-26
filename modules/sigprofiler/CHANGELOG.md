# Changelog

All notable changes to the `sigprofiler` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-08-16

This release was authored by Prasath Pararajalingam.

- NMF occurs in two stages (estimate and extract). Estimate runs NMF on a range of ranks with low numer of replications to find the optimal rank. Extract uses the optimal rank to run NMF on optimal rank +/- 2 using a high number of replications. This should provide a speed-up for matrices containing many tumour samples.

