# krak

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/krak/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/krak/badges/version.svg)](http://bioconda.github.io/recipes/krak/README.html)
[![Build Status](https://github.com/clintval/krak/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/clintval/krak/actions/workflows/rust.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/krak/badge.svg?branch=main)](https://coveralls.io/github/clintval/krak?branch=main)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

An addicting set of Kraken-enhancing tools.

![Monkey Pod Tree](.github/img/cover.jpg)

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    -c conda-forge -c bioconda \
    krak --help
```

## Introduction

This project provides tools for integrating [Kraken](https://github.com/DerrickWood/kraken) and [Kraken2](https://github.com/DerrickWood/kraken2) taxonomic classifications into FASTX (FASTQ or FASTA), SAM, BAM, and CRAM workflows.
Kraken creates _k_-mer-based taxonomic labels from genomic sequences such as single-end or paired-end sequencing reads.
`krak` bridges the gap between genomic sequences stored in alignment files and their taxonomic classifications so that you can filter, rescue, or otherwise act on the classification of each genomic sequence.

All subcommands support both Kraken v1 and Kraken v2 output formats.

## Quick Start

Classify, annotate, and filter a queryname-grouped BAM:

```bash
krak annotate \
    -i input.bam \
    -d /kraken-db \
    -a <(krak prep input.bam | kraken2 --db /kraken-db --output - -) \
| krak filter -t 9606 -o output.bam
```

## Features

- `krak prep`: for preparing FASTA/FASTQ/SAM/BAM/CRAM for Kraken input
- `krak annotate`: annotates SAM/BAM/CRAM records with Kraken classifications
- `krak n2ref`: converts aligned Ns to reference bases
- `krak filter`: filters a FASTA/FASTQ/SAM/BAM/CRAM based on Kraken output
- `krak report2tsv`: converts a Kraken report to a tab-separated text file

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
