# LME MaPSY

## Table of Contents

1. [Overview](#overview)
2. [Project Structure](#project-structure)
3. [Dependencies](#dependencies)
4. [Usage](#usage)
    - [First Steps](#first-steps)
    - [Step-by-Step Analysis](#step-by-step-analysis)
5. [Classes and Methods](#classes-and-methods)
    - [ReadsMeta](#readsmeta)
    - [IdentifyAssociations](#identifyassociations)
    - [ResolveAssociations](#resolveassociations)
    - [ProcessIncData](#processincdata)
    - [QuantifyInclusion](#quantifyinclusion)
6. [License](#license)
7. [Authors](#authors)

## Overview

This repository contains the code used to process and analyse massively parallel splicing assay data for the study:

Bonnal, S., Bajew, S., Martinez-Corral, R., et al. **Core splicing architecture and early spliceosomal recognition determine microexon sensitivity to SRRM3/4.** _Nature Structural & Molecular Biology_ (2025). [doi.org/10.1038/s41594-025-01634-1](https://doi.org/10.1038/s41594-025-01634-1)

## Project Structure

```
.
├── lme-mapsy
│   ├── __init__.py
│   ├── reads_meta.py
│   ├── identify_associations.py
│   ├── resolve_associations.py
│   ├── process_inc_data.py
│   └── quantify_inclusion.py
├── LICENSE
└── README.md
```

## Dependencies

```
python >= 3.9
pandas >= 1.3.0
regex >= 2021.8.3
biopython >= 1.79
numpy >= 1.21.0
gzip
```

## Usage

### First Steps

### Step-by-Step Analysis

1. **Parsing Reads**

2. **Identify Barcode-Variant Associations**

3. **Resolve Associations**

4. **Quantify Inclusion**

5. **Process Inclusion Data**

## Classes and Methods

### ReadsMeta

Class for handling sequencing reads and metadata.

#### Methods:
- `parse_reads(reads: str, reverse: bool = False) -> dict`
- `parse_barcodes(reads: dict, regex_str: str) -> tuple[dict, pandas.DataFrame]`
- `filter_reads(valid_barcodes: dict, reads: dict) -> dict`
- `parse_variants_metadata(metadata: str) -> tuple[dict, dict]`

<hr>

### IdentifyAssociations

Child class for parsing sequencing reads from LME MaPSy libraries.

#### Methods:
- `match_variants(filtered_reads: dict, oligo_sequences: dict, variant_seq_lengths: dict) -> pandas.DataFrame`

<hr>

### ResolveAssociations

Class to handle the resolution of barcode-variant associations.

#### Methods:
- `__init__(input_dir: str, max_mismatch_pct: int, min_nreads: int, nreads_second_highest: int, pct_second_best: int, min_rep: int)`
- `call_mismatches(data: pandas.DataFrame, pct_mms: int) -> tuple[pd.DataFrame, pd.DataFrame]`
- `__read_input_data(associations_dir: str) -> pd.DataFrame`
- `__filter_associations(accepted_data: pd.DataFrame, min_nreads: int) -> tuple[pd.DataFrame, pd.DataFrame]`
- `__resolve_associations(correct_data: pd.DataFrame, misassignment_data: pd.DataFrame, nreads_second_highest: int, pct_second_best: int, min_rep: int) -> tuple[pd.DataFrame, pd.DataFrame]`

<hr>

### QuantifyInclusion

Class for quantifying the inclusion levels and filtering barcode outliers from LME MaPSy libraries.

#### Methods:
- `__init__(meta: ReadsMeta, associations: str, condition: str, correct_aberrant: bool, replicates: bool)`
- `__process_associations(data_dir: str) -> dict`
- `__parse_variants_metadata(metadata: str) -> tuple[dict, dict]`
- `__quantify_inclusion(reads: dict, exon_regex: str, associations: dict, trusted_barcodes: dict, variants_seq: dict, variants: dict, condition: str) -> tuple[pandas.DataFrame, pandas.DataFrame]`
- `__correct_inclusion(inclusion_table: pandas.DataFrame, summary_table: pandas.DataFrame, associations: dict) -> tuple[pandas.DataFrame, pandas.DataFrame]`

<hr>

### ProcessIncData

Contains methods to filter out barcode outliers and create a meta inclusion table.

#### Methods:
- `__init__(barcode_inclusion_table: pandas.DataFrame, min_nreads: int, max_proportion: float, min_barcodes: int, min_samples: int, factor: float, variants_info: str)`
- `__select_min_barcodes(row: pandas.Series, min_bc: int, min_cond: int) -> str`
- `__filter_outliers_row(row: pandas.Series, factor: float) -> str`
- `__filter_outliers(corr_inclusion_table: pandas.DataFrame, nreads: int, max_prop: float, min_barcodes: int, min_conditions: int, fc: float) -> tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]`
- `__prepare_meta_inclusion(inc_table: pandas.DataFrame, info: str) -> pandas.DataFrame`


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

- **[Simon Bajew, Ph.D.](https://github.com/simon-bt)**
- **[Manuel Irimia, Ph.D.](https://github.com/mirimia)**
