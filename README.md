# GConstraint VEP Plugin

An Ensembl VEP plugin that adds gene constraint annotations from gnomAD version 4.1.1.

## Description

This plugin annotates transcript-level variants with gnomAD v4.1.1 gene constraint metrics. It uses a tabix-indexed constraint metrics TSV file to retrieve values for the selected transcript.

## Requirements

- Ensembl VEP
- `tabix` installed and available in `PATH`
- gnomAD v4.1.1 constraint metrics TSV file, tabix-indexed

## Installation

1. Place `GConstraint.pm` in your VEP plugins directory, for example:
   ```bash
   mv GConstraint.pm ~/.vep/Plugins/
   ```
2. Make sure the constraint metrics file is downloaded, sorted, and tabix-indexed.

## Input file preparation

The gnomAD constraint metrics file can be downloaded from:

- https://gnomad.broadinstitute.org/downloads#v4

To prepare the file for tabix, use:

```bash
zcat gnomad.v4.1.1.constraint_metrics.tsv.bgz | (head -n 1 && tail -n +2 | sort -t$'\t' -k 9,9 -k 10,10n ) > gnomad.v4.1.1.constraint_metrics_sorted.tsv
sed '1s/.*/#&/' gnomad.v4.1.1.constraint_metrics_sorted.tsv > gnomad.v4.1.1.constraint_metrics_final.tsv
bgzip gnomad.v4.1.1.constraint_metrics_final.tsv
tabix -f -s 9 -b 10 -e 11 gnomad.v4.1.1.constraint_metrics_final.tsv.gz
```

## Usage

Annotate the default columns:

```bash
vep -i variations.vcf --plugin GConstraint,file=/path/to/gnomad_constraint_final.tsv.gz
```

Annotate a custom set of columns:

```bash
vep -i variations.vcf --plugin GConstraint,file=/path/to/gnomad_constraint_final.tsv.gz,cols=syn.z_score,mis.z_score,lof.pLI
```

## Supported columns

Default columns returned by the plugin:

- `syn.z_score`
- `syn.oe_ci.lower`
- `syn.oe_ci.upper`
- `mis.z_score`
- `mis.oe_ci.lower`
- `mis.oe_ci.upper`
- `lof.z_score`
- `lof.oe_ci.lower`
- `lof.oe_ci.upper`
- `lof.pLI`

Additional supported columns:

- `gene_flags`
- `constraint_flags`

## License

Apache License 2.0

## Contact

For questions, contact: romualdofarm@gmail.com
