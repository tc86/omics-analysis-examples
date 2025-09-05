# PyClone Snakemake / Input Builder

Utility to generate a **PyClone-VI input TSV** from per-sample VCFs and a CNV table.

## Usage
```bash
python make_pyclone_input.py   --vcf-dir vcfs/   --cnv-info cnv_info.tsv   --output pyclone_input.tsv
```

## Notes
- Requires: Python ≥3.8, pandas
- Columns in CNV table: `sample_id chrom start end major_cn minor_cn normal_cn [tumour_content]`
- Chrom labels normalized (chr1 → 1)

## Pyclone-VI information
For installation guide and input/output information please visit the [Pyclone-VI](https://github.com/Roth-Lab/pyclone-vi) github page. 