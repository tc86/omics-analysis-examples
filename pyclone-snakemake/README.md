# PyClone Snakemake Workflow

A minimal Snakemake workflow to demonstrate running PyClone on somatic variant data.

**Pipeline structure (planned):**
- Input: VCFs + copy number table
- Preprocessing & filtering
- PyClone run
- Summarization & visualization

**Features:**
- Modular rules
- Configurable with `config.yaml`
- Small toy dataset for quick testing

**Dependencies:** Snakemake, PyClone, Python ≥3.8

> Intended as a lightweight showcase of workflow design principles, not a production-ready pipeline.
