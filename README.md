# CLSA GWAS for T2D Complications

Workspace for curating CLSA T2D complication phenotypes and running GWAS (Regenie) using CLSA genotype data.

## Inputs used

Admin instructions (this repo):

- `0.admin/1.email.pdf`
- `0.admin/T2D-complications_meta-analysis_AnalysisPlan_Phenodefinitions_22122025.docx`

CLSA phenotype data:

- Baseline: `/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_Baseline`
- FUP1: `/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_FUP1`
- FUP2: `/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_FUP2`

CLSA genotype / sample QC data:

- `/home/csu/projects/rrg-vmooser/CERC_Private/Geno/CLSA/Genomics3_clsa`

Reference project used for layout and scripts:

- `/home/csu/scratch/38.clsa-go3`

## Current workflow

1. Curate phenotypes and write Regenie input files

```bash
python3 01.phenotype/01.curate_t2d_complications.py
```

2. Run (once) imputed BGEN -> PGEN conversion/filter prep

```bash
cd 02.regenie/02.filter_impute
sbatch 01.filter.sh
```

3. Submit Regenie step 1 for all ready phenotypes listed in `01.phenotype/output/regenie_targets.txt`

```bash
cd 02.regenie/01.step1
./run_01.regenie_step1.sh
```

4. Submit Regenie step 2 (updated stricter settings: `--firth`, `--minMAC 20`, `--minINFO 0.8` when INFO exists)

```bash
cd 02.regenie/03.step2
./run_01.regenie_step2.sh
```

## Notes

- The current CLSA mapping supports a self-report-based T2D-complication pipeline with some known phenotype limitations (see `docs/phenotype_definition_summary.md`).
- Neurological diabetic complication (`neuropathy`) is not currently mapped from the reviewed CLSA variables and is marked unavailable in the manifest.
- The phenotype script defaults to the `EUR` ancestry cluster (`pca.cluster.id == 4`) to match the prior OA workflow; change this via CLI args for other ancestry groups.
