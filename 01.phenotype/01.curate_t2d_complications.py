#!/usr/bin/env python3
"""Curate CLSA T2D complication phenotypes for GWAS (Regenie-ready outputs).

This script harmonizes CLSA baseline/FUP1/FUP2 self-reported phenotype data,
joins sample QC/genetics metadata (CLSA v3 SQC), and writes binary phenotype,
covariate, and keep files for T2D complication GWAS.

Design choices (current implementation):
- T2D ascertainment uses self-report diabetes + type II evidence (baseline/FUP1).
  No ICD linkage is used; broader ascertainment is possible if linked data become
  available.
- Complications are self-report based (plus dialysis as renal severity proxy).
- Ophthalmic phenotype uses explicit diabetic retinopathy variables only (FUP1/FUP2);
  broader ophthalmic outcomes (cataract, glaucoma, macular edema) are not captured
  in the CLSA variables examined.
- Renal phenotype is a self-report kidney disease/failure + dialysis proxy, not strict
  DKD/CKD based on lab criteria (albuminuria/eGFR).
- Cardiovascular phenotype includes hypertension (CCC_HBP) in addition to heart
  disease, angina, MI, CAB/revascularization, and PVD.
- Neurological diabetic complication phenotype is not available in CLSA variables
  examined and is therefore flagged unavailable in the manifest. MICRO is partial
  (no neurological component).
- Sensitivity phenotypes require complication age >= T2D age + 2 years, when both ages
  are available. Cases with unknown timing remain eligible only for main analyses.
- Main analyses exclude cases where known complication age precedes known T2D age;
  this is stricter than a literal "no time limit" interpretation but avoids
  biologically implausible reverse-timing cases.
- Covariates: age, t2d_duration, sex (omitted in sex-stratified runs), PCs, center
  dummies, and batch dummies. One reference level is dropped for center and batch
  to avoid collinearity with the intercept.

Outputs (default):
  01.phenotype/output/
    - phenotype_manifest.tsv
    - regenie_targets.txt
    - *_pheno.tsv
    - *_covar.tsv
    - *_FIDIID_noheader.tsv
    - participant_harmonized.tsv.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# Defaults in this environment
DEFAULT_BASELINE_CSV = "/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_Baseline/2109005_McGill_DTaliun_Baseline_CoPv7.csv"
DEFAULT_FUP1_CSV = "/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_FUP1/2109005_McGill_DTaliun_FUP1_CoPv3_2.csv"
DEFAULT_FUP2_CSV = "/home/csu/projects/rrg-vmooser/CERC_Private/Pheno/CLSA/2109005_McGill_DTaliun_FUP2/2109005_McGill_DTaliun_FUP2_CoPv3.csv"
DEFAULT_SQC = "/home/csu/projects/rrg-vmooser/CERC_Private/Geno/CLSA/Genomics3_clsa/clsa_sqc_v3.txt"
DEFAULT_OUTDIR = os.path.join(os.path.dirname(__file__), "output")


BASELINE_COLS = [
    "entity_id",
    "ADM_GWAS3_COM",
    "SEX_ASK_COM",
    "AGE_NMBR_COM",
    "GEOSTRATA_COM",
    "DIA_DIAB_COM",
    "DIA_TYPE_COM",
    "DIA_AGE_NB_COM",
    "DIA_MED_COM",
    "DIA_MEDCUR_COM",
    "DIA_MEDAGE_NB_COM",
    "DIA_PRGDIA_COM",
    "CCC_HEART_COM",
    "CCC_PVD_COM",
    "CCC_HBP_COM",
    "CCC_KIDN_COM",
    "CCC_DITYP_COM",
    "CCC_ANGI_COM",
    "CCC_CVA_COM",
    "CCC_AMI_COM",
    "CCC_TIA_COM",
    "CCC_CVAFX_COM",
    "IHD_CAB_COM",
    "IHD_ANGIAGE_NB_COM",
    "IHD_AMIAGE_NB_COM",
    "STR_CVAAGE_NB_COM",
    "STR_TIAAGE_NB_COM",
    "BLD_HBA1c_COM",
    "BLD_CREAT_COM",
    "BLD_ALB_COM",
]

FUP1_COLS = [
    "entity_id",
    "AGE_NMBR_COF1",
    "DIA_DIAB_COF1",
    "DIA_TYPE_COF1",
    "DIA_MED_COF1",
    "DIA_DIABRT_COF1",
    "CCC_DIAB_DRAGE_NB_COF1",
    "CCC_HEART_COF1",
    "CCC_HEARTAGE_NB_COF1",
    "CCC_PVD_COF1",
    "CCC_PVDAGE_NB_COF1",
    "CCC_HBP_COF1",
    "CCC_KIDN_COF1",
    "CCC_KIDNAGE_NB_COF1",
    "CCC_DITYP_COF1",
    "CCC_ANGI_COF1",
    "CCC_CVA_COF1",
    "CCC_AMI_COF1",
    "CCC_TIA_COF1",
    "CCC_CVAFX_COF1",
    "IHD_CAB_COF1",
    "IHD_ANGIAGE_NB_COF1",
    "IHD_AMIAGE_NB_COF1",
    "STR_CVAAGE_NB_COF1",
    "STR_TIAAGE_NB_COF1",
]

FUP2_COLS = [
    "entity_id",
    "SEX_ASK_COM",
    "AGE_NMBR_COF2",
    "DIA_DIAB_COF2",
    "DIA_AGE_NB_COF2",
    "DIA_MED_COF2",
    "DIA_DIABRT_COF2",
    "DIA_DIABRTAGE_NB_COF2",
    "CCC_HEART_COF2",
    "CCC_HEARTAGE_NB_COF2",
    "CCC_KIDN_COF2",
    "CCC_KIDNAGE_NB_COF2",
    "CCC_DITYP_COF2",
    "CCC_ANGI_COF2",
    "CCC_CVA_COF2",
    "CCC_TIA_COF2",
    "CCC_AMI_COF2",
    "CCC_CVAFX_COF2",
    "IHD_CAB_COF2",
    "IHD_ANGIAGE_NB_COF2",
    "IHD_AMIAGE_NB_COF2",
    "STR_CVAAGE_NB_COF2",
    "STR_TIAAGE_NB_COF2",
]

# SQC columns: file is whitespace-delimited, not tab-delimited.
SQC_COLS = [
    "ADM_GWAS_COM",
    "batch",
    "chromosomal.sex",
    "pca.cluster.id",
    "in.kinship",
    "in.hetmiss",
    "ePC1",
    "ePC2",
    "ePC3",
    "ePC4",
    "ePC5",
    "ePC6",
    "ePC7",
    "ePC8",
    "ePC9",
    "ePC10",
]

PC_COLS = [f"PC{i}" for i in range(1, 11)]

MISSING_STRINGS = {"", "NA", "NaN", "nan"}


@dataclass
class Participant:
    entity_id: str
    gwasid: str
    sex_mf: str
    age_baseline: Optional[float]
    center: str
    sqc: Dict[str, str]
    baseline: Dict[str, str]
    fup1: Dict[str, str]
    fup2: Dict[str, str]
    sex_regenie: Optional[int] = None  # female=0, male=1
    t2d_flag: bool = False
    t1d_flag: bool = False
    t2d_age: Optional[float] = None
    last_obs_age: Optional[float] = None
    t2d_duration: Optional[float] = None
    control_eligible: bool = False
    # complication flags/ages
    ophth_yes: bool = False
    ophth_age: Optional[float] = None
    renal_yes: bool = False
    renal_age: Optional[float] = None
    cardio_yes: bool = False
    cardio_age: Optional[float] = None
    cerebro_yes: bool = False
    cerebro_age: Optional[float] = None
    neuro_available: bool = False


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline-csv", default=DEFAULT_BASELINE_CSV)
    p.add_argument("--fup1-csv", default=DEFAULT_FUP1_CSV)
    p.add_argument("--fup2-csv", default=DEFAULT_FUP2_CSV)
    p.add_argument("--sqc", default=DEFAULT_SQC)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--ancestry-label", default="EUR")
    p.add_argument("--pca-cluster-id", type=int, default=4)
    p.add_argument("--min-control-t2d-years", type=float, default=2.0)
    p.add_argument("--min-case-gap-years", type=float, default=2.0)
    p.add_argument("--min-cases-for-target", type=int, default=20)
    p.add_argument(
        "--allow-t2d-without-type2",
        action="store_true",
        help="Include diabetes-yes participants without explicit type II evidence if they have diabetes meds and no type I evidence.",
    )
    return p.parse_args()


def load_csv_selected(path: str, columns: Sequence[str], key: str) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        missing_cols = [c for c in columns if c not in reader.fieldnames]
        if missing_cols:
            raise ValueError(f"{path}: missing columns: {missing_cols}")
        for row in reader:
            k = row.get(key, "")
            if not k:
                continue
            out[k] = {c: row.get(c, "") for c in columns}
    return out


def load_sqc(path: str) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().strip().split()
        idx = {name: i for i, name in enumerate(header)}
        missing_cols = [c for c in SQC_COLS if c not in idx]
        if missing_cols:
            raise ValueError(f"{path}: missing SQC columns: {missing_cols}")
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < len(header):
                continue
            row = {c: parts[idx[c]] for c in SQC_COLS}
            gwasid = row["ADM_GWAS_COM"]
            out[gwasid] = row
    return out


def is_missing(v: Optional[str]) -> bool:
    if v is None:
        return True
    return str(v).strip() in MISSING_STRINGS


def to_num(v: Optional[str]) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip()
    if s in MISSING_STRINGS:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def valid_age(v: Optional[str]) -> Optional[float]:
    x = to_num(v)
    if x is None:
        return None
    if x in {-99999, -88888, -88880, -8888, -8, 9998, 9999}:
        return None
    # Some CLSA age fields are integer years. Keep a generous plausible range.
    if x < 0 or x > 130:
        return None
    return x


def valid_lab(v: Optional[str]) -> Optional[float]:
    x = to_num(v)
    if x is None:
        return None
    if x in {-99999, -88888, -88880, -8888, -8888.0, -8}:
        return None
    return x


def parse_int_str(v: Optional[str]) -> Optional[int]:
    x = to_num(v)
    if x is None:
        return None
    if not math.isfinite(x):
        return None
    return int(x)


def yn_code(v: Optional[str], yes_codes: Iterable[str] = ("1",), no_codes: Iterable[str] = ("2",)) -> Optional[bool]:
    s = (v or "").strip()
    if s in MISSING_STRINGS:
        return None
    if s in set(yes_codes):
        return True
    if s in set(no_codes):
        return False
    return None


def min_age(*vals: Optional[float]) -> Optional[float]:
    xs = [x for x in vals if x is not None]
    return min(xs) if xs else None


def max_age(*vals: Optional[float]) -> Optional[float]:
    xs = [x for x in vals if x is not None]
    return max(xs) if xs else None


def any_true(items: Iterable[Optional[bool]]) -> bool:
    return any(x is True for x in items)


def all_covar_present(p: Participant) -> bool:
    if not p.gwasid or p.sex_regenie is None or p.age_baseline is None or not p.center:
        return False
    required = ["batch", "chromosomal.sex", "pca.cluster.id", "in.kinship", "in.hetmiss"] + [f"ePC{i}" for i in range(1, 11)]
    for c in required:
        if is_missing(p.sqc.get(c, "")):
            return False
        if p.sqc.get(c, "") in {"-99999", "-88888", "-88880"}:
            return False
    return True


def compute_t2d_flags(p: Participant, allow_without_type2: bool) -> None:
    b, f1, f2 = p.baseline, p.fup1, p.fup2

    type_vals = [b.get("DIA_TYPE_COM", ""), f1.get("DIA_TYPE_COF1", "")]
    diabetes_vals = [b.get("DIA_DIAB_COM", ""), f1.get("DIA_DIAB_COF1", ""), f2.get("DIA_DIAB_COF2", "")]
    med_vals = [b.get("DIA_MED_COM", ""), f1.get("DIA_MED_COF1", ""), f2.get("DIA_MED_COF2", "")]

    p.t1d_flag = any((v or "").strip() == "1" for v in type_vals)
    type2_any = any((v or "").strip() == "2" for v in type_vals)
    diabetes_yes_any = any_true(yn_code(v) for v in diabetes_vals)
    med_yes_any = any_true(yn_code(v) for v in med_vals)
    hba1c = valid_lab(b.get("BLD_HBA1c_COM", ""))
    hba1c_dm = hba1c is not None and hba1c >= 6.5

    # Strict default: explicit Type II evidence in baseline/FUP1.
    # Optional lenient extension: diabetes+meds and no T1D, or baseline WHO HbA1c criterion.
    if type2_any:
        p.t2d_flag = True and not p.t1d_flag
    elif allow_without_type2 and not p.t1d_flag and (diabetes_yes_any and (med_yes_any or hba1c_dm)):
        p.t2d_flag = True
    else:
        p.t2d_flag = False

    p.t2d_age = min_age(
        valid_age(b.get("DIA_AGE_NB_COM", "")),
        valid_age(f2.get("DIA_AGE_NB_COF2", "")),
    )
    p.last_obs_age = max_age(
        valid_age(b.get("AGE_NMBR_COM", "")),
        valid_age(f1.get("AGE_NMBR_COF1", "")),
        valid_age(f2.get("AGE_NMBR_COF2", "")),
    )
    if p.t2d_age is not None and p.last_obs_age is not None:
        p.t2d_duration = p.last_obs_age - p.t2d_age


def compute_complications(p: Participant) -> None:
    b, f1, f2 = p.baseline, p.fup1, p.fup2

    # Ophthalmic: CLSA explicit diabetic retinopathy (FUP1/FUP2).
    p.ophth_yes = any_true([
        yn_code(f1.get("DIA_DIABRT_COF1", "")),
        yn_code(f2.get("DIA_DIABRT_COF2", "")),
    ])
    p.ophth_age = min_age(
        valid_age(f1.get("CCC_DIAB_DRAGE_NB_COF1", "")),
        valid_age(f2.get("DIA_DIABRTAGE_NB_COF2", "")),
    )

    # Renal: kidney disease/failure and/or dialysis (broad, T2D-restricted at case/control stage).
    dialysis_yes = any_true([
        yn_code(b.get("CCC_DITYP_COM", ""), yes_codes=("1", "2"), no_codes=("3",)),
        yn_code(f1.get("CCC_DITYP_COF1", ""), yes_codes=("1", "2"), no_codes=("3",)),
        yn_code(f2.get("CCC_DITYP_COF2", ""), yes_codes=("1", "2"), no_codes=("3",)),
    ])
    kidney_yes = any_true([
        yn_code(b.get("CCC_KIDN_COM", "")),
        yn_code(f1.get("CCC_KIDN_COF1", "")),
        yn_code(f2.get("CCC_KIDN_COF2", "")),
    ])
    p.renal_yes = kidney_yes or dialysis_yes
    p.renal_age = min_age(
        valid_age(f1.get("CCC_KIDNAGE_NB_COF1", "")),
        valid_age(f2.get("CCC_KIDNAGE_NB_COF2", "")),
    )

    # Cardiovascular: heart disease, angina, MI, CAB/revascularization, PAD/PVD, hypertension.
    heart_yes = any_true([
        yn_code(b.get("CCC_HEART_COM", "")),
        yn_code(f1.get("CCC_HEART_COF1", "")),
        yn_code(f2.get("CCC_HEART_COF2", ""), yes_codes=("1", "11"), no_codes=("2",)),
    ])
    angi_yes = any_true([
        yn_code(b.get("CCC_ANGI_COM", "")),
        yn_code(f1.get("CCC_ANGI_COF1", "")),
        yn_code(f2.get("CCC_ANGI_COF2", "")),
    ])
    ami_yes = any_true([
        yn_code(b.get("CCC_AMI_COM", "")),
        yn_code(f1.get("CCC_AMI_COF1", "")),
        yn_code(f2.get("CCC_AMI_COF2", "")),
    ])
    cab_yes = any_true([
        yn_code(b.get("IHD_CAB_COM", "")),
        yn_code(f1.get("IHD_CAB_COF1", "")),
        yn_code(f2.get("IHD_CAB_COF2", "")),
    ])
    pvd_yes = any_true([
        yn_code(b.get("CCC_PVD_COM", "")),
        yn_code(f1.get("CCC_PVD_COF1", "")),
    ])
    hbp_yes = any_true([
        yn_code(b.get("CCC_HBP_COM", "")),
        yn_code(f1.get("CCC_HBP_COF1", "")),
    ])
    p.cardio_yes = heart_yes or angi_yes or ami_yes or cab_yes or pvd_yes or hbp_yes
    p.cardio_age = min_age(
        valid_age(b.get("IHD_ANGIAGE_NB_COM", "")),
        valid_age(b.get("IHD_AMIAGE_NB_COM", "")),
        valid_age(f1.get("CCC_HEARTAGE_NB_COF1", "")),
        valid_age(f1.get("CCC_PVDAGE_NB_COF1", "")),
        valid_age(f1.get("IHD_ANGIAGE_NB_COF1", "")),
        valid_age(f1.get("IHD_AMIAGE_NB_COF1", "")),
        valid_age(f2.get("CCC_HEARTAGE_NB_COF2", "")),
        valid_age(f2.get("IHD_ANGIAGE_NB_COF2", "")),
        valid_age(f2.get("IHD_AMIAGE_NB_COF2", "")),
    )

    # Cerebrovascular: stroke/CVA, TIA, and residual effects of stroke/TIA.
    cva_yes = any_true([
        yn_code(b.get("CCC_CVA_COM", "")),
        yn_code(f1.get("CCC_CVA_COF1", "")),
        yn_code(f2.get("CCC_CVA_COF2", "")),
    ])
    tia_yes = any_true([
        yn_code(b.get("CCC_TIA_COM", "")),
        yn_code(f1.get("CCC_TIA_COF1", "")),
        yn_code(f2.get("CCC_TIA_COF2", "")),
    ])
    cvafx_yes = any_true([
        yn_code(b.get("CCC_CVAFX_COM", "")),
        yn_code(f1.get("CCC_CVAFX_COF1", "")),
        yn_code(f2.get("CCC_CVAFX_COF2", "")),
    ])
    p.cerebro_yes = cva_yes or tia_yes or cvafx_yes
    p.cerebro_age = min_age(
        valid_age(b.get("STR_CVAAGE_NB_COM", "")),
        valid_age(b.get("STR_TIAAGE_NB_COM", "")),
        valid_age(f1.get("STR_CVAAGE_NB_COF1", "")),
        valid_age(f1.get("STR_TIAAGE_NB_COF1", "")),
        valid_age(f2.get("STR_CVAAGE_NB_COF2", "")),
        valid_age(f2.get("STR_TIAAGE_NB_COF2", "")),
    )

    p.neuro_available = False


@dataclass
class TraitDef:
    name: str
    label: str
    sensitivity: bool
    available: bool


def derive_trait_statuses(
    p: Participant, min_control_years: float, min_case_gap_years: float
) -> Dict[str, Optional[int]]:
    """Return binary trait values (1=case, 0=control, None=exclude) for all traits."""
    out: Dict[str, Optional[int]] = {}

    any_comp = p.ophth_yes or p.renal_yes or p.cardio_yes or p.cerebro_yes

    # Controls: T2D, complication-free, and T2D duration >= threshold.
    p.control_eligible = bool(
        p.t2d_flag
        and (not p.t1d_flag)
        and (not any_comp)
        and (p.t2d_duration is not None)
        and (p.t2d_duration >= min_control_years)
    )

    def main_case(comp_yes: bool, comp_age: Optional[float]) -> bool:
        if not p.t2d_flag or p.t1d_flag or not comp_yes:
            return False
        # If both ages are known, enforce ordering (complication at/after T2D diagnosis).
        if p.t2d_age is not None and comp_age is not None and comp_age < p.t2d_age:
            return False
        return True

    def sens_case(comp_yes: bool, comp_age: Optional[float]) -> bool:
        if not main_case(comp_yes, comp_age):
            return False
        if p.t2d_age is None or comp_age is None:
            return False
        return (comp_age - p.t2d_age) >= min_case_gap_years

    ophth_main = main_case(p.ophth_yes, p.ophth_age)
    ophth_sens = sens_case(p.ophth_yes, p.ophth_age)
    renal_main = main_case(p.renal_yes, p.renal_age)
    renal_sens = sens_case(p.renal_yes, p.renal_age)
    cardio_main = main_case(p.cardio_yes, p.cardio_age)
    cardio_sens = sens_case(p.cardio_yes, p.cardio_age)
    cerebro_main = main_case(p.cerebro_yes, p.cerebro_age)
    cerebro_sens = sens_case(p.cerebro_yes, p.cerebro_age)

    micro_main = ophth_main or renal_main  # Neurological unavailable in current CLSA mapping
    micro_sens = ophth_sens or renal_sens
    macro_main = cardio_main or cerebro_main
    macro_sens = cardio_sens or cerebro_sens
    micromacro_main = micro_main or macro_main
    micromacro_sens = micro_sens or macro_sens

    trait_cases = {
        "OPHTH_MAIN": ophth_main,
        "RENAL_MAIN": renal_main,
        "CARDIO_MAIN": cardio_main,
        "CEREBRO_MAIN": cerebro_main,
        "MICRO_MAIN": micro_main,
        "MACRO_MAIN": macro_main,
        "MICROMACRO_MAIN": micromacro_main,
        "OPHTH_SENS2Y": ophth_sens,
        "RENAL_SENS2Y": renal_sens,
        "CARDIO_SENS2Y": cardio_sens,
        "CEREBRO_SENS2Y": cerebro_sens,
        "MICRO_SENS2Y": micro_sens,
        "MACRO_SENS2Y": macro_sens,
        "MICROMACRO_SENS2Y": micromacro_sens,
    }

    for trait, is_case in trait_cases.items():
        if is_case:
            out[trait] = 1
        elif p.control_eligible:
            out[trait] = 0
        else:
            out[trait] = None

    # Explicitly unavailable neurological traits to keep manifest aligned with requested phenotypes.
    out["NEURO_MAIN"] = None
    out["NEURO_SENS2Y"] = None

    return out


TRAITS_ORDER = [
    "OPHTH_MAIN",
    "NEURO_MAIN",
    "RENAL_MAIN",
    "CARDIO_MAIN",
    "CEREBRO_MAIN",
    "MICRO_MAIN",
    "MACRO_MAIN",
    "MICROMACRO_MAIN",
    "OPHTH_SENS2Y",
    "NEURO_SENS2Y",
    "RENAL_SENS2Y",
    "CARDIO_SENS2Y",
    "CEREBRO_SENS2Y",
    "MICRO_SENS2Y",
    "MACRO_SENS2Y",
    "MICROMACRO_SENS2Y",
]


def sort_key_gwasid(x: str) -> Tuple[int, str]:
    try:
        return (0, f"{int(float(x)):012d}")
    except Exception:
        return (1, x)


def sanitize_header_token(s: str) -> str:
    """Make a whitespace-safe column token for covariate headers."""
    x = re.sub(r"[^A-Za-z0-9]+", "_", s.strip())
    x = re.sub(r"_+", "_", x).strip("_")
    return x or "UNK"


def write_tsv(path: str, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_keep_noheader(path: str, rows: List[Dict[str, object]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in rows:
            w.writerow([row["FID"], row["IID"]])


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    baseline = load_csv_selected(args.baseline_csv, BASELINE_COLS, key="entity_id")
    fup1 = load_csv_selected(args.fup1_csv, FUP1_COLS, key="entity_id")
    fup2 = load_csv_selected(args.fup2_csv, FUP2_COLS, key="entity_id")
    sqc = load_sqc(args.sqc)

    participants: List[Participant] = []

    for entity_id, b in baseline.items():
        gwasid = (b.get("ADM_GWAS3_COM", "") or "").strip()
        if not gwasid:
            continue
        if gwasid not in sqc:
            continue

        sex_mf = (b.get("SEX_ASK_COM", "") or "").strip()
        sex_regenie: Optional[int]
        if sex_mf == "M":
            sex_regenie = 1
        elif sex_mf == "F":
            sex_regenie = 0
        else:
            sex_regenie = None

        p = Participant(
            entity_id=entity_id,
            gwasid=gwasid,
            sex_mf=sex_mf,
            age_baseline=valid_age(b.get("AGE_NMBR_COM", "")),
            center=(b.get("GEOSTRATA_COM", "") or "").strip(),
            sqc=sqc[gwasid],
            baseline=b,
            fup1=fup1.get(entity_id, {}),
            fup2=fup2.get(entity_id, {}),
            sex_regenie=sex_regenie,
        )

        # Match OA pipeline QC conventions.
        if not all_covar_present(p):
            continue
        if parse_int_str(p.sqc.get("pca.cluster.id")) != args.pca_cluster_id:
            continue
        if parse_int_str(p.sqc.get("in.kinship")) != 0:
            continue
        if parse_int_str(p.sqc.get("in.hetmiss")) != 0:
            continue
        chrom_sex = parse_int_str(p.sqc.get("chromosomal.sex"))
        if chrom_sex is None:
            continue
        # baseline sex encoded M/F -> regenie 1/0; chromosomal.sex is 1=male,2=female.
        expected_chrom = 1 if p.sex_regenie == 1 else 2 if p.sex_regenie == 0 else None
        if expected_chrom is None or chrom_sex != expected_chrom:
            continue

        compute_t2d_flags(p, allow_without_type2=args.allow_t2d_without_type2)
        compute_complications(p)
        participants.append(p)

    # One-hot levels after QC filtering; drop first level as reference to avoid collinearity.
    center_levels = sorted({p.center for p in participants if p.center})
    center_ref = center_levels[0] if center_levels else None
    center_dummy_levels = [c for c in center_levels if c != center_ref]
    center_col_map = {c: f"center_{sanitize_header_token(c)}" for c in center_dummy_levels}
    batch_levels = sorted({parse_int_str(p.sqc.get("batch")) for p in participants if parse_int_str(p.sqc.get("batch")) is not None})
    batch_ref = batch_levels[0] if batch_levels else None
    batch_dummy_levels = [b for b in batch_levels if b != batch_ref]

    # Build a harmonized participant-level table (useful for auditing).
    harmonized_rows: List[Dict[str, object]] = []
    trait_status_by_gwasid: Dict[str, Dict[str, Optional[int]]] = {}

    for p in participants:
        traits = derive_trait_statuses(
            p,
            min_control_years=args.min_control_t2d_years,
            min_case_gap_years=args.min_case_gap_years,
        )
        trait_status_by_gwasid[p.gwasid] = traits
        row: Dict[str, object] = {
            "entity_id": p.entity_id,
            "gwasid": p.gwasid,
            "sex_mf": p.sex_mf,
            "sex": p.sex_regenie,
            "age_baseline": p.age_baseline,
            "center": p.center,
            "batch": parse_int_str(p.sqc.get("batch")),
            "pca.cluster.id": parse_int_str(p.sqc.get("pca.cluster.id")),
            "t1d_flag": int(p.t1d_flag),
            "t2d_flag": int(p.t2d_flag),
            "t2d_age": p.t2d_age,
            "last_obs_age": p.last_obs_age,
            "t2d_duration": p.t2d_duration,
            "control_eligible": int(p.control_eligible),
            "ophth_yes": int(p.ophth_yes),
            "ophth_age": p.ophth_age,
            "renal_yes": int(p.renal_yes),
            "renal_age": p.renal_age,
            "cardio_yes": int(p.cardio_yes),
            "cardio_age": p.cardio_age,
            "cerebro_yes": int(p.cerebro_yes),
            "cerebro_age": p.cerebro_age,
            "hba1c_baseline": valid_lab(p.baseline.get("BLD_HBA1c_COM", "")),
            "creat_baseline": valid_lab(p.baseline.get("BLD_CREAT_COM", "")),
            "alb_baseline": valid_lab(p.baseline.get("BLD_ALB_COM", "")),
        }
        for trait in TRAITS_ORDER:
            row[trait] = traits.get(trait)
        harmonized_rows.append(row)

    harmonized_rows.sort(key=lambda r: sort_key_gwasid(str(r["gwasid"])))
    harm_cols = list(harmonized_rows[0].keys()) if harmonized_rows else []
    with gzip.open(os.path.join(args.outdir, "participant_harmonized.tsv.gz"), "wt", newline="", encoding="utf-8") as gz:
        if harm_cols:
            w = csv.DictWriter(gz, fieldnames=harm_cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in harmonized_rows:
                w.writerow(row)

    # Prepare covariate base rows by participant.
    covar_base: Dict[str, Dict[str, object]] = {}
    for p in participants:
        row: Dict[str, object] = {
            "FID": p.gwasid,
            "IID": p.gwasid,
            "age": p.age_baseline,
            "sex": p.sex_regenie,
            "t2d_duration": p.t2d_duration if p.t2d_duration is not None else "NA",
        }
        for i, pc in enumerate(PC_COLS, start=1):
            row[pc] = p.sqc.get(f"ePC{i}")
        for center in center_dummy_levels:
            row[center_col_map[center]] = 1 if p.center == center else 0
        batch_val = parse_int_str(p.sqc.get("batch"))
        for b in batch_dummy_levels:
            row[f"genobatch{b}"] = 1 if batch_val == b else 0
        covar_base[p.gwasid] = row

    trait_available = {t: (not t.startswith("NEURO_")) for t in TRAITS_ORDER}

    manifest_rows: List[Dict[str, object]] = []
    regenie_targets: List[str] = []

    strata = {
        "both": lambda p: True,
        "male": lambda p: p.sex_regenie == 1,
        "female": lambda p: p.sex_regenie == 0,
    }

    # Build fast participant index by gwasid.
    participant_by_gwasid = {p.gwasid: p for p in participants}

    for sex_stratum, sex_fn in strata.items():
        subset_gwasids = sorted([p.gwasid for p in participants if sex_fn(p)], key=sort_key_gwasid)
        # Controls should be identical across traits within a stratum.
        stratum_controls = [g for g in subset_gwasids if participant_by_gwasid[g].control_eligible]

        for trait in TRAITS_ORDER:
            is_sens = trait.endswith("_SENS2Y")
            base_trait = trait.replace("_MAIN", "").replace("_SENS2Y", "")
            available = trait_available[trait]

            case_ids: List[str] = []
            for g in subset_gwasids:
                val = trait_status_by_gwasid[g].get(trait)
                if val == 1:
                    case_ids.append(g)

            n_controls = len(stratum_controls)
            n_cases = len(case_ids)
            n_total = n_controls + n_cases
            ready = int(available and n_cases >= args.min_cases_for_target and n_controls > 0)

            analysis_id = f"{args.ancestry_label}_{sex_stratum}_{trait}"
            manifest_rows.append(
                {
                    "analysis_id": analysis_id,
                    "ancestry": args.ancestry_label,
                    "sex_stratum": sex_stratum,
                    "trait": trait,
                    "trait_base": base_trait,
                    "sensitivity": int(is_sens),
                    "available": int(available),
                    "n_controls": n_controls,
                    "n_cases": n_cases,
                    "n_total": n_total,
                    "ready_for_gwas": ready,
                    "notes": (
                        "Neurological diabetic complication phenotype not identified in current CLSA variable mapping"
                        if not available
                        else (
                            "Micro excludes neurological (unavailable in CLSA)" if base_trait == "MICRO"
                            else "Retinopathy-only proxy; broader ophthalmic (cataract/glaucoma/macular edema) not captured" if base_trait == "OPHTH"
                            else "Self-report kidney disease/failure + dialysis proxy; not strict DKD/CKD by lab criteria" if base_trait == "RENAL"
                            else ""
                        )
                    ),
                }
            )

            if not ready:
                continue

            # Build phenotype/covar/keep rows.
            pheno_rows: List[Dict[str, object]] = []
            covar_rows: List[Dict[str, object]] = []
            keep_rows: List[Dict[str, object]] = []
            include_ids = set(stratum_controls) | set(case_ids)

            for g in sorted(include_ids, key=sort_key_gwasid):
                pheno_val = trait_status_by_gwasid[g][trait]
                if pheno_val not in (0, 1):
                    continue
                pheno_rows.append({"FID": g, "IID": g, trait: pheno_val})
                covar_rows.append(covar_base[g])
                keep_rows.append({"FID": g, "IID": g})

            pheno_path = os.path.join(args.outdir, f"{analysis_id}_pheno.tsv")
            covar_path = os.path.join(args.outdir, f"{analysis_id}_covar.tsv")
            keep_path = os.path.join(args.outdir, f"{analysis_id}_FIDIID_noheader.tsv")

            write_tsv(pheno_path, pheno_rows, ["FID", "IID", trait])
            base_covars = ["FID", "IID", "age", "t2d_duration"]
            if sex_stratum == "both":
                base_covars.append("sex")
            covar_fields = base_covars + PC_COLS + [center_col_map[c] for c in center_dummy_levels] + [f"genobatch{b}" for b in batch_dummy_levels]
            write_tsv(covar_path, covar_rows, covar_fields)
            write_keep_noheader(keep_path, keep_rows)
            regenie_targets.append(analysis_id)

    # Write manifest and targets.
    manifest_path = os.path.join(args.outdir, "phenotype_manifest.tsv")
    if manifest_rows:
        manifest_fields = [
            "analysis_id",
            "ancestry",
            "sex_stratum",
            "trait",
            "trait_base",
            "sensitivity",
            "available",
            "n_controls",
            "n_cases",
            "n_total",
            "ready_for_gwas",
            "notes",
        ]
        write_tsv(manifest_path, manifest_rows, manifest_fields)

    targets_path = os.path.join(args.outdir, "regenie_targets.txt")
    with open(targets_path, "w", encoding="utf-8") as fh:
        for target in sorted(set(regenie_targets)):
            fh.write(target + "\n")

    # Minimal console summary.
    print(f"Participants after QC/ancestry filter: {len(participants)}")
    if center_ref is not None:
        print(f"Center reference level (dropped from dummies): {center_ref}")
    if batch_ref is not None:
        print(f"Batch reference level (dropped from dummies): {batch_ref}")
    print(f"Wrote: {manifest_path}")
    print(f"Wrote: {targets_path} ({len(set(regenie_targets))} analyses)")


if __name__ == "__main__":
    main()
