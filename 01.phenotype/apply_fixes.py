#!/usr/bin/env python3
"""Apply reviewer fixes to 01.curate_t2d_complications.py"""
import os

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "01.curate_t2d_complications.py")

with open(SCRIPT, "r") as f:
    c = f.read()

applied = []

# 1. Add CCC_HBP_COF1 to FUP1_COLS
old = '    "CCC_PVDAGE_NB_COF1",\n    "CCC_KIDN_COF1",'
new = '    "CCC_PVDAGE_NB_COF1",\n    "CCC_HBP_COF1",\n    "CCC_KIDN_COF1",'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("1: Added CCC_HBP_COF1 to FUP1_COLS")

# 2a. Update cardio comment
old = '    # Cardiovascular: heart disease, angina, MI, CAB/revascularization, PAD/PVD.'
new = '    # Cardiovascular: heart disease, angina, MI, CAB/revascularization, PAD/PVD, hypertension.'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("2a: Updated cardio comment")

# 2b. Add hbp_yes + update cardio_yes
old = '    pvd_yes = any_true([\n        yn_code(b.get("CCC_PVD_COM", "")),\n        yn_code(f1.get("CCC_PVD_COF1", "")),\n        # No PVD variable appears in FUP2 CoPv3 dictionary.\n    ])\n    p.cardio_yes = heart_yes or angi_yes or ami_yes or cab_yes or pvd_yes'
new = '    pvd_yes = any_true([\n        yn_code(b.get("CCC_PVD_COM", "")),\n        yn_code(f1.get("CCC_PVD_COF1", "")),\n    ])\n    hbp_yes = any_true([\n        yn_code(b.get("CCC_HBP_COM", "")),\n        yn_code(f1.get("CCC_HBP_COF1", "")),\n    ])\n    p.cardio_yes = heart_yes or angi_yes or ami_yes or cab_yes or pvd_yes or hbp_yes'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("2b: Added hbp_yes, updated p.cardio_yes")

# 3. Drop reference levels for center and batch
old = '    # One-hot levels after QC filtering.\n    center_levels = sorted({p.center for p in participants if p.center})\n    center_col_map = {c: f"center_{sanitize_header_token(c)}" for c in center_levels}\n    batch_levels = sorted({parse_int_str(p.sqc.get("batch")) for p in participants if parse_int_str(p.sqc.get("batch")) is not None})'
new = '    # One-hot levels after QC filtering; drop first level as reference to avoid collinearity.\n    center_levels = sorted({p.center for p in participants if p.center})\n    center_ref = center_levels[0] if center_levels else None\n    center_dummy_levels = [c for c in center_levels if c != center_ref]\n    center_col_map = {c: f"center_{sanitize_header_token(c)}" for c in center_dummy_levels}\n    batch_levels = sorted({parse_int_str(p.sqc.get("batch")) for p in participants if parse_int_str(p.sqc.get("batch")) is not None})\n    batch_ref = batch_levels[0] if batch_levels else None\n    batch_dummy_levels = [b for b in batch_levels if b != batch_ref]'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("3: Drop reference level for center/batch dummies")

# 4. Add t2d_duration to covar_base + use dummy_levels
old = '            "sex": p.sex_regenie,\n        }\n        for i, pc in enumerate(PC_COLS, start=1):\n            row[pc] = p.sqc.get(f"ePC{i}")\n        for center in center_levels:\n            row[center_col_map[center]] = 1 if p.center == center else 0\n        batch_val = parse_int_str(p.sqc.get("batch"))\n        for b in batch_levels:\n            row[f"genobatch{b}"] = 1 if batch_val == b else 0'
new = '            "sex": p.sex_regenie,\n            "t2d_duration": p.t2d_duration if p.t2d_duration is not None else "NA",\n        }\n        for i, pc in enumerate(PC_COLS, start=1):\n            row[pc] = p.sqc.get(f"ePC{i}")\n        for center in center_dummy_levels:\n            row[center_col_map[center]] = 1 if p.center == center else 0\n        batch_val = parse_int_str(p.sqc.get("batch"))\n        for b in batch_dummy_levels:\n            row[f"genobatch{b}"] = 1 if batch_val == b else 0'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("4: Added t2d_duration to covar_base, use dummy_levels")

# 5. Update covar_fields
old = '            covar_fields = ["FID", "IID", "age", "sex"] + PC_COLS + [center_col_map[c] for c in center_levels] + [f"genobatch{b}" for b in batch_levels]'
new = '            base_covars = ["FID", "IID", "age", "t2d_duration"]\n            if sex_stratum == "both":\n                base_covars.append("sex")\n            covar_fields = base_covars + PC_COLS + [center_col_map[c] for c in center_dummy_levels] + [f"genobatch{b}" for b in batch_dummy_levels]'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("5: Updated covar_fields (t2d_duration, conditional sex, dummy_levels)")

# 6. Update manifest notes
old = '                    "notes": (\n                        "Neurological diabetic complication phenotype not identified in current CLSA variable mapping"\n                        if not available\n                        else ("Micro excludes neurological (unavailable)" if base_trait == "MICRO" else "")\n                    ),'
new = '                    "notes": (\n                        "Neurological diabetic complication phenotype not identified in current CLSA variable mapping"\n                        if not available\n                        else (\n                            "Micro excludes neurological (unavailable in CLSA)" if base_trait == "MICRO"\n                            else "Retinopathy-only proxy; broader ophthalmic (cataract/glaucoma/macular edema) not captured" if base_trait == "OPHTH"\n                            else "Self-report kidney disease/failure + dialysis proxy; not strict DKD/CKD by lab criteria" if base_trait == "RENAL"\n                            else ""\n                        )\n                    ),'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("6: Updated manifest notes for OPHTH, RENAL, MICRO")

# 7. Add reference level info to console summary
old = '    # Minimal console summary.\n    print(f"Participants after QC/ancestry filter: {len(participants)}")\n    print(f"Wrote: {manifest_path}")\n    print(f"Wrote: {targets_path} ({len(set(regenie_targets))} analyses)")'
new = '    # Minimal console summary.\n    print(f"Participants after QC/ancestry filter: {len(participants)}")\n    if center_ref is not None:\n        print(f"Center reference level (dropped from dummies): {center_ref}")\n    if batch_ref is not None:\n        print(f"Batch reference level (dropped from dummies): {batch_ref}")\n    print(f"Wrote: {manifest_path}")\n    print(f"Wrote: {targets_path} ({len(set(regenie_targets))} analyses)")'
if old in c:
    c = c.replace(old, new, 1)
    applied.append("7: Added reference level info to console summary")

with open(SCRIPT, "w") as f:
    f.write(c)

print(f"Applied {len(applied)} of 7 change groups:")
for a in applied:
    print(f"  OK: {a}")
if len(applied) < 7:
    print("WARNING: Some changes could not be applied (pattern not found)")
