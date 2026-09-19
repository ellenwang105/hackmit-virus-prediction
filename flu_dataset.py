"""Build the influenza antibody-antigen dataset from a SAbDab search export."""

import glob
import re

import pandas as pd

OUTPUT_PATH = "flu_dataset.csv"
MAX_RESOLUTION = 4.0
VALID_HA_SUBTYPES = set(range(1, 19))

FULL_SUBTYPE_RE = re.compile(r"\bH(\d{1,2})N(\d{1,2})\b", re.I)
BARE_HA_RE = re.compile(r"\bH(\d{1,2})\b", re.I)
# "parainfluenza" contains "influenza" and would otherwise be pulled in as a flu antigen
INFLUENZA_RE = r"(?<!para)influenza"
# engineered constructs such as cH5/1 carry sequence from two subtypes, so a single
# subtype label would leak the donor subtype across a leave-one-subtype-out split
CHIMERA_RE = re.compile(r"chimeric|chimera|\bcH\d{1,2}/\d{1,2}\b", re.I)
# pandemic years map unambiguously to an HA subtype; recorded as inferred so they can be excluded
HISTORICAL_YEARS = {"1918": "H1", "1957": "H2", "1968": "H3"}


def load_latest_export():
    candidates = sorted(glob.glob("sabdab_search_summary_*.csv"))
    if not candidates:
        raise FileNotFoundError("No sabdab_search_summary_*.csv file found in this directory.")
    path = candidates[-1]
    print(f"Loading {path}")
    return pd.read_csv(path, low_memory=False), path


def classify_target(antigen_name):
    name = str(antigen_name).lower()
    has_ha = "hemagglutinin" in name or "haemagglutinin" in name
    has_na = "neuraminidase" in name
    if has_ha and has_na:
        return "HA+NA"
    if has_ha:
        return "HA"
    if has_na:
        return "NA"
    return "other"


def resolve_subtype(row):
    """Return (raw_label, ha_subtype, source) — source records how confident the label is."""
    species = str(row["antigen_species"])
    compound = str(row["compound"])

    if CHIMERA_RE.search(compound):
        return "chimeric", "chimeric", "chimeric_construct"

    if re.search(r"influenza\s+b", species, re.I) or re.search(r"influenza\s+b", compound, re.I):
        return "B", "B", "species" if re.search(r"influenza\s+b", species, re.I) else "compound"

    for text, source in ((species, "species"), (compound, "compound")):
        match = FULL_SUBTYPE_RE.search(text)
        if match and int(match.group(1)) in VALID_HA_SUBTYPES:
            return match.group(0).upper(), f"H{int(match.group(1))}", source

    match = BARE_HA_RE.search(compound)
    if match and int(match.group(1)) in VALID_HA_SUBTYPES:
        return match.group(0).upper(), f"H{int(match.group(1))}", "compound_bare_H"

    for year, subtype in HISTORICAL_YEARS.items():
        if year in compound:
            return year, subtype, "historical_year"

    return "unknown", "unknown", "unresolved"


df, source_path = load_latest_export()
print(f"start: {len(df)} rows / {df['PDB'].nunique()} PDBs")

df = df[df["antigen_type"].str.contains("protein", case=False, na=False)]
df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
df = df[df["resolution"] <= MAX_RESOLUTION]
df = df.dropna(subset=["antigen_chain", "Hchain"]).copy()
print(f"after quality filters: {len(df)} rows / {df['PDB'].nunique()} PDBs")

df = df[df["antigen_species"].str.contains(INFLUENZA_RE, case=False, na=False, regex=True)].copy()
print(f"after influenza filter: {len(df)} rows / {df['PDB'].nunique()} PDBs")

df["date"] = pd.to_datetime(df["date"], errors="coerce")
df["year"] = df["date"].dt.year
df["antigen_chains"] = df["antigen_chain"].str.split(r"\s*\|\s*")
df["target"] = df["antigen_name"].map(classify_target)

resolved = df.apply(resolve_subtype, axis=1, result_type="expand")
df[["subtype_raw", "ha_subtype", "subtype_source"]] = resolved

columns = [
    "INSTANCE", "PDB", "SABDAB_ID", "Hchain", "Lchain", "antigen_chain", "antigen_chains",
    "target", "ha_subtype", "subtype_raw", "subtype_source",
    "antigen_name", "antigen_species", "compound",
    "date", "year", "resolution", "method",
    "VH", "CDR-H1", "CDR-H2", "CDR-H3", "VL", "CDR-L1", "CDR-L2", "CDR-L3",
]
df = df[columns].sort_values(["ha_subtype", "date"])

df.to_csv(OUTPUT_PATH, index=False)

print()
print(f"--- target ({len(df)} rows / {df['PDB'].nunique()} PDBs) ---")
print(df.groupby("target").agg(rows=("PDB", "size"), pdbs=("PDB", "nunique")).sort_values("rows", ascending=False))
print()
print("--- HA subtype ---")
print(df.groupby("ha_subtype").agg(rows=("PDB", "size"), pdbs=("PDB", "nunique")).sort_values("rows", ascending=False))
print()
print("--- how each subtype label was obtained ---")
print(df["subtype_source"].value_counts())
print()
ha = df[df["target"].isin(["HA", "HA+NA"])]
print(f"HA-only subset: {len(ha)} rows / {ha['PDB'].nunique()} PDBs")
print(f"  resolved subtype: {(ha['ha_subtype'] != 'unknown').sum()} rows")
print(f"  H5/H7 holdout candidates: {ha['ha_subtype'].isin(['H5', 'H7']).sum()} rows")
print()
print(f"wrote {len(df)} rows to {OUTPUT_PATH}")
