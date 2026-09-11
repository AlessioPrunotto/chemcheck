"""Input loading: CSV/TSV/SDF/Parquet/Excel/JSONL + split handling."""
from __future__ import annotations

import os

import pandas as pd

SMILES_CANDIDATES = [
    "smiles", "canonical_smiles", "SMILES", "Smiles", "SMILES_CANONICAL",
    "mol_smiles", "structure", "Structure", "molecule", "compound_smiles",
]

SPLIT_CANDIDATES = ["split", "Split", "SPLIT", "set", "fold", "subset"]
ID_CANDIDATES = ["id", "ID", "Id", "name", "Name", "mol_id", "compound_id", "cid", "_sdf_name", "title"]
LABEL_CANDIDATES = ["activity", "label", "target", "value", "pIC50", "pIC50_value", "y"]


def _rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except Exception:
        return False


def _guess_smiles_column(df: pd.DataFrame) -> str | None:
    for c in SMILES_CANDIDATES:
        if c in df.columns:
            return c
    # heuristic: string column with best SMILES parse rate on a sample
    best, best_rate = None, 0.0
    try:
        from rdkit import Chem
    except Exception:
        # without RDKit, fall back to first string column
        for c in df.columns:
            if df[c].dtype == object:
                return c
        return None
    for c in df.columns:
        if df[c].dtype != object:
            continue
        sample = df[c].dropna().astype(str).head(20).tolist()
        if not sample:
            continue
        ok = sum(1 for s in sample if Chem.MolFromSmiles(s) is not None)
        rate = ok / max(1, len(sample))
        if rate > best_rate:
            best, best_rate = c, rate
    return best if best_rate >= 0.3 else None


def _read_single(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".sd"):
        return _read_sdf(path)
    if ext in (".csv", ".txt", ".tsv"):
        if ext == ".tsv":
            return pd.read_csv(path, sep="\t")
        if ext == ".txt":
            with open(path) as fh:
                head = fh.read(4096)
            if "\t" in head and "," not in head.splitlines()[0]:
                return pd.read_csv(path, sep="\t")
        return pd.read_csv(path)
    if ext in (".parquet", ".pq"):
        return pd.read_parquet(path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if ext in (".jsonl", ".ndjson"):
        return pd.read_json(path, lines=True)
    if ext == ".json":
        try:
            return pd.read_json(path, lines=True)
        except ValueError:
            return pd.read_json(path)
    raise ValueError(f"Unsupported input format: {ext} ({path})")


def _read_sdf(path: str) -> pd.DataFrame:
    try:
        from rdkit import Chem
    except Exception as e:
        raise ImportError(
            "RDKit is required but not installed. Install with: "
            "conda install -c conda-forge rdkit"
        ) from e
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    rows = []
    for i, mol in enumerate(suppl):
        if mol is None:
            rows.append({"_sdf_name": f"mol_{i}", "smiles": ""})
            continue
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            smi = ""
        row = {"_sdf_name": mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}",
               "smiles": smi}
        try:
            for k in mol.GetPropNames():
                row[k] = mol.GetProp(k)
        except Exception:
            pass
        rows.append(row)
    return pd.DataFrame(rows)


def _normalize(df: pd.DataFrame, smiles_col: str | None, id_col: str | None,
               label_col: str | None, split_col: str | None,
               split_value: str | None) -> pd.DataFrame:
    out = df.copy()
    if smiles_col is None:
        smiles_col = _guess_smiles_column(out)
    if smiles_col is None or smiles_col not in out.columns:
        raise ValueError(
            f"Could not detect SMILES column (tried {SMILES_CANDIDATES}). "
            "Pass --smiles-col explicitly."
        )
    if id_col is not None and id_col not in out.columns:
        raise ValueError(f"--id-col '{id_col}' not found in columns {list(out.columns)[:20]}")
    if label_col is not None and label_col not in out.columns:
        raise ValueError(f"--label-col '{label_col}' not found in columns {list(out.columns)[:20]}")
    if split_col is not None and split_col not in out.columns:
        raise ValueError(f"--split-col '{split_col}' not found in columns {list(out.columns)[:20]}")

    if id_col is None:
        for c in ID_CANDIDATES:
            if c in out.columns:
                id_col = c
                break
    out["_row"] = range(len(out))
    out["_id"] = out[id_col].astype(str) if id_col else [f"row_{i}" for i in range(len(out))]
    out["_smiles"] = out[smiles_col].astype(str).str.strip()
    out["_split"] = out[split_col].astype(str) if split_col else (split_value or None)
    out["_label"] = out[label_col] if label_col else None
    out.attrs["smiles_col"] = smiles_col
    out.attrs["id_col"] = id_col
    out.attrs["label_col"] = label_col
    out.attrs["split_col"] = split_col
    return out


def load_table(paths: list[str], smiles_col: str | None = None,
               id_col: str | None = None, label_col: str | None = None,
               split_col: str | None = None) -> pd.DataFrame:
    """Load one file, or two files as train/test splits."""
    if len(paths) == 1:
        df = _read_single(paths[0])
        # autodetect split col if not given
        if split_col is None:
            for c in SPLIT_CANDIDATES:
                if c in df.columns:
                    split_col = c
                    break
        norm = _normalize(df, smiles_col, id_col, label_col, split_col, None)
        norm.attrs["source"] = paths[0]
        return norm
    elif len(paths) == 2:
        a = _normalize(_read_single(paths[0]), smiles_col, id_col, label_col, None, "train")
        b = _normalize(_read_single(paths[1]), smiles_col, id_col, label_col, None, "test")
        # align label col presence: if one side lacks labels it's fine
        merged = pd.concat([a, b], ignore_index=True)
        merged["_row"] = range(len(merged))
        merged.attrs["source"] = f"{paths[0]} + {paths[1]}"
        merged.attrs["split_col"] = split_col or "(two-file train/test)"
        return merged
    raise ValueError("Pass 1 dataset file, or 2 files as train/test.")
