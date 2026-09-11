"""Per-molecule parsing, standardization, and cached descriptors.

One pass over the table builds MoleculeRecord objects reused by all checks.
"""
from __future__ import annotations

from dataclasses import dataclass, field


def _require_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
    except Exception as e:
        raise ImportError(
            "RDKit is required but not installed. Install with: "
            "conda install -c conda-forge rdkit"
        ) from e


@dataclass
class MoleculeRecord:
    idx: int
    row_id: str
    raw_smiles: str
    split: str | None
    label: object
    mol: object = None  # RDKit Mol or None
    valid: bool = False
    sanitize_error: str = ""
    canon_smi: str | None = None          # isomeric canonical
    connectivity_smi: str | None = None   # non-isomeric (stereo-agnostic)
    parent_smi: str | None = None         # desalted, isomeric
    parent_connectivity: str | None = None
    tautomer_smi: str | None = None
    scaffold: str | None = None
    fp: object = None
    mw: float | None = None
    logp: float | None = None
    n_atoms: int = 0
    formal_charge: int = 0
    n_fragments: int = 0
    elements: set = field(default_factory=set)
    has_isotope: bool = False
    has_radical: bool = False
    n_stereo_defined: int = 0
    n_stereo_unassigned: int = 0
    n_tautomers: int = 1
    max_ring_size: int = 0


def _classify_sanitize_error(msg: str) -> str:
    m = (msg or "").lower()
    if "valence" in m or "explicit valence" in m:
        return "valence"
    if "kekul" in m or "aromaticity" in m or "aromat" in m:
        return "aromaticity"
    if "charge" in m:
        return "charge"
    return "other"


def build_records(df) -> list[MoleculeRecord]:
    _require_rdkit()
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    # shared helpers (constructed once)
    salt_remover = SaltRemover.SaltRemover()
    taut_enumerator = rdMolStandardize.TautomerEnumerator()
    # cap enumeration cost
    try:
        taut_enumerator.SetMaxTautomers(32)
    except Exception:
        pass

    records: list[MoleculeRecord] = []
    for _, r in df.iterrows():
        idx = int(r["_row"])
        rec = MoleculeRecord(
            idx=idx, row_id=str(r["_id"]), raw_smiles=str(r["_smiles"]),
            split=str(r["_split"]) if r["_split"] is not None else None,
            label=r["_label"],
        )
        smi = rec.raw_smiles
        if smi in ("", "nan", "None", "none"):
            rec.sanitize_error = "empty SMILES"
            records.append(rec)
            continue
        # parse without sanitize to classify errors stepwise
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            # retry with sanitize to get here; definitely invalid
            rec.sanitize_error = "RDKit could not parse SMILES"
            records.append(rec)
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            # try to sanitize with catch to assign kekulize etc; record error kind
            rec.sanitize_error = f"{type(e).__name__}: {e} [{_classify_sanitize_error(str(e))}]"
            # keep molecule attempt via default parse for partial info? mark invalid
            records.append(rec)
            continue
        rec.mol = mol
        rec.valid = True
        try:
            rec.canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            no_iso = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
            # strip stereo atoms remnants: also remove / \ @ by construction above
            rec.connectivity_smi = no_iso
            rec.n_atoms = mol.GetNumAtoms()
            rec.formal_charge = Chem.GetFormalCharge(mol)
            rec.elements = {a.GetSymbol() for a in mol.GetAtoms()}
            rec.has_isotope = any(a.GetIsotope() != 0 for a in mol.GetAtoms())
            rec.has_radical = any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms())
            frags = Chem.GetMolFrags(mol, asMols=False)
            rec.n_fragments = len(frags)
            # parent (desalted): largest fragment by atom count
            try:
                parent = salt_remover.StripMol(mol, dontRemoveEverything=True)
                rec.parent_smi = Chem.MolToSmiles(parent, isomericSmiles=True, canonical=True)
                rec.parent_connectivity = Chem.MolToSmiles(parent, isomericSmiles=False, canonical=True)
            except Exception:
                rec.parent_smi = rec.canon_smi
                rec.parent_connectivity = rec.connectivity_smi
            # scaffold
            try:
                scaf = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(rec.canon_smi or smi)
                rec.scaffold = scaf if scaf else "(acyclic)"
            except Exception:
                rec.scaffold = None
            # tautomer canonical
            try:
                taut = taut_enumerator.Canonicalize(mol)
                rec.tautomer_smi = Chem.MolToSmiles(taut, isomericSmiles=True, canonical=True)
            except Exception:
                rec.tautomer_smi = rec.canon_smi
            try:
                rec.n_tautomers = len(taut_enumerator.Enumerate(mol))
            except Exception:
                rec.n_tautomers = 1
            # stereo: assigned vs unassigned
            try:
                Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
                centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                rec.n_stereo_defined = sum(1 for _, v in centers if v in ("R", "S"))
                rec.n_stereo_unassigned = sum(1 for _, v in centers if v == "?")
                # also catch E/Z unassigned via bond stereo? thin: check unspecified double bonds
                for b in mol.GetBonds():
                    if b.GetBondType() == Chem.BondType.DOUBLE:
                        st = b.GetStereo()
                        if st == Chem.BondStereo.STEREONONE:
                            # only count if both sides have >1 neighbor (potentially stereo)
                            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
                            if a1.GetDegree() > 1 and a2.GetDegree() > 1:
                                # check substituents differ — thin heuristic, count once per mol max
                                pass
            except Exception:
                pass
            # descriptors + fp + rings
            try:
                rec.mw = Descriptors.MolWt(mol)
                rec.logp = Descriptors.MolLogP(mol)
            except Exception:
                pass
            try:
                from rdkit.Chem import AllChem
                rec.fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            except Exception:
                rec.fp = None
            try:
                ri = mol.GetRingInfo()
                rec.max_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
            except Exception:
                pass
        except Exception:
            # never let derived-property failure invalidate the molecule
            pass
        records.append(rec)
    return records
