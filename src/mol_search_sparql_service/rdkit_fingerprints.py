import csv
import os
import tempfile
from typing import Any, Callable
from dataclasses import dataclass, replace
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import (
    rdFingerprintGenerator,
    RDKFingerprint,
    MACCSkeys,
    PatternFingerprint,
)
from rdkit.Chem.AtomPairs import Pairs, Torsions


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class FingerprintExplainability:
    """Explainability information for a fingerprint type, used for documentation and UI hints."""

    level: str
    mechanism: str
    limitations: str
    typical_explanations: list[str]


@dataclass
class FingerprintConfig:
    """Configuration for a single fingerprint type."""

    short_name: str
    python_method: Any  # Callable — RDKit generator / fingerprint function
    default_options: dict[str, Any]
    stereo_options: dict[str, Any]
    description: str
    explainability: FingerprintExplainability


@dataclass
class CompoundEntry:
    """A single compound with its precomputed fingerprint."""

    id: str
    smiles: str
    db_name: str = "unknown"
    fp: Any = None  # RDKit fingerprint object (ExplicitBitVect / IntSparseIntVect)


@dataclass
class Dataset:
    """All precomputed data for one fingerprint type."""

    data: list[CompoundEntry]
    fps: list[Any]  # parallel list of raw RDKit FP objects
    db_indices: dict[str, list[int]]  # db_name → list of indices into data/fps


@dataclass
class SimilarityResult:
    """One hit from a similarity search."""

    compound: CompoundEntry
    similarity: float


@dataclass
class SubstructureResult:
    """One hit from a substructure search (compound data + match count)."""

    id: str
    smiles: str
    db_name: str
    fp: Any = None
    match_count: int = 0


# ---------------------------------------------------------------------------
# Fingerprint registry
# ---------------------------------------------------------------------------

FINGERPRINTS: dict[str, FingerprintConfig] = {
    "morgan_ecfp": FingerprintConfig(
        short_name="ECFP",
        python_method=rdFingerprintGenerator.GetMorganGenerator,
        default_options={
            "radius": 2,
            "fpSize": 2048,
            "includeChirality": False,
            "useBondTypes": True,
            "countSimulation": False,
        },
        stereo_options={"includeChirality": True},
        description=(
            "Extended Connectivity Fingerprint (ECFP). "
            "Encodes atom-centered circular environments up to a given radius. "
            "Widely used for similarity search, clustering, and QSAR."
        ),
        explainability=FingerprintExplainability(
            level="high",
            mechanism=(
                "Each bit corresponds to one or more atom-centered environments "
                "(atom index + radius). Bit-to-substructure mapping is available "
                "via additionalOutput / bitInfo."
            ),
            limitations=(
                "Bits are hashed; collisions are possible. "
                "One bit may correspond to multiple distinct substructures."
            ),
            typical_explanations=[
                "Highlighted atom environments",
                "Similarity maps",
                "Per-atom importance aggregation",
            ],
        ),
    ),
    "morgan_fcfp": FingerprintConfig(
        short_name="FCFP",
        python_method=rdFingerprintGenerator.GetMorganGenerator,
        default_options={
            "radius": 2,
            "fpSize": 2048,
            "includeChirality": False,
            "useBondTypes": True,
            "countSimulation": False,
        },
        stereo_options={"includeChirality": True},
        description=(
            "Functional-Class Fingerprint (FCFP). "
            "Morgan fingerprint using pharmacophoric atom features instead of "
            "exact atom types."
        ),
        explainability=FingerprintExplainability(
            level="high",
            mechanism=(
                "Same as ECFP, but environments are defined over functional "
                "roles (HBD, HBA, aromatic, charged, etc.)."
            ),
            limitations="Chemical specificity is reduced compared to ECFP.",
            typical_explanations=[
                "Functional similarity",
                "Scaffold hopping rationales",
            ],
        ),
    ),
    "rdk_topological": FingerprintConfig(
        short_name="RDK",
        python_method=RDKFingerprint,
        default_options={
            "minPath": 1,
            "maxPath": 7,
            "fpSize": 2048,
            "useHs": True,
            "branchedPaths": True,
        },
        stereo_options={},
        description=(
            "RDKit topological (path-based) fingerprint. "
            "Encodes linear bond paths similar to Daylight fingerprints."
        ),
        explainability=FingerprintExplainability(
            level="high",
            mechanism=(
                "Each bit corresponds to one or more explicit bond paths. "
                "Exact atom and bond indices can be recovered via bitInfo."
            ),
            limitations=(
                "Sensitive to small structural changes; "
                "less robust for scaffold hopping."
            ),
            typical_explanations=[
                "Exact substructure paths",
                "Bond-path highlighting",
            ],
        ),
    ),
    "atom_pair": FingerprintConfig(
        short_name="AP",
        python_method=Pairs.GetAtomPairFingerprint,  # type: ignore[attr-defined]
        default_options={},
        stereo_options={},
        description=(
            "Atom Pair fingerprint. Encodes pairs of atoms along with their "
            "topological distance."
        ),
        explainability=FingerprintExplainability(
            level="medium",
            mechanism=(
                "Each feature represents a pair of atoms at a given distance. "
                "Explanations identify which atom pairs contributed."
            ),
            limitations=(
                "No connected subgraph; explanations are relational rather than "
                "structural."
            ),
            typical_explanations=[
                "Activity cliff analysis",
                "Long-range interaction reasoning",
            ],
        ),
    ),
    "topological_torsion": FingerprintConfig(
        short_name="TT",
        python_method=Torsions.GetTopologicalTorsionFingerprint,  # type: ignore[attr-defined]
        default_options={},
        stereo_options={},
        description="Topological Torsion fingerprint. Encodes sequences of four bonded atoms.",
        explainability=FingerprintExplainability(
            level="medium",
            mechanism=(
                "Each feature corresponds to a specific 4-atom sequence (A–B–C–D)."
            ),
            limitations="Local view only; torsions are hashed in bit-vector form.",
            typical_explanations=[
                "Linker characterization",
                "Conformation-sensitive similarity",
            ],
        ),
    ),
    "maccs": FingerprintConfig(
        short_name="MACCS",
        python_method=MACCSkeys.GenMACCSKeys,  # type: ignore[attr-defined]
        default_options={},
        stereo_options={},
        description=(
            "MACCS structural keys (166 bits). "
            "Each bit corresponds to a predefined chemical pattern."
        ),
        explainability=FingerprintExplainability(
            level="very high",
            mechanism=(
                "Each bit has a fixed semantic meaning defined in the MACCS "
                "specification."
            ),
            limitations="Low resolution; many subtle SAR effects are not captured.",
            typical_explanations=[
                "Human-readable feature presence",
                "Medicinal chemistry reports",
            ],
        ),
    ),
    "pattern": FingerprintConfig(
        short_name="Pattern",
        python_method=PatternFingerprint,
        default_options={"fpSize": 2048},
        stereo_options={},
        description=("RDKit Pattern fingerprint. Designed for substructure screening."),
        explainability=FingerprintExplainability(
            level="low",
            mechanism=(
                "Bits correspond to various small substructures/paths. "
                "Mainly used for pre-filtering substructure matches."
            ),
            limitations="High collision rate; screening only.",
            typical_explanations=["Substructure screening"],
        ),
    ),
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def get_fingerprint(
    mol: Chem.Mol, name: str = "morgan_ecfp", stereo: bool = False
) -> Any:
    """Generates a fingerprint for a molecule using the specified configuration name.

    Returns:
        An RDKit fingerprint object (`ExplicitBitVect` or `IntSparseIntVect`).
    """
    if name not in FINGERPRINTS:
        raise ValueError(f"Unknown fingerprint type: {name}")

    cfg = FINGERPRINTS[name]
    opts = cfg.default_options.copy()
    if stereo:
        opts.update(cfg.stereo_options)

    # Handle FCFP (feature invariants)
    if name == "morgan_fcfp":
        opts["atomInvariantsGenerator"] = (
            rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        )

    func: Callable[..., Any] = cfg.python_method

    if name in ["morgan_ecfp", "morgan_fcfp"]:
        generator = func(**opts)
        return generator.GetFingerprint(mol)
    else:
        # For others (including pattern), options are passed directly to the function along with mol
        return func(mol, **opts)


def safe_mol_from_smiles(smiles: str, cid: str = "unknown") -> Chem.Mol | None:
    """
    Safely parses a SMILES string into an RDKit Mol object, avoiding strict
    sanitization bugs on organo-metalic compounds by carefully applying SanitizeFlags.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None

    problems = Chem.DetectChemistryProblems(mol)
    if problems:
        for p in problems:
            print(f"{cid}\terror\t{p.GetType()}: {p.Message()}")

    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            ^ Chem.SanitizeFlags.SANITIZE_CLEANUP,
        )
        return mol
    except Exception as e:
        print(f"{cid}\terror\t{str(e)}")
        return None


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------


class MolSearchEngine:
    def __init__(self) -> None:
        # keys: fingerprint type name (e.g. 'morgan_ecfp', 'pattern')
        # values: compiled Dataset for that fingerprint type
        self.datasets: dict[str, Dataset] = {}

    def add_data(self, data: list[CompoundEntry], fp_type: str) -> None:
        """Populate the engine with a list of CompoundEntry objects.

        Each entry must have a precomputed 'fp' and optionally a 'db_name'.
        """
        # Optimization: Pre-calculate indices for each db_name
        db_indices: dict[str, list[int]] = {}
        fps: list[Any] = []

        for idx, entry in enumerate(data):
            # Extract FP for bulk operations
            fps.append(entry.fp)

            # Index by db_name
            if entry.db_name not in db_indices:
                db_indices[entry.db_name] = []
            db_indices[entry.db_name].append(idx)

        self.datasets[fp_type] = Dataset(data=data, fps=fps, db_indices=db_indices)

    def _get_indices(
        self, fp_type: str, db_names: list[str] | None = None
    ) -> range | list[int]:
        """Helper to get valid indices based on db_names filter.

        Returns:
            A range (all) or a filtered list of integer indices.
        """
        dataset = self.datasets.get(fp_type)
        if not dataset:
            raise ValueError(f"Dataset {fp_type} not loaded.")

        if not db_names:
            # Return range covering all indices
            return range(len(dataset.data))

        indices: list[int] = []
        for db in db_names:
            if db in dataset.db_indices:
                indices.extend(dataset.db_indices[db])
        return indices

    def search_similarity(
        self,
        query_smiles: str,
        limit: int = 5,
        db_names: list[str] | None = None,
        fp_type: str = "morgan_ecfp",
        use_chirality: bool = False,
        min_score: float = 0.0,
    ) -> list[SimilarityResult]:
        """Executes a similarity search."""
        if fp_type not in self.datasets:
            raise ValueError(f"Error: Dataset {fp_type} not loaded.")

        dataset = self.datasets[fp_type]

        query_mol = safe_mol_from_smiles(query_smiles, cid="query")
        if not query_mol:
            print(f"Error: Invalid query SMILES: {query_smiles}")
            return []

        query_fp = get_fingerprint(query_mol, fp_type, stereo=use_chirality)

        # Get indices to search
        indices = self._get_indices(fp_type, db_names)

        # Retrieve FPs for the target indices
        if not db_names:
            target_fps: list[Any] = dataset.fps
            target_data: list[CompoundEntry] = dataset.data
        else:
            # Construct subset list
            target_fps = [dataset.fps[i] for i in indices]
            target_data = [dataset.data[i] for i in indices]

        if not target_fps:
            return []

        sims: list[float] = DataStructs.BulkTanimotoSimilarity(query_fp, target_fps)

        results: list[SimilarityResult] = []
        for i, sim in enumerate(sims):
            if sim >= min_score:
                results.append(
                    SimilarityResult(compound=target_data[i], similarity=sim)
                )

        results.sort(key=lambda x: x.similarity, reverse=True)
        return results[:limit]

    def search_substructure(
        self,
        query_smiles: str,
        limit: int = 5,
        use_chirality: bool = False,
        db_names: list[str] | None = None,
        fp_type: str = "pattern",
        min_match_count: int = 1,
    ) -> list[SubstructureResult]:
        """Executes a substructure search (Screening + Verification)."""
        if fp_type not in self.datasets:
            raise ValueError(f"Error: Dataset {fp_type} not loaded.")

        dataset = self.datasets[fp_type]

        query_mol = safe_mol_from_smiles(query_smiles, cid="query")
        if not query_mol:
            print(f"Error: Invalid query SMILES: {query_smiles}")
            return []

        query_fp = get_fingerprint(query_mol, fp_type)

        # Get candidate indices
        indices = self._get_indices(fp_type, db_names)

        candidates_data: list[CompoundEntry] = []

        # Screening
        # Optimization: We avoid creating intermediate lists of ALL compounds
        data = dataset.data
        fps = dataset.fps

        for i in indices:
            if DataStructs.AllProbeBitsMatch(query_fp, fps[i]):
                candidates_data.append(data[i])

        # Verification
        results: list[SubstructureResult] = []
        for entry in candidates_data:
            try:
                # Optimized verification: Parsing SMILES is the slow part.
                target_mol = safe_mol_from_smiles(entry.smiles, cid=entry.id)
                if target_mol is None:
                    continue
                matches = target_mol.GetSubstructMatches(
                    query_mol, useChirality=use_chirality
                )

                if matches and len(matches) >= min_match_count:
                    results.append(
                        SubstructureResult(
                            id=entry.id,
                            smiles=entry.smiles,
                            db_name=entry.db_name,
                            fp=entry.fp,
                            match_count=len(matches),
                        )
                    )

                    if limit and len(results) >= limit:
                        break
            except Exception:
                continue
        return results

    def load_from_sparql(self, endpoint: str, query: str) -> None:
        """
        Fetch compound data from a SPARQL endpoint, store it in a temp TSV file,
        and load it into the engine. The temp file path is stored in the
        COMPOUNDS_FILE environment variable so that additional Uvicorn workers can
        reload the same data on startup without re-querying the endpoint.
        """
        import requests

        print(f"Fetching data from {endpoint}...")
        try:
            resp = requests.post(
                endpoint,
                data={"query": query},
                headers={"Accept": "text/tab-separated-values"},
            )
            resp.raise_for_status()
        except Exception as e:
            raise RuntimeError(f"Failed to fetch data from SPARQL endpoint: {e}") from e

        fd, temp_path = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(fd, "wb") as f:
            f.write(resp.content)

        # Signal workers to delete the temp file on shutdown
        os.environ["DELETE_COMPOUNDS_FILE"] = "1"
        self.load_file(temp_path)

    def load_file(self, compounds_file: str) -> None:
        # Persist the path so that additional Uvicorn workers can pick it up
        os.environ["COMPOUNDS_FILE"] = compounds_file
        print(f"Reading compounds from {compounds_file}...")
        compounds: list[CompoundEntry] = []
        with open(compounds_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # New format: ?db, ?chem, ?smiles
                try:
                    # Extract ID from ?chem (<URI>)
                    cid = row.get("?chem", "").strip("<>")
                    if not cid:
                        continue

                    # Extract SMILES from ?smiles ("SMILES")
                    smiles = row.get("?smiles", "").strip('"')

                    # Extract DB from ?db (<URI>)
                    db = row.get("?db", "").strip("<>") if "?db" in row else "unknown"

                    compounds.append(CompoundEntry(id=cid, smiles=smiles, db_name=db))
                except Exception:
                    continue

        # Compile In-Memory
        print(f"Compiling fingerprints dynamically ({len(compounds)} compounds)...")

        for fp_name in FINGERPRINTS.keys():
            print(f"  - Compiling {fp_name}...")
            try:
                data = compile_fingerprints_in_memory(compounds, fp_name)
                self.add_data(data, fp_name)
            except Exception as e:
                print(f"    Error compiling {fp_name}: {e}")

        print("Compilation complete.")


def compile_fingerprints_in_memory(
    compounds: list[CompoundEntry], fp_type: str
) -> list[CompoundEntry]:
    """Compiles fingerprints for a list of `CompoundEntry` objects in-memory without writing to disk."""
    data: list[CompoundEntry] = []
    for entry in compounds:
        mol = safe_mol_from_smiles(entry.smiles, cid=entry.id)
        if mol:
            fp = get_fingerprint(mol, fp_type)
            # Create a new entry with the fingerprint attached
            data.append(replace(entry, fp=fp))
    return data


# Initialize Engine globally for easy sharing
engine = MolSearchEngine()
