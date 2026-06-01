import contextlib
import csv
import hashlib
import os
import pickle
import tempfile
from typing import Any, Callable
from dataclasses import dataclass, field, replace
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import (
    rdFingerprintGenerator,
    RDKFingerprint,
    MACCSkeys,
    PatternFingerprint,
)

# NOTE: we need to silence RDKit warnings that magically poped up from nowhere
# even if it was working before, and no libs version have been changed
@contextlib.contextmanager
def _silence_stderr():
    """Redirect fd 2 to /dev/null to suppress C++-level stderr output."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = os.dup(2)
    os.dup2(devnull, 2)
    try:
        yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)
        os.close(devnull)


# ---------------------------------------------------------------------------
# Fingerprint cache
# ---------------------------------------------------------------------------

# Default directory (relative to where the CLI is executed) for caching the
# computed fingerprints so they don't have to be recomputed on every restart.
DEFAULT_CACHE_DIR = ".mol-search-service"

# Bump this whenever the fingerprint computation logic or FINGERPRINTS config
# changes in a way that invalidates previously cached fingerprints.
CACHE_VERSION = 1


def _file_content_hash(path: str) -> str:
    """Stream a file through sha256 and return its hex digest."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _cache_path(cache_dir: str, compounds_file: str, fp_types: list[str]) -> str:
    """Build the cache file path keyed on file content, fp types and cache version."""
    content_hash = _file_content_hash(compounds_file)
    key_src = f"{content_hash}|{CACHE_VERSION}|{','.join(sorted(fp_types))}"
    key = hashlib.sha256(key_src.encode("utf-8")).hexdigest()
    return os.path.join(cache_dir, f"{key}.pkl")


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
    """A single compound's metadata."""

    id: str
    smiles: str
    db_name: str = "unknown"


@dataclass
class Dataset:
    """All precomputed data for one fingerprint type."""

    fps: list[Any]  # parallel list of raw RDKit FP objects


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
    # One entry per distinct matched fragment (deduplicated). Parallel lists:
    # matched_smiles[i] and matched_smarts[i] describe the same matched atoms,
    # rendered from the target molecule so its stereochemistry is preserved.
    matched_smiles: list[str] = field(default_factory=list)
    matched_smarts: list[str] = field(default_factory=list)


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
        python_method=rdFingerprintGenerator.GetAtomPairGenerator,
        default_options={"fpSize": 2048},
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
        python_method=rdFingerprintGenerator.GetTopologicalTorsionGenerator,
        default_options={"fpSize": 2048},
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

# Dynamically expand FINGERPRINTS to include explicit chiral variants right after their base versions
_ordered_fps = {}
for name, cfg in FINGERPRINTS.items():
    _ordered_fps[name] = cfg
    if cfg.stereo_options:
        chiral_cfg = replace(
            cfg,
            short_name=f"{cfg.short_name}_C",
            default_options={**cfg.default_options, **cfg.stereo_options},
            stereo_options={},
            description=f"{cfg.description} (Computed with stereochemistry enabled).",
        )
        _ordered_fps[f"{name}_chiral"] = chiral_cfg
FINGERPRINTS.clear()
FINGERPRINTS.update(_ordered_fps)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def get_fingerprint(
    mol: Chem.Mol, name: str = "morgan_ecfp"
) -> Any:
    """Generates a fingerprint for a molecule using the specified configuration name.

    Returns:
        An RDKit fingerprint object (`ExplicitBitVect` or `IntSparseIntVect`).
    """
    if name not in FINGERPRINTS:
        raise ValueError(f"Unknown fingerprint type: {name}")

    cfg = FINGERPRINTS[name]
    opts = cfg.default_options.copy()

    # Handle FCFP (feature invariants)
    if name.startswith("morgan_fcfp"):
        opts["atomInvariantsGenerator"] = (
            rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        )

    func: Callable[..., Any] = cfg.python_method

    if func.__name__.endswith("Generator"):
        generator = func(**opts)
        return generator.GetFingerprint(mol)
    else:
        # For others (including pattern), options are passed directly to the function along with mol
        # atom_pair and topological_torsion print C++-level deprecation spam on every call
        with _silence_stderr():
            return func(mol, **opts)


def safe_mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """
    Safely parses a SMILES string into an RDKit Mol object, avoiding strict
    sanitization bugs on organo-metalic compounds by carefully applying SanitizeFlags.

    Returns None on any parse/sanitization failure. RDKit's noisy C++-level
    warnings are suppressed; the caller decides how (and whether) to report the
    failure, so it can be summarized rather than printed per compound.
    """
    with _silence_stderr():
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                ^ Chem.SanitizeFlags.SANITIZE_CLEANUP,
            )
        except Exception:
            return None

        # Perceive stereochemistry from the parsed structure. With sanitize=False,
        # the SMILES directional bonds (/ and \) are stored but double-bond E/Z
        # stereo is never assigned (SanitizeMol does not do this step), leaving such
        # bonds as STEREONONE. Without this, use_chirality=True cannot distinguish
        # E/Z geometry. Guarded separately so a stereo-perception failure does not
        # discard an otherwise-valid molecule.
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass
    return mol


def safe_mol_from_smarts(smarts: str) -> Chem.Mol | None:
    """Parse a SMARTS pattern into an RDKit query Mol.

    Unlike SMILES, a SMARTS is parsed as a query graph: it keeps query
    features (wildcards, any-bond `~`, recursive SMARTS, degree/H/charge
    constraints) and is NOT sanitized or re-aromatized. Matching therefore
    follows the pattern literally — e.g. an aliphatic `C` will not match an
    aromatic carbon. Returns None on a parse failure.
    """
    with _silence_stderr():
        return Chem.MolFromSmarts(smarts)


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------


class MolSearchEngine:
    def __init__(self) -> None:
        # keys: fingerprint type name (e.g. 'morgan_ecfp', 'pattern')
        # values: compiled Dataset for that fingerprint type
        self.datasets: dict[str, Dataset] = {}
        self.core_data: list[CompoundEntry] = []
        self.db_indices: dict[str, list[int]] = {}

    def _get_indices(
        self, fp_type: str, db_names: list[str] | None = None
    ) -> range | list[int]:
        """Helper to get valid indices based on db_names filter."""
        dataset = self.datasets.get(fp_type)
        if not dataset:
            raise ValueError(f"Dataset {fp_type} not loaded.")

        if not db_names:
            # Return range covering all indices
            return range(len(self.core_data))

        indices: list[int] = []
        for db in db_names:
            if db in self.db_indices:
                indices.extend(self.db_indices[db])
        return indices

    def get_databases(self) -> list[str]:
        """Return a list of all database names present in the engine."""
        return list(self.db_indices.keys())

    def search_similarity(
        self,
        query_smiles: str,
        limit: int = 5,
        db_names: list[str] | None = None,
        fp_type: str = "morgan_ecfp",
        min_score: float = 0.0,
    ) -> list[SimilarityResult]:
        """Executes a similarity search and returns the top hits."""
        if fp_type not in self.datasets:
            raise ValueError(f"Error: Dataset {fp_type} not loaded.")

        dataset = self.datasets[fp_type]

        query_mol = safe_mol_from_smiles(query_smiles)
        if not query_mol:
            print(f"Error: Invalid query SMILES: {query_smiles}")
            return []

        query_fp = get_fingerprint(query_mol, fp_type)

        # Get candidate indices
        indices = self._get_indices(fp_type, db_names)

        # Retrieve FPs for the target indices
        if not db_names:
            target_fps: list[Any] = dataset.fps
            target_data: list[CompoundEntry] = self.core_data
        else:
            # Construct subset list
            target_fps = [dataset.fps[i] for i in indices]
            target_data = [self.core_data[i] for i in indices]

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
        query: str,
        limit: int = 5,
        db_names: list[str] | None = None,
        fp_type: str = "pattern",
        min_match_count: int = 1,
        use_chirality: bool = False,
        query_type: str = "smiles",
    ) -> list[SubstructureResult]:
        """Executes a substructure search (Screening + Verification).

        Args:
            query: The query pattern, interpreted according to ``query_type``.
            query_type: ``"smiles"`` (default) parses ``query`` with
                ``MolFromSmiles`` — sanitized, aromatized, stereo perceived.
                ``"smarts"`` parses it with ``MolFromSmarts`` as a query graph,
                keeping query features and matching literally (no
                re-aromatization).
            use_chirality: If True, both tetrahedral (R/S) and double-bond (E/Z)
                stereochemistry are taken into account during matching.
                Defaults to False for maximum recall. Most meaningful for SMILES
                queries; SMARTS encodes its own stereo constraints in the pattern.
        """
        if fp_type not in self.datasets:
            raise ValueError(f"Error: Dataset {fp_type} not loaded.")

        if query_type not in ("smiles", "smarts"):
            raise ValueError(
                f"Error: query_type must be 'smiles' or 'smarts', got {query_type!r}."
            )

        dataset = self.datasets[fp_type]

        if query_type == "smarts":
            query_mol = safe_mol_from_smarts(query)
        else:
            query_mol = safe_mol_from_smiles(query)
        if not query_mol:
            print(f"Error: Invalid query {query_type.upper()}: {query}")
            return []

        query_fp = get_fingerprint(query_mol, fp_type)

        # Get candidate indices
        indices = self._get_indices(fp_type, db_names)

        candidates_data: list[tuple[CompoundEntry, Any]] = []

        # Screening
        # Optimization: We avoid creating intermediate lists of ALL compounds
        fps = dataset.fps

        for i in indices:
            if fps[i] is None:
                continue
            if DataStructs.AllProbeBitsMatch(query_fp, fps[i]):
                candidates_data.append((self.core_data[i], fps[i]))

        # Verification
        results: list[SubstructureResult] = []
        for entry, fp in candidates_data:
            try:
                # Optimized verification: Parsing SMILES is the slow part.
                target_mol = safe_mol_from_smiles(entry.smiles)
                if target_mol is None:
                    continue
                matches = target_mol.GetSubstructMatches(
                    query_mol, useChirality=use_chirality
                )

                if matches and len(matches) >= min_match_count:
                    # Render each matched atom set back to SMILES and SMARTS from
                    # the target molecule, so the result carries the target's
                    # stereochemistry even when matching ignored it. Deduplicate
                    # fragments that render identically.
                    seen: set[tuple[str, str]] = set()
                    matched_smiles: list[str] = []
                    matched_smarts: list[str] = []
                    for atom_ids in matches:
                        atoms = list(atom_ids)
                        try:
                            smi = Chem.MolFragmentToSmiles(
                                target_mol, atomsToUse=atoms, isomericSmiles=True
                            )
                            sma = Chem.MolFragmentToSmarts(target_mol, atomsToUse=atoms)
                        except Exception:
                            continue
                        key = (smi, sma)
                        if key in seen:
                            continue
                        seen.add(key)
                        matched_smiles.append(smi)
                        matched_smarts.append(sma)

                    results.append(
                        SubstructureResult(
                            id=entry.id,
                            smiles=entry.smiles,
                            db_name=entry.db_name,
                            fp=fp,
                            match_count=len(matches),
                            matched_smiles=matched_smiles,
                            matched_smarts=matched_smarts,
                        )
                    )

                    if limit and len(results) >= limit:
                        break
            except Exception:
                continue
        return results

    def load_from_sparql(
        self,
        endpoint: str,
        query: str,
        fp_types: list[str] | None = None,
        cache_dir: str | None = DEFAULT_CACHE_DIR,
    ) -> None:
        """
        Fetch compound data from a SPARQL endpoint, store it in a temp TSV file,
        and load it into the engine. The temp file path is stored in the
        COMPOUNDS_FILE environment variable so that additional Uvicorn workers can
        reload the same data on startup without re-querying the endpoint.
        """
        import requests

        print(f"Fetching data from {endpoint}...")
        try:
            with requests.post(
                endpoint,
                data={"query": query},
                headers={"Accept": "text/tab-separated-values"},
                stream=True
            ) as resp:
                resp.raise_for_status()

                # Stream directly to temp TSV file
                fd, temp_path = tempfile.mkstemp(suffix=".tsv")
                with os.fdopen(fd, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)
        except Exception as e:
            raise RuntimeError(f"Failed to fetch data from SPARQL endpoint: {e}") from e

        # Signal workers to delete the temp file on shutdown
        os.environ["DELETE_COMPOUNDS_FILE"] = "1"
        self.load_file(temp_path, fp_types=fp_types, cache_dir=cache_dir)

    def load_file(
        self,
        compounds_file: str,
        fp_types: list[str] | None = None,
        cache_dir: str | None = DEFAULT_CACHE_DIR,
    ) -> None:
        """
        Load compounds from a TSV file.
        Format (by column order):
        1. chem IRI (e.g. <http://...>)
        2. SMILES string
        3. db name (optional)

        If `cache_dir` is set, computed fingerprints are cached there (keyed on
        file content + fingerprint types) and reused on subsequent loads instead
        of being recomputed. Pass `cache_dir=None` to disable caching.
        """
        import re

        # Regex for basic validation
        # IRI: wrapped in <> or starting with http
        iri_regex = re.compile(r"^<[^>]+>$|^https?://[^\s]+$")

        # Persist the path so that additional Uvicorn workers can pick it up
        os.environ["COMPOUNDS_FILE"] = compounds_file
        print(f"Reading and parsing compounds from {compounds_file}...")

        target_fps = fp_types if fp_types else list(FINGERPRINTS.keys())
        invalid_fps = [fp for fp in target_fps if fp not in FINGERPRINTS]
        if invalid_fps:
            available = ", ".join(sorted(FINGERPRINTS.keys()))
            raise ValueError(
                f"Unknown fingerprint type(s): {', '.join(invalid_fps)}\n"
                f"Available types: {available}"
            )

        # 'pattern' is always required for substructure search — add it silently if missing
        if "pattern" not in target_fps:
            print("  - Note: 'pattern' fingerprint automatically added (required for substructure search).")
            target_fps = list(target_fps) + ["pattern"]

        valid_fps = target_fps

        # Try to load precomputed fingerprints from the on-disk cache
        cache_file: str | None = None
        if cache_dir:
            try:
                cache_file = _cache_path(cache_dir, compounds_file, valid_fps)
            except OSError as e:
                print(f"  - Warning: could not hash compounds file for caching: {e}")
                cache_file = None

            if cache_file and os.path.exists(cache_file):
                try:
                    with open(cache_file, "rb") as cf:
                        state = pickle.load(cf)
                    self.core_data = state["core_data"]
                    self.db_indices = state["db_indices"]
                    self.datasets = state["datasets"]
                    print(
                        f"Loaded {len(self.core_data)} compounds with precomputed "
                        f"fingerprints from cache: {cache_file}"
                    )
                    return
                except Exception as e:
                    print(f"  - Warning: failed to read fingerprint cache ({e}). Recomputing.")

        for fp_name in valid_fps:
            self.datasets[fp_name] = Dataset(fps=[])

        self.core_data = []
        self.db_indices = {}

        valid_count = 0
        skipped_count = 0

        with open(compounds_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")

            # Peek at first row to detect header
            try:
                first_row = next(reader)
            except StopIteration:
                return

            # Header detection: if first col starts with '?' or doesn't look like an IRI
            is_header = False
            if first_row and (first_row[0].startswith("?") or not iri_regex.match(first_row[0])):
                is_header = True
                print(f"  - Header detected and skipped: {first_row}")

            # Process rows
            current_row = first_row if not is_header else None

            while True:
                if current_row:
                    row = current_row
                    current_row = None
                else:
                    try:
                        row = next(reader)
                    except StopIteration:
                        break

                if len(row) < 2:
                    raise ValueError(
                        f"Row {reader.line_num} has insufficient columns (found {len(row)}, expected at least 2)"
                    )

                # 1. chem IRI
                cid_raw = row[0].strip()
                # 2. SMILES
                smiles_raw = row[1].strip().strip('"')
                # 3. db (optional)
                # We do NOT strip("<>") here to preserve the distinction between URIs and literals.
                db_raw = row[2].strip() if len(row) > 2 else "unknown"

                # Validate with regex
                if not iri_regex.match(cid_raw):
                    raise ValueError(
                        f"Invalid IRI format on row {reader.line_num}: '{cid_raw}'. IRIs must be wrapped in <> or start with http(s)://"
                    )

                cid = cid_raw.strip("<>")
                entry = CompoundEntry(id=cid, smiles=smiles_raw, db_name=db_raw)

                # Defer SMILES validation to RDKit. Compounds that fail to parse
                # are skipped (and thus never enter the cache); we only count them
                # here and report a single summary line at the end to avoid flooding
                # the output with one warning per bad compound.
                mol = safe_mol_from_smiles(smiles_raw)
                if mol is None:
                    skipped_count += 1
                    continue

                # Append metadata
                idx = len(self.core_data)
                self.core_data.append(entry)

                if entry.db_name not in self.db_indices:
                    self.db_indices[entry.db_name] = []
                self.db_indices[entry.db_name].append(idx)

                # Compute and append fingerprints immediately to avoid holding mol in memory
                for fp_name in valid_fps:
                    try:
                        fp = get_fingerprint(mol, fp_name)
                        self.datasets[fp_name].fps.append(fp)
                    except Exception as e:
                        print(f"  - Warning: Failed to compute {fp_name} for {cid}: {e}")
                        self.datasets[fp_name].fps.append(None)

                valid_count += 1
                if valid_count % 100000 == 0:
                    print(f"  - Processed {valid_count} valid compounds...")

        if skipped_count:
            print(f"  - Skipped {skipped_count} compound(s) that RDKit could not parse.")
        print(f"Compilation complete. {valid_count} compounds loaded into the engine.")

        # Persist computed fingerprints to the cache for faster future restarts
        if cache_file and cache_dir:
            try:
                os.makedirs(cache_dir, exist_ok=True)
                state = {
                    "version": CACHE_VERSION,
                    "fp_types": sorted(valid_fps),
                    "core_data": self.core_data,
                    "db_indices": self.db_indices,
                    "datasets": self.datasets,
                }
                # Write to a temp file then atomically replace to avoid corrupt caches
                tmp_path = f"{cache_file}.tmp"
                with open(tmp_path, "wb") as cf:
                    pickle.dump(state, cf, protocol=pickle.HIGHEST_PROTOCOL)
                os.replace(tmp_path, cache_file)
                print(f"Cached fingerprints to {cache_file}")
            except Exception as e:
                print(f"  - Warning: failed to write fingerprint cache: {e}")


# Initialize Engine globally for easy sharing
engine = MolSearchEngine()
