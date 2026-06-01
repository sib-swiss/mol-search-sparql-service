import pytest
from mol_search_sparql_service.rdkit_fingerprints import engine, FINGERPRINTS

# Setup a small test engine
@pytest.fixture(scope="module", autouse=True)
def setup_engine():
    import tempfile
    import os
    
    # Create a small TSV for testing
    tsv_content = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/1>\tCCO\tdb1\n"
        "<http://ex.org/2>\tc1ccccc1\tdb1\n"
        "<http://ex.org/3>\tCC(=O)O\tdb2\n"
        "<http://ex.org/4>\t[NH3+][C@@H](C)C(=O)[O-]\tdb2\n" # Alanine
    )
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        tmp.write(tsv_content)
        tmp_path = tmp.name
        
    try:
        engine.load_file(tmp_path)
        yield
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

# Test combinations for similarity search
@pytest.mark.parametrize("fp_type", list(FINGERPRINTS.keys()))
@pytest.mark.parametrize("limit", [1, 5])
@pytest.mark.parametrize("min_score", [0.0, 0.5])
def test_similarity_search_parameters(fp_type, limit, min_score):
    query_smiles = "[NH3+][C@@H](C)C(=O)[O-]" # Alanine
    results = engine.search_similarity(
        query_smiles,
        limit=limit,
        fp_type=fp_type,
        min_score=min_score
    )
    
    # Validation
    assert isinstance(results, list)
    assert len(results) <= limit
    for r in results:
        assert r.similarity >= min_score
        assert r.similarity <= 1.0000001 # Floating point precision
        assert hasattr(r.compound, "id")
        assert hasattr(r.compound, "smiles")
        
        # Specific value validation for identical molecule
        # Since the engine now automatically builds chiral datasets,
        # querying the identical molecule should yield ~1.0 regardless of the use_chirality flag.
        if r.compound.smiles == query_smiles:
            assert r.similarity == pytest.approx(1.0, rel=1e-5)

# Test specific known similarity values
def test_specific_similarity_values():
    # Alanine vs Alanine
    alanine_smiles = "[NH3+][C@@H](C)C(=O)[O-]"
    results = engine.search_similarity(alanine_smiles, fp_type="morgan_ecfp")
    alanine_score = next((r.similarity for r in results if "4" in r.compound.id), 0.0) # ex.org/4 is Alanine
    
    assert alanine_score == pytest.approx(1.0)
    
    # Ethanol vs Ethanol
    results = engine.search_similarity("CCO", fp_type="morgan_ecfp")
    ethanol_score = next((r.similarity for r in results if "1" in r.compound.id), 0.0) # ex.org/1 is Ethanol
    assert ethanol_score == pytest.approx(1.0)

# Test combinations for substructure search
@pytest.mark.parametrize("limit", [1, 2])
def test_substructure_search_parameters(limit):
    query_smarts = "C" # Any carbon
    results = engine.search_substructure(
        query_smarts,
        limit=limit
    )
    
    # Validation
    assert isinstance(results, list)
    assert len(results) <= limit
    for r in results:
        assert r.match_count >= 1
        assert isinstance(r.id, str)
        assert isinstance(r.smiles, str)

# Test filtering by db_names
def test_db_filtering():
    # Filter for db1 (Ethanol and Benzene)
    results = engine.search_similarity("CCO", db_names=["db1"])
    for r in results:
        assert r.compound.db_name == "db1"
        
    # Filter for db2 (Acetic acid and Alanine)
    results = engine.search_similarity("CCO", db_names=["db2"])
    for r in results:
        assert r.compound.db_name == "db2"

# Test invalid parameters
def test_invalid_parameters():
    with pytest.raises(ValueError):
        engine.search_similarity("CCO", fp_type="invalid_fp")
        
    # Invalid SMILES should return empty list (graceful failure)
    results = engine.search_similarity("INVALID_SMILES")
    assert results == []

def test_get_databases():
    dbs = engine.get_databases()
    assert set(dbs) == {"db1", "db2"}


def test_substructure_chirality_tetrahedral():
    """use_chirality=True must distinguish R/S tetrahedral stereocenters."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/l-ala>\t[NH3+][C@@H](C)C(=O)[O-]\ttest\n"   # L-Ala (S)
        "<http://ex.org/d-ala>\t[NH3+][C@H](C)C(=O)[O-]\ttest\n"    # D-Ala (R)
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        query = "[NH3+][C@@H](C)C(=O)[O-]"  # L-Ala query

        # Without chirality: both enantiomers match
        results = eng.search_substructure(query, use_chirality=False)
        ids = {r.id for r in results}
        assert "http://ex.org/l-ala" in ids
        assert "http://ex.org/d-ala" in ids

        # With chirality: only L-Ala matches
        results = eng.search_substructure(query, use_chirality=True)
        ids = {r.id for r in results}
        assert "http://ex.org/l-ala" in ids
        assert "http://ex.org/d-ala" not in ids
    finally:
        os.unlink(path)


def test_substructure_chirality_ez():
    """use_chirality=True must distinguish E/Z double-bond geometry."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/trans>\tF/C=C/F\ttest\n"    # (E)
        "<http://ex.org/cis>\tF/C=C\\F\ttest\n"     # (Z)
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        query = "F/C=C/F"  # (E) query

        # Without chirality: both isomers match
        results = eng.search_substructure(query, use_chirality=False)
        ids = {r.id for r in results}
        assert "http://ex.org/trans" in ids
        assert "http://ex.org/cis" in ids

        # With chirality: only (E) matches
        results = eng.search_substructure(query, use_chirality=True)
        ids = {r.id for r in results}
        assert "http://ex.org/trans" in ids
        assert "http://ex.org/cis" not in ids
    finally:
        os.unlink(path)


def test_matched_fragments_populated_and_parallel():
    """Each result carries parallel matched_smiles / matched_smarts lists."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/butene>\tC/C=C/C\ttest\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        results = eng.search_substructure("C=CC", use_chirality=False)
        assert len(results) == 1
        r = results[0]
        # Lists are parallel and non-empty.
        assert len(r.matched_smiles) == len(r.matched_smarts)
        assert len(r.matched_smiles) >= 1
        # match_count counts raw matches; deduped fragments never exceed it.
        assert len(r.matched_smiles) <= r.match_count
        assert all(isinstance(s, str) and s for s in r.matched_smiles)
        assert all(isinstance(s, str) and s for s in r.matched_smarts)
    finally:
        os.unlink(path)


def test_matched_fragments_preserve_stereo_when_search_ignores_it():
    """A stereo-free query still returns matched fragments carrying the
    target's E/Z geometry."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/trans>\tC/C=C/C\ttest\n"   # (E)
        "<http://ex.org/cis>\tC/C=C\\C\ttest\n"    # (Z)
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        # Query has no stereo, so both isomers match...
        results = {r.id: r for r in eng.search_substructure("C=CC", use_chirality=False)}
        assert set(results) == {"http://ex.org/trans", "http://ex.org/cis"}

        # ...but each matched fragment keeps the directional bond from its target.
        trans_frags = results["http://ex.org/trans"].matched_smiles
        cis_frags = results["http://ex.org/cis"].matched_smiles
        assert any("/" in s and "\\" not in s for s in trans_frags)
        assert any("\\" in s for s in cis_frags)
        # SMARTS fragments also encode the geometry.
        assert any("/" in s or "\\" in s for s in results["http://ex.org/trans"].matched_smarts)
    finally:
        os.unlink(path)


def test_matched_fragments_preserve_tetrahedral_stereo():
    """A stereo-free query still returns matched fragments carrying the
    target's tetrahedral (R/S) configuration."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/l-ala>\t[NH3+][C@@H](C)C(=O)[O-]\ttest\n"  # L-Ala (S)
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        # Query has no stereo, but the fragment must keep the stereocenter.
        r = eng.search_substructure("NC(C)C(=O)O", use_chirality=False)[0]
        assert any("@" in s for s in r.matched_smiles)
        assert any("@" in s for s in r.matched_smarts)
    finally:
        os.unlink(path)


def test_matched_fragments_deduplicated():
    """Fragments that render identically are collapsed; distinct ones are kept."""
    import tempfile, os
    from mol_search_sparql_service.rdkit_fingerprints import MolSearchEngine

    # Bibenzyl has two benzene rings: a benzene query matches both, but the two
    # fragments render identically -> collapsed to a single deduped fragment.
    tsv = (
        "?chem\t?smiles\t?db\n"
        "<http://ex.org/bibenzyl>\tc1ccccc1CCc1ccccc1\ttest\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv)
        path = f.name

    try:
        eng = MolSearchEngine()
        eng.load_file(path)

        r = eng.search_substructure("c1ccccc1", use_chirality=False)[0]
        assert r.match_count == 2  # two rings matched
        assert len(r.matched_smiles) == 1  # both render the same -> deduped
        # No duplicates remain among the (smiles, smarts) pairs.
        pairs = list(zip(r.matched_smiles, r.matched_smarts))
        assert len(pairs) == len(set(pairs))
    finally:
        os.unlink(path)

