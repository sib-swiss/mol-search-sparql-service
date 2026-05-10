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

