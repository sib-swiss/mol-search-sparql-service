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
@pytest.mark.parametrize("use_chirality", [True, False])
@pytest.mark.parametrize("limit", [1, 5])
@pytest.mark.parametrize("min_score", [0.0, 0.5])
def test_similarity_search_parameters(fp_type, use_chirality, limit, min_score):
    query_smiles = "[NH3+][C@@H](C)C(=O)[O-]" # Alanine
    results = engine.search_similarity(
        query_smiles,
        limit=limit,
        fp_type=fp_type,
        use_chirality=use_chirality,
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
        # NOTE: We only expect 1.0 if use_chirality matches the database (which is False)
        if r.compound.smiles == query_smiles:
            if not use_chirality:
                assert r.similarity == pytest.approx(1.0, rel=1e-5)
            else:
                # If chirality is enabled for query but not database, score will be < 1.0
                # BUT only for fingerprint types that actually use the chirality flag (Morgan)
                if fp_type in ["morgan_ecfp", "morgan_fcfp"]:
                    assert r.similarity < 1.0
                else:
                    # Others don't use the flag, so they still match 1.0
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
@pytest.mark.parametrize("use_chirality", [True, False])
@pytest.mark.parametrize("limit", [1, 2])
def test_substructure_search_parameters(use_chirality, limit):
    query_smarts = "C" # Any carbon
    results = engine.search_substructure(
        query_smarts,
        limit=limit,
        use_chirality=use_chirality
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
