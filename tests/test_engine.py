import pytest
from mol_search_sparql_service.rdkit_fingerprints import SearchEngine

@pytest.fixture(scope="module")
def engine():
    eng = SearchEngine()
    eng.load_and_compile("compounds.tsv")
    return eng

def test_engine_similarity_search(engine):
    test_mol = "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]"

    # 1. Test defaults
    results = engine.search_similarity(test_mol)
    assert len(results) > 0
    assert "similarity" in results[0]

    # 2. Test limit
    results = engine.search_similarity(test_mol, limit=3)
    assert len(results) <= 3

    # 3. Test different fingerprint type
    results = engine.search_similarity(test_mol, fp_type="morgan_fcfp")
    assert len(results) > 0

    results = engine.search_similarity(test_mol, fp_type="topological_torsion")
    assert len(results) > 0


def test_engine_substructure_search(engine):
    smart = "c1ccccc1"

    # 1. Test defaults
    results = engine.search_substructure(smart)
    assert len(results) > 0
    assert "match_count" in results[0]

    # 2. Test limit
    results = engine.search_substructure(smart, limit=2)
    assert len(results) <= 2

def test_invalid_fingerprint(engine):
    with pytest.raises(ValueError):
        engine.search_similarity("C", fp_type="nonexistent_fp")
