def test_import():
    try:
        import mol_search_sparql_service
        assert True
    except Exception as e:
        print(e)
        assert False, "Failed to import mol_search_sparql_service"
