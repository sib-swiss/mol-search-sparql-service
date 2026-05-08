import pytest
import subprocess
import tempfile
import os
from pathlib import Path


def run_cli(*args):
    """Run the CLI and return the result."""
    result = subprocess.run(
        ["uv", "run", "mol-search-sparql-service"] + list(args),
        capture_output=True,
        text=True,
    )
    return result


def test_unknown_fingerprint_type():
    """Test that unknown fingerprint type produces informative error."""
    result = run_cli("-f", "compounds.tsv", "-t", "invalid_fingerprint_type")
    assert result.returncode != 0
    assert "Unknown fingerprint type" in result.stderr
    assert "invalid_fingerprint_type" in result.stderr
    # Should list available types
    assert "morgan_ecfp" in result.stderr
    assert "pattern" in result.stderr


def test_multiple_unknown_fingerprint_types():
    """Test error message lists all unknown fingerprint types."""
    result = run_cli(
        "-f",
        "compounds.tsv",
        "-t",
        "unknown1,unknown2,pattern,unknown3",
    )
    assert result.returncode != 0
    assert "Unknown fingerprint type" in result.stderr
    assert "unknown1" in result.stderr
    assert "unknown2" in result.stderr
    assert "unknown3" in result.stderr


def test_missing_compounds_file():
    """Test that missing compounds file produces informative error."""
    result = run_cli("-f", "nonexistent_compounds_file.tsv")
    assert result.returncode != 0
    assert "Compounds file not found" in result.stderr
    assert "nonexistent_compounds_file.tsv" in result.stderr


def test_missing_sparql_query_file():
    """Test that missing SPARQL query file produces informative error."""
    result = run_cli(
        "-s",
        "nonexistent_query.rq",
        "-e",
        "http://example.com/sparql",
    )
    assert result.returncode != 0
    assert "SPARQL query file not found" in result.stderr
    assert "nonexistent_query.rq" in result.stderr


def test_invalid_port_too_high():
    """Test that port number above 65535 produces error."""
    result = run_cli("-f", "compounds.tsv", "-p", "99999")
    assert result.returncode != 0
    assert "Port must be between 1 and 65535" in result.stderr
    assert "99999" in result.stderr


def test_invalid_port_too_low():
    """Test that port number below 1 produces error."""
    result = run_cli("-f", "compounds.tsv", "-p", "0")
    assert result.returncode != 0
    assert "Port must be between 1 and 65535" in result.stderr


def test_invalid_port_negative():
    """Test that negative port number produces error."""
    result = run_cli("-f", "compounds.tsv", "-p", "-1")
    assert result.returncode != 0
    assert "Port must be between 1 and 65535" in result.stderr


def test_valid_fingerprint_subset():
    """Test that valid fingerprint subset is accepted (doesn't error)."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as f:
        # Write minimal valid TSV
        f.write("?chem\t?smiles\t?db\n")
        f.write("<http://example.com/mol1>\tCCO\texample\n")
        temp_file = f.name

    try:
        # This should not error during CLI parsing (it will error later when trying to start server)
        result = run_cli("-f", temp_file, "-t", "morgan_ecfp,pattern")
        # We expect it to fail when trying to start the server (port issue), not CLI parsing
        # So we check that the error is NOT about fingerprints
        assert "Unknown fingerprint type" not in result.stderr
    finally:
        os.unlink(temp_file)


def test_valid_port_range():
    """Test that valid port numbers (1-65535) don't error in validation."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as f:
        # Write minimal valid TSV
        f.write("?chem\t?smiles\t?db\n")
        f.write("<http://example.com/mol1>\tCCO\texample\n")
        temp_file = f.name

    try:
        for port in [1, 80, 8080, 65535]:
            result = run_cli("-f", temp_file, "-p", str(port))
            # Should not error due to port validation
            assert "Port must be between 1 and 65535" not in result.stderr
    finally:
        os.unlink(temp_file)
