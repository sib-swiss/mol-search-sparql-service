import pytest
import subprocess
import tempfile
import os
import signal


def run_cli(*args):
    """Run the CLI, wait for it to exit, and return the result.
    Only use for commands that are expected to fail fast (validation errors).
    """
    result = subprocess.run(
        ["uv", "run", "mol-search-sparql-service"] + list(args),
        capture_output=True,
        text=True,
    )
    return result


def run_cli_validation_only(*args, timeout: float = 5.0):
    """Run the CLI and kill it after `timeout` seconds.

    Validation errors are printed and the process exits before data loading,
    so a short timeout is sufficient to capture them. If the process is still
    alive at the deadline it means validation passed — we kill it and return
    whatever stderr was collected so far.

    Returns (returncode_or_None, stderr).
    """
    # `uv run` launches the actual server as a grandchild process. Killing only
    # the `uv` parent (proc.kill()) would orphan the server, which keeps the
    # inherited stdout/stderr pipes open and makes the follow-up communicate()
    # block forever. start_new_session=True puts the whole tree in its own
    # process group so we can signal every descendant at once.
    proc = subprocess.Popen(
        ["uv", "run", "mol-search-sparql-service"] + list(args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )

    def _kill_process_group() -> None:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            # Process group already gone; nothing to clean up.
            pass

    try:
        _, stderr = proc.communicate(timeout=timeout)
        return proc.returncode, stderr
    except subprocess.TimeoutExpired:
        # Validation passed (process is still alive). Kill the entire group so
        # the orphaned server releases the pipes and communicate() can return.
        _kill_process_group()
        _, stderr = proc.communicate()
        return None, stderr  # None means "still running at timeout = validation passed"


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
    """Valid fingerprint subset must not produce a fingerprint validation error."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("?chem\t?smiles\t?db\n")
        f.write("<http://example.com/mol1>\tCCO\texample\n")
        temp_file = f.name

    try:
        rc, stderr = run_cli_validation_only("-f", temp_file, "-t", "morgan_ecfp,pattern")
        assert "Unknown fingerprint type" not in stderr
    finally:
        os.unlink(temp_file)


def test_valid_port_range():
    """Valid port numbers must not produce a port-range validation error."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("?chem\t?smiles\t?db\n")
        f.write("<http://example.com/mol1>\tCCO\texample\n")
        temp_file = f.name

    try:
        for port in [1, 80, 8080, 65535]:
            rc, stderr = run_cli_validation_only("-f", temp_file, "-p", str(port))
            assert "Port must be between 1 and 65535" not in stderr, (
                f"Port {port} should be in valid range but got: {stderr}"
            )
    finally:
        os.unlink(temp_file)
