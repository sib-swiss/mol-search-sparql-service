"""Documentation generator for the Chemistry Search SPARQL Service.

Generates markdown docs from function signatures and docstrings to
automatically update the README.md.

Usage:
    uv run src/mol_search_sparql_service/gen_docs.py
"""

import re
from pathlib import Path

from mol_search_sparql_service.sparql_service import generate_docs

# Resolve repo root relative to this script so it works regardless of cwd.
REPO_ROOT = Path(__file__).parent.parent.parent
README = REPO_ROOT / "README.md"

START_MARKER = "<!-- AUTOGEN_DOCS_START -->"
END_MARKER = "<!-- AUTOGEN_DOCS_END -->"


def update_readme() -> None:
    """Inject the generated docs into the README between the marker comments."""
    docs = generate_docs()
    new_block = f"{START_MARKER}\n{docs}\n{END_MARKER}"

    text = README.read_text()
    pattern = re.compile(
        re.escape(START_MARKER) + r".*?" + re.escape(END_MARKER),
        re.DOTALL,
    )
    if not pattern.search(text):
        raise ValueError(
            f"Could not find {START_MARKER!r} and {END_MARKER!r} in {README}. "
            "Please ensure these markers are present to indicate where the docs should be injected."
        )

    README.write_text(pattern.sub(new_block, text))
    print(f"✨ Updated {README.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    update_readme()
