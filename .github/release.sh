#!/bin/bash

# Validate version bump type
BUMP_TYPE=$1
if [ -z "$BUMP_TYPE" ] || { [ "$BUMP_TYPE" != "fix" ] && [ "$BUMP_TYPE" != "minor" ] && [ "$BUMP_TYPE" != "major" ]; }; then
    echo "Error: Version bump type is required and must be one of: fix, minor, major"
    echo "Usage: $0 <fix|minor|major>"
    exit 1
fi

uvx hatch version $BUMP_TYPE

VERSION=$(uvx hatch version)

# https://git-cliff.org/docs/usage/initializing/
uvx git-cliff --config unconventional -o CHANGELOG.md --tag v$VERSION
git add CHANGELOG.md src/*/__init__.py
git commit -m "chore: bump to v$VERSION"
git push
git tag -a "v$VERSION" -m "Release v$VERSION"
git push origin "v$VERSION"
