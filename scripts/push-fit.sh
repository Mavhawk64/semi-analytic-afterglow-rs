#!/usr/bin/env bash

set -e

BRANCH=$(git branch --show-current)

echo "Pushing to GitHub.com..."
git push origin "$BRANCH"

echo "Syncing to FIT..."

git push fit \
    "$BRANCH":sync-"$BRANCH" \
    --force-with-lease

echo
echo "Done."
echo
echo "GitHub.com branch:"
echo "  $BRANCH"
echo
echo "FIT sync branch:"
echo "  sync-$BRANCH"
