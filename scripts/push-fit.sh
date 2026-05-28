#!/usr/bin/env bash

set -euo pipefail

PERSONAL_BRANCH=$(git branch --show-current)
FIT_EMAIL="mberkland2023@my.fit.edu"
FIT_NAME="Maverick Berkland"
FIT_BRANCH="sync-${PERSONAL_BRANCH}"
TMP_BRANCH="fit-${PERSONAL_BRANCH}"

git fetch origin
git fetch fit

git checkout -B "${TMP_BRANCH}" fit/main

COMMITS=$(git rev-list --reverse "fit/main..${PERSONAL_BRANCH}")

if [[ -z "${COMMITS}" ]]; then
    echo "No commits to sync."
    git checkout "${PERSONAL_BRANCH}"
    exit 0
fi

for COMMIT in ${COMMITS}; do
    GIT_AUTHOR_NAME="${FIT_NAME}" \
    GIT_AUTHOR_EMAIL="${FIT_EMAIL}" \
    GIT_COMMITTER_NAME="${FIT_NAME}" \
    GIT_COMMITTER_EMAIL="${FIT_EMAIL}" \
    git cherry-pick "${COMMIT}"
done

git push fit "${TMP_BRANCH}:${FIT_BRANCH}" --force-with-lease

git checkout "${PERSONAL_BRANCH}"

echo
echo "Done."
echo "GitHub.com branch: ${PERSONAL_BRANCH}"
echo "FIT sync branch: ${FIT_BRANCH}"