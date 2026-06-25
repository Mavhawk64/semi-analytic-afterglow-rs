#!/usr/bin/env bash

set -euo pipefail

PERSONAL_BRANCH=$(git branch --show-current)
FIT_EMAIL="mberkland2023@my.fit.edu"
FIT_NAME="Maverick Berkland"
FIT_BRANCH="sync-${PERSONAL_BRANCH}"
TMP_BRANCH="fit-${PERSONAL_BRANCH}"

if [[ "${PERSONAL_BRANCH}" == fit-* || "${PERSONAL_BRANCH}" == sync-* ]]; then
    echo "Do not run this from a FIT sync branch."
    exit 1
fi

git fetch origin
git fetch fit

git checkout -B "${TMP_BRANCH}" fit/main

COMMITS=$(git rev-list --reverse --no-merges "origin/main..${PERSONAL_BRANCH}")

if [[ -z "${COMMITS}" ]]; then
    echo "No commits to sync."
    git checkout "${PERSONAL_BRANCH}"
    exit 0
fi

for COMMIT in ${COMMITS}; do
    echo "Cherry-picking ${COMMIT}..."

    set +e
    GIT_AUTHOR_NAME="${FIT_NAME}" \
    GIT_AUTHOR_EMAIL="${FIT_EMAIL}" \
    GIT_COMMITTER_NAME="${FIT_NAME}" \
    GIT_COMMITTER_EMAIL="${FIT_EMAIL}" \
    git cherry-pick "${COMMIT}"
    STATUS=$?
    set -e

    if [[ ${STATUS} -ne 0 ]]; then
        if git diff --quiet && git diff --cached --quiet; then
            echo "Skipping empty cherry-pick: ${COMMIT}"
            git cherry-pick --skip
        else
            echo "Cherry-pick conflict at ${COMMIT}."
            exit 1
        fi
    fi
done

git push fit "${TMP_BRANCH}:${FIT_BRANCH}" --force-with-lease

git checkout "${PERSONAL_BRANCH}"

echo
echo "Done."
echo "GitHub.com branch: ${PERSONAL_BRANCH}"
echo "FIT sync branch: ${FIT_BRANCH}"
