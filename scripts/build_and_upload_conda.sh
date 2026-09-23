#!/bin/bash
# Build the hmmploidy conda package and upload it to anaconda.org under the
# `samuele.soraggi` user channel.
#
# Prerequisites (one-off):
#   conda install -n base -c conda-forge conda-build anaconda-client
#   anaconda login             # authenticate once, stores a token locally
#
# Usage:
#   ./scripts/build_and_upload_conda.sh            # build + upload
#   ./scripts/build_and_upload_conda.sh --build-only
set -euo pipefail

ANACONDA_USER="samuele.soraggi"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RECIPE_DIR="${REPO_ROOT}/recipe"

command -v conda-build >/dev/null 2>&1 || {
    echo "conda-build not found. Install with: conda install -n base -c conda-forge conda-build" >&2
    exit 1
}

echo "==> Building hmmploidy from ${RECIPE_DIR}"
OUTPUT_PATH="$(conda-build "${RECIPE_DIR}" -c conda-forge --output)"
conda-build "${RECIPE_DIR}" -c conda-forge

echo "==> Built package: ${OUTPUT_PATH}"

if [[ "${1:-}" == "--build-only" ]]; then
    echo "Skipping upload (--build-only). Upload later with:"
    echo "  anaconda upload -u ${ANACONDA_USER} ${OUTPUT_PATH}"
    exit 0
fi

command -v anaconda >/dev/null 2>&1 || {
    echo "anaconda-client not found. Install with: conda install -n base -c conda-forge anaconda-client" >&2
    exit 1
}

echo "==> Uploading to anaconda.org/${ANACONDA_USER}"
anaconda upload -u "${ANACONDA_USER}" "${OUTPUT_PATH}"
