#!/usr/bin/env bash
# Build the web app on Vercel.
#
# The app needs web/public/data (JSON per antigen chain, the phylogeny) and
# web/public/structures (the .cif.gz files the viewer loads). Neither is in git:
# both are generated from results/ and structures/, which are. So the build
# regenerates them, and a deploy always reflects the committed results rather
# than whatever someone last exported by hand.
#
# Needs python3 with venv (the Vercel build image has it) and node.
set -euo pipefail
cd "$(dirname "$0")/.."

python3 -m venv .build-venv
# shellcheck disable=SC1091
. .build-venv/bin/activate
pip install --quiet --disable-pip-version-check pandas numpy biopython

python scripts/08_export_web_data.py
python scripts/09_phylogeny.py

cd web
npm run build
