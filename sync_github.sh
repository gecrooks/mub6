#!/bin/sh
# Sync MUB6 commits from the ops repo to github.com/gecrooks/mub6.
# Run from anywhere. The subtree split is deterministic, so this is
# idempotent and safe to re-run.
set -e
cd "$(dirname "$0")/.."
git subtree split --prefix=MUB6 -b mub6-only >/dev/null
cd "$HOME/Work/mub6"
git pull -q "$HOME/Work/ops" mub6-only
git push -q origin main
echo "synced: $(git log --oneline -1)"
