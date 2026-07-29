#!/bin/sh
# Sync MUB6 commits from the ops repo into the local ~/Work/mub6
# mirror. Does NOT push — Gavin pushes manually (his rule: agents
# never change remote state). Run from anywhere; idempotent.
set -e
cd "$(dirname "$0")/.."
git subtree split --prefix=MUB6 -b mub6-only >/dev/null
cd "$HOME/Work/mub6"
git pull -q "$HOME/Work/ops" mub6-only
echo "local mirror updated: $(git log --oneline -1)"
echo "to publish, run:  cd ~/Work/mub6 && git push origin main"
