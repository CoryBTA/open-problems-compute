#!/usr/bin/env bash
# One-time setup of the always-on worker on a fresh Ubuntu VM (Oracle Cloud Always Free ARM works; any Linux with python3.11+).
# Usage: bash setup-vm.sh   (run as the normal user; uses sudo for apt and systemd)
set -euo pipefail
REPO_SSH="git@github.com:CoryBTA/open-problems-compute.git"
DIR="$HOME/millennium-research"
sudo apt-get update -y && sudo apt-get install -y git python3-venv python3-dev build-essential
if [ ! -f "$HOME/.ssh/millennium_deploy" ]; then
  ssh-keygen -t ed25519 -N "" -C "millennium-vm-worker" -f "$HOME/.ssh/millennium_deploy"
  printf 'Host github.com\n  IdentityFile ~/.ssh/millennium_deploy\n  IdentitiesOnly yes\n' >> "$HOME/.ssh/config"
  chmod 600 "$HOME/.ssh/config"
  echo; echo "Add this key at https://github.com/CoryBTA/open-problems-compute/settings/keys (tick 'Allow write access'), then rerun this script:"; echo
  cat "$HOME/.ssh/millennium_deploy.pub"; exit 0
fi
ssh-keyscan github.com >> "$HOME/.ssh/known_hosts" 2>/dev/null
[ -d "$DIR/.git" ] || git clone --branch cloud "$REPO_SSH" "$DIR"
cd "$DIR"
python3 -m venv .venv
.venv/bin/pip install -q -r requirements.txt -r requirements-extra.txt
sed -e "s#__USER__#$USER#" -e "s#__REPO__#$DIR#g" scripts/vm/millennium-worker.service | sudo tee /etc/systemd/system/millennium-worker.service >/dev/null
sudo systemctl daemon-reload && sudo systemctl enable --now millennium-worker
echo "Started. Follow it with: journalctl -u millennium-worker -f"
