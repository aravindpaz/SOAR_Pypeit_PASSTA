#!/usr/bin/env bash
# Download TellPCA_3000_26000_R15000.fits from the PypeIt Telluric Google Drive
# and register it with PypeIt (avoids S3 + short Astropy timeouts in pypeit_install_telluric).
#
# Default file:
#   https://drive.google.com/file/d/1DvIGYiFJG9oWkYDKRP8x8ronD6fkW9_c/view?usp=drive_link
#
# Requires: Python 3, network access. Uses ``gdown`` with ``uc?id=`` on older
# releases, or ``--fuzzy`` when the installed ``gdown`` supports it. Installs
# ``gdown`` via pip if missing. Requires ``pypeit_install_telluric`` on PATH.
# After a successful install, the download under WORKDIR is removed so only the
# PypeIt cache copy (telluric/atm_grids) remains.

set -euo pipefail

# Google Drive file id for TellPCA_3000_26000_R15000.fits (PypeIt Telluric folder)
GDRIVE_FILE_ID="${GDRIVE_FILE_ID:-1DvIGYiFJG9oWkYDKRP8x8ronD6fkW9_c}"
OUT_NAME="${OUT_NAME:-TellPCA_3000_26000_R15000.fits}"
WORKDIR="${WORKDIR:-.}"

usage() {
  cat <<'EOF'
Download TellPCA_3000_26000_R15000.fits from the PypeIt Telluric Google Drive,
run pypeit_install_telluric --local_file, then delete the file under WORKDIR so
only the PypeIt cache copy remains.

Usage:
  install_telluric_pca_from_gdrive.sh

Environment overrides:
  GDRIVE_FILE_ID   Drive file id (default: TellPCA_3000_26000_R15000)
  OUT_NAME         Save as this filename under WORKDIR
  WORKDIR          Download directory (default: .)
  PYTHON           Python for gdown / pip (default: python3)

Examples:
  conda activate passta
  ./scripts/install_telluric_pca_from_gdrive.sh

  WORKDIR=/tmp ./scripts/install_telluric_pca_from_gdrive.sh
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

PYTHON="${PYTHON:-python3}"

if ! command -v pypeit_install_telluric &>/dev/null; then
  echo "error: pypeit_install_telluric not on PATH (activate the env where PypeIt is installed)." >&2
  exit 1
fi

mkdir -p "$WORKDIR"
cd "$WORKDIR"
DEST="$(pwd)/${OUT_NAME}"

ensure_gdown() {
  if command -v gdown &>/dev/null; then
    echo "Using gdown: $(command -v gdown)"
    return 0
  fi
  if "$PYTHON" -m gdown --help &>/dev/null; then
    echo "Using: $PYTHON -m gdown"
    return 0
  fi
  echo "Installing gdown for: $PYTHON"
  "$PYTHON" -m pip install -q gdown
}

run_gdown() {
  # Older gdown releases do not support --fuzzy; uc?id= works for typical file sizes.
  local uc_url="https://drive.google.com/uc?id=${GDRIVE_FILE_ID}"
  local view_url="https://drive.google.com/file/d/${GDRIVE_FILE_ID}/view?usp=drive_link"
  if command -v gdown &>/dev/null; then
    if gdown -h 2>&1 | grep -q '[[:space:]]--fuzzy'; then
      gdown --fuzzy "$view_url" -O "$OUT_NAME"
    else
      gdown "$uc_url" -O "$OUT_NAME"
    fi
  else
    if "$PYTHON" -m gdown -h 2>&1 | grep -q '[[:space:]]--fuzzy'; then
      "$PYTHON" -m gdown --fuzzy "$view_url" -O "$OUT_NAME"
    else
      "$PYTHON" -m gdown "$uc_url" -O "$OUT_NAME"
    fi
  fi
}

ensure_gdown
echo "Downloading to: $DEST"
run_gdown

if [[ ! -f "$OUT_NAME" ]]; then
  echo "error: download did not create $OUT_NAME" >&2
  exit 1
fi

# Reject accidental HTML login / virus-scan pages (tiny non-FITS files)
size="$(wc -c < "$OUT_NAME" | awk '{print $1}')"
if [[ "$size" -lt 100000 ]]; then
  echo "error: $OUT_NAME is only ${size} bytes; expected a multi-MB FITS. Check Drive access (sign in / permissions)." >&2
  exit 1
fi

echo "Installing into PypeIt cache..."
pypeit_install_telluric --local_file "$DEST"
rm -f "$DEST"
echo "Removed download copy: $DEST"
echo "Done. Telluric PCA is only in the PypeIt cache; you can run pypeit_tellfit / goodman-pypeit-reduce without fetching from S3."
