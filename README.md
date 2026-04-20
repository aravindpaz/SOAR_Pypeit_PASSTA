# SOAR PypeIt PASSTA

Python tools for **PASSTA** (PI: Hosseinzadeh) quicklook spectroscopy on **SOAR/Goodman**: pulling frames from the Las Cumbres archive, organizing raw FITS, running **PypeIt**, and optionally classifying spectra with **SNID SAGE** and posting summaries to **Slack**.

This repository is a companion to general SOAR/PypeIt examples such as [soar_pypeit_example](https://github.com/charliekilpatrick/soar_pypeit_example).

## Features

- **Archive download** (`goodman-spec-download`): Query the LCO archive for PASSTA Goodman science, standards, lamp flats, and arcs; write them under a dated `rawdata` tree. Shared calibrators are downloaded once per night; frame fetches can use multiple threads.
- **PypeIt pipeline** (`goodman-pypeit-reduce`): Symlink data per grating mode and object, run `pypeit_setup` / `run_pypeit`, build or reuse a sensitivity-function database, then flux-calibrate, coadd, and telluric-correct. Optional parallel workers speed independent object reductions and SNID runs.
- **Local file prep** (`goodman-file-organize`): Copy FITS from a staging directory into the same dated layout and run `fpack` where appropriate.
- **LCO client** (`passta.lcogt`): Reusable `LCOGT` class for archive queries, downloads, and (broader) LCO observation-request construction.

## Requirements

- **Python** 3.10+ (3.11+ recommended for current **PypeIt** wheels)
- **NumPy** 2.4.1 or newer is required by recent **PypeIt** releases; if you also use **numba**, ensure your environment resolves both constraints (for example use a dedicated venv for PASSTA reductions).
- **Environment variables** (for archive download): `LCOGT_USERNAME`, `LCOGT_PASSWORD`
- **External tools** on `PATH` where used:
  - **CFITSIO** `fpack` / `funpack` (file organize / some LCO downloads)
  - **PypeIt** CLIs: `pypeit_setup`, `run_pypeit`, `pypeit_sensfunc`, `pypeit_flux_calib`, `pypeit_coadd_1dspec`, `pypeit_tellfit`
  - Optional: if **`pypeit_tellfit`** times out fetching PCA telluric data from object storage, run **`scripts/install_telluric_pca_from_gdrive.sh`** (after `conda activate …`) to pull **`TellPCA_3000_26000_R15000.fits`** from the PypeIt Telluric [Google Drive](https://drive.google.com/drive/folders/1FFRWjUZ58HiDuDD33MYqBzMWDQanBRRy), **`pypeit_install_telluric --local_file`** it into the cache, then remove the download under **`WORKDIR`** (default repo root) so only the cached copy remains.
  - Optional: **`sage`** (SNID SAGE) for `--snid`
  - Optional: **Slack** bot token file for `--slack` (see below)

## Installation

### Conda (recommended)

From the **repository root** (the directory that contains `environment.yml` and `pyproject.toml`):

```bash
cd /path/to/SOAR_Pypeit_PASSTA
conda env create -f environment.yml
conda activate passta
```

The file **`environment.yml`** creates the **`passta`** environment on **`conda-forge`** (Python and **pip** only), then uses **pip** to install this package in **editable** mode (`pip install -e .`) so versions match **`pyproject.toml`** (including **PypeIt** and the other Python dependencies). **`fpack`** / **`funpack`** are not installed by this file; add them separately (system packages, a HEASARC **CFITSIO** build on **`PATH`**, or another channel recipe that ships the binaries). Version metadata uses **setuptools-scm**, so keep the **`.git`** directory when cloning.

To refresh the install after `git pull`:

```bash
conda activate passta
python -m pip install -U pip
python -m pip install -U -e .
```

For a **non-editable** install inside the same environment:

```bash
conda activate passta
python -m pip install .
```

### pip and venv only

If you prefer not to use conda, create a venv and install from the repo root:

```bash
cd /path/to/SOAR_Pypeit_PASSTA
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install -U pip
python -m pip install -e .
```

Install **`cfitsio`** / **`fpack`** separately (for example from your OS package manager or a conda-only throwaway env) if you need **`goodman-file-organize`** or archive paths that call **`funpack`**.

### Entry points

After installation, three console scripts are available:

| Script | Module |
|--------|--------|
| `goodman-spec-download` | `passta.spec_download` |
| `goodman-pypeit-reduce` | `passta.pypeit_reduce` |
| `goodman-file-organize` | `passta.file_organize` |

Import the library in Python:

```python
from passta import LCOGT, __version__
```

## Usage

### 1. Download data for a night

```bash
export LCOGT_USERNAME=... LCOGT_PASSWORD=...
goodman-spec-download 2025-04-15 --propid SOAR2025B-004 --outdir /path/to/reductions
```

Use a **comma-separated** `--propid` list so a single process handles every program for that night (fewer archive queries and one pass over shared standards/flats):

```bash
goodman-spec-download 2025-04-15 \
  --propid SOAR2025B-004,SOAR2026A-002 \
  --outdir /path/to/reductions
```

Parallel HTTP downloads (bounded by the archive; default is a few threads) are controlled with `--download-workers` (set to `1` to force strictly sequential downloads).

This creates `/path/to/reductions/20250415/rawdata` (and an empty `workspace` tree) populated from the archive.

### 2. Run PypeIt and post-processing

```bash
goodman-pypeit-reduce 2025-04-15 /path/to/reductions
```

Optional:

```bash
goodman-pypeit-reduce 2025-04-15 /path/to/reductions --snid --slack
```

```bash
goodman-pypeit-reduce 2025-04-15 /path/to/reductions --slack-dry-run
```

- **`--snid`**: Runs `sage` on each telluric-corrected coadd under the reduction tree.
- **`--slack`**: Posts to Slack when combined with `--snid` (requires `SLACK_TOKEN_FILE_PASSTA` pointing to a file containing the bot token).
- **`--slack-dry-run`**: Runs SNID and builds the same Slack message as a live post, writes `spectra/slack/*.slack.debug`, and logs the text; does not call the Slack API (implies `--snid`). If both `--slack` and `--slack-dry-run` are passed, only the dry-run runs.
- **`--jobs N`**: Run `pypeit_setup`, `run_pypeit`, and SNID in parallel across independent object directories (default `1`). Try `2`–`4` on a multi-core machine; avoid oversubscribing shared NFS storage.

Astropy’s observatory site cache is refreshed **once** at the start of each reduction run (no longer on every import or every object directory).

### Logging

All three CLIs accept **`--log-level`** (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`; default `INFO`) and **`-v` / `--verbose`** (sets `DEBUG`). Logs go to **stderr** with timestamps and logger names, for example:

```bash
goodman-pypeit-reduce 2025-04-15 /path/to/reductions --log-level WARNING
goodman-spec-download 2025-04-15 --outdir /path/to/reductions -v
```

When importing `passta` as a library, configure logging in your application (the package does not call `basicConfig` except from these CLI entry points).

Sensitivity functions shipped with the package live under the installed `passta/caldb/` directory (`caldb.txt` plus FITS). To use your own database directory:

```bash
goodman-pypeit-reduce 2025-04-15 /path/to/reductions --caldb /path/to/my_caldb
```

### 3. Organize pre-downloaded FITS

```bash
goodman-file-organize 2025-04-15 /path/to/local_fits --outdir /path/to/reductions
```

## Performance and cron

Typical end-to-end time is dominated by **PypeIt** (and **SNID** if enabled), not Python overhead. Practical wins:

1. **One download command** with all proposal IDs (`--propid A,B`) instead of two sequential `timeout 600` runs.
2. **`--download-workers`** (defaults to a small positive number) for overlapping LCO frame downloads; lower it if the archive returns HTTP 429 throttling.
3. **`goodman-pypeit-reduce --jobs 2`** (or `3`–`4`) when you have several object directories for the same night—each reduction is independent until shared `caldb` updates (sensitivity functions still run sequentially inside each reduction as today).
4. **Timeouts**: a single longer download timeout may replace two shorter ones; keep headroom for large nights.
5. **Skip SNID** on the fastest path (`--snid` off) when classifications are not needed that morning.

Example (Linux server, conda env named `pypeit`, installed entry points on `PATH`):

```bash
date=$(date -d "16 hours ago" +%Y-%m-%d)
source "$CONDA_SRC/bin/activate" pypeit
timeout 900 goodman-spec-download "$date" \
  --outdir /data/ckilpatrick/PASSTA \
  --propid SOAR2025B-004,SOAR2026A-002
timeout 3600 goodman-pypeit-reduce "$date" /data/ckilpatrick/PASSTA \
  --snid --slack --jobs 3
```

Further gains require **PypeIt** tuning (fewer calibration steps, smaller extraction windows, GPU builds if available) or running reductions on local SSD rather than network disks.

## Layout conventions

Under a chosen `--outdir` / `outdir`:

- **`YYYYMMDD/rawdata`**: Raw or fpacked Goodman spectra and calibrations.
- **`YYYYMMDD/workspace`**: Per-grating, per-object symlink trees used for PypeIt.
- **`spectra/`**: Symlinks to final `_tellcorr.fits` products from `goodman-pypeit-reduce`.

## Repository layout

```text
passta/
  __init__.py          # Package version and exports
  lcogt.py             # LCO archive / request client (class LCOGT)
  spec_download.py     # PASSTA Goodman archive download CLI
  pypeit_reduce.py     # PypeIt driver and post-processing CLI
  file_organize.py     # Local FITS organize + fpack CLI
  slack_utils.py       # Slack upload helpers
  caldb/
    caldb.txt          # Shipped sensitivity-function index (ECSV)
    *.fits             # Sensitivity functions listed in caldb.txt
pyproject.toml
README.md
```

## License and attribution

Code authorship is noted in module docstrings (notably the LCO client). Set `license` in `pyproject.toml` to match your distribution policy if you publish the package.
