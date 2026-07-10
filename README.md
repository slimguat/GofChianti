# GofChianti

**Lightweight Python access to precomputed CHIANTI contribution functions
`G(T, nₑ)` — no IDL, no ChiantiPy at runtime.**

GofChianti gives you direct, offline-capable access to CHIANTI contribution
functions for selected spectral lines. The heavy computation is done once by the
maintainers using the official CHIANTI/SSW IDL routines; the results are stored
compactly and shipped as a small dataset. As a user you simply ask for a line,
and GofChianti downloads (and caches) the precomputed `G(T, nₑ)`, then lets you
interpolate it at any density/temperature — optionally scaled by any elemental
abundance set.

## Install

```bash
pip install gofchianti
```

## Quick start

```python
import astropy.units as u
import gofchianti as gc

# 1. What is available?
df = gc.available_lines()            # pandas DataFrame (ion, λ, version, f, A, ...)
df = gc.available_lines(version="11.0.2")

# 2. Load a line (downloaded + cached on first use, offline afterwards).
cf = gc.get_line("Fe_12", 195.119)

# 3. Evaluate the *bare* contribution function G(nₑ, T).
#    Inputs are astropy Quantities; the result is a Quantity too.
g = cf.get_gofnt(density=1e9 * u.cm**-3, temperature=1.5e6 * u.K)

# 3b. Give any two of (density, temperature, pressure). Pressure may be a
#     reduced pressure nₑ*T (cm^-3 K, CHIANTI-style) or a thermal pressure (Pa).
g = cf.get_gofnt(pressure=1e15 * u.cm**-3 * u.K, temperature=1.5e6 * u.K)
g = cf.get_gofnt(pressure=0.02 * u.Pa, density=1e9 * u.cm**-3)

# 4. Abundance-scaled values: multiply by Fe/H from a CHIANTI abundance set.
cf = gc.get_line("Fe_12", 195.119, abundance="sun_photospheric_2021_asplund")
g = cf.get_gofnt(density=1e9 * u.cm**-3, temperature=1.5e6 * u.K)  # now × Fe/H
```

## Abundances

`G` is stored **without** any elemental abundance. To get abundance-scaled
values, attach an abundance set:

```python
gc.available_abundances()                       # names shipped with the dataset
cf.set_abundance("sun_coronal_2021_chianti")    # by name (fetched/cached)
cf.set_abundance("/path/to/custom.abund")       # a local .abund file
cf.set_abundance(gc.Abundance.from_dex({"Fe": 8.0, "H": 12.0}))  # custom inline
cf.set_abundance(None)                          # back to bare G
```

## Offline use

```python
gc.download_all()        # fetch the entire dataset into the cache, once
```

After that everything works with no network access.

## Cache & configuration

Downloaded files live in an OS-standard cache directory
(`~/.cache/gofchianti` on Linux) and everything works offline after the first
download. Every setting is available both as a function call and as an
environment variable. All variables are prefixed `GOFCHIANTI_` so they are easy
to find in your shell profile:

| Setting | Function | Environment variable |
| --- | --- | --- |
| Cache directory | `gc.set_cache_dir(path)` | `GOFCHIANTI_CACHE` |
| Local dataset dir (used instead of downloading) | `gc.set_dataset_dir(path)` | `GOFCHIANTI_DATASET_DIR` |
| Download base URL | `gc.set_base_url(url)` | `GOFCHIANTI_BASE_URL` |

By default the dataset is downloaded from the IAS SPICE data server:

```
https://spice.osups.universite-paris-saclay.fr/spice-data/contribution_functions/
```

Point `GOFCHIANTI_BASE_URL` (or `gc.set_base_url(...)`) elsewhere to use a
mirror or a local copy. Clear the cache with `gc.clear_cache()`.

## For maintainers

Regenerating and publishing the dataset is a **maintainer-only** task. None of
this is part of the installed package — end users only ever `pip install
gofchianti` and call the API above.

### Prerequisites

- **IDL + SSW/CHIANTI** to compute the raw `G(T, nₑ)` tables. The exact IDL
  routines used to produce the `*_gofnt_v-*.dat` files are vendored, for
  reference, under [`maintainers/idl/`](maintainers/idl/)
  (`compute_gofnt.pro`, `chi_find_transition.pro`).
- A local CHIANTI **abundance** directory, e.g.
  `/usr/local/ssw/packages/chianti/dbase/abundance`.
- **[`rsync`](https://rsync.samba.org/)** with SSH access to the web server
  that hosts the dataset (the default publishing backend).
- Optionally the **[`gh`](https://cli.github.com/) CLI**, authenticated
  (`gh auth login`), if you also publish a GitHub release (secondary backend).
- The dev/maintainer dependencies:

  ```bash
  pip install -e ".[dev]"     # tests
  pip install -e ".[maintainer]"  # ChiantiPy, only if you regenerate inputs
  ```

### 1. Regenerate the dataset

Convert the IDL `.dat` output into the shippable dataset (per-line `.npz`,
`catalog.parquet`, the bundled abundance files and a hashed `manifest.json`):

```bash
python maintainers/convert_dat_to_npz.py \
    --dat-dir ../gofnt \
    --abund-src /usr/local/ssw/packages/chianti/dbase/abundance \
    --out-dir ./dataset
```

The same tool is exposed as a console script after install:

```bash
gofchianti-build-dataset --dat-dir ../gofnt --out-dir ./dataset
```

The build also refreshes the catalogue bundled inside the package
(`src/gofchianti/data/catalog.parquet`) so `available_lines()` works offline.
Pass `--package-data-dir ""` to skip that bundling. The `dataset/` directory
itself is git-ignored — it is distributed as release assets, not committed.

### 2. Publish the dataset

Publishing copies the flat asset set (per-line `.npz`, the `.abund` files,
`catalog.parquet` and `manifest.json`) to wherever end users download it from.
`rsync` over SSH is the **default** backend; a GitHub release via `gh` is a
secondary option. Add `--publish` to publish for real, or `--dry-run` to print
and validate the actions without touching the remote.

#### rsync over SSH (default)

The destination is an `rsync`/`ssh` target `user@host:/path/`, given with
`--dest` or the `GOFCHIANTI_UPLOAD_DEST` environment variable. Files are sent
*flat* so they line up with the flat download URL, and `rsync` only transfers
what changed (safe to re-run):

```bash
# Validate the asset list + command without touching the remote:
python maintainers/convert_dat_to_npz.py --dat-dir ../gofnt --out-dir ./dataset \
    --dry-run --verbose 1

# Publish for real (destination via flag or GOFCHIANTI_UPLOAD_DEST):
python maintainers/convert_dat_to_npz.py --dat-dir ../gofnt --out-dir ./dataset \
    --publish \
    --dest user@host:/var/www/spice-data/contribution_functions/ \
    --verbose 1
```

**Authentication.** For security the tool never stores or forwards a password
itself — authentication is delegated to `ssh`:

- **SSH key** — pass `--ssh-key /path/to/key` (or set `GOFCHIANTI_SSH_KEY`). A
  key already loaded in `ssh-agent` needs no argument and allows fully
  unattended runs.
- **Interactive** — with no key, `ssh` prompts for the password or key
  passphrase directly, in real time; the secret is typed into `ssh` and never
  passes through this tool.

> Password-from-file / password-from-env piping (e.g. `sshpass`) is
> intentionally **not** supported: it would expose the secret in the process
> list and on disk. Use an SSH key or `ssh-agent` for automation.

The public download URL matching the example destination above is
`https://spice.osups.universite-paris-saclay.fr/spice-data/contribution_functions/`.

#### GitHub release (secondary)

```bash
python maintainers/convert_dat_to_npz.py --dat-dir ../gofnt --out-dir ./dataset \
    --publish --target github --repo OWNER/NAME --verbose 1
```

The tag defaults to `dataset-v<dataset_version>` and the upload is idempotent
(re-uploads with `--clobber`). Use `--target both` to publish to the web server
and a GitHub release in one run. End users only fetch from
`releases/latest/download/` if you also point `GOFCHIANTI_BASE_URL` there.

### 3. Run the tests

```bash
pytest
```

The suite runs fully offline against the freshly built `dataset/`. See
[`maintainers/convert_dat_to_npz.py`](maintainers/convert_dat_to_npz.py) for the
converter's full docstring and the IDL output quirks it handles.
