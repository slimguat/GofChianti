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
import gofchianti as gc

# 1. What is available?
df = gc.available_lines()            # pandas DataFrame (ion, λ, version, f, A, ...)
df = gc.available_lines(version="11.0.2")

# 2. Load a line (downloaded + cached on first use, offline afterwards).
cf = gc.get_line("Fe_12", 195.119)

# 3. Evaluate the *bare* contribution function G(T, nₑ).
g = cf.get_gofnt(density=1e9, temperature=1.5e6)

# 4. Abundance-scaled values: multiply by Fe/H from a CHIANTI abundance set.
cf = gc.get_line("Fe_12", 195.119, abundance="sun_photospheric_2021_asplund")
g = cf.get_gofnt(1e9, 1.5e6)         # now × Fe/H
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

## Cache

Downloaded files live in an OS-standard cache directory
(`~/.cache/gofchianti` on Linux). Override it with the `GOFCHIANTI_CACHE`
environment variable or `gc.set_cache_dir(...)`. Clear it with
`gc.clear_cache()`.

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
- The **[`gh`](https://cli.github.com/) CLI**, authenticated
  (`gh auth login`), for publishing releases.
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

### 2. Publish a release

Build and upload every asset to a GitHub release via `gh` (idempotent;
re-uploads with `--clobber`). The tag defaults to `dataset-v<dataset_version>`:

```bash
# Validate the asset list and commands without touching the remote:
python maintainers/convert_dat_to_npz.py --dat-dir ../gofnt --out-dir ./dataset \
    --dry-run-upload --verbose 1

# Publish for real:
python maintainers/convert_dat_to_npz.py --dat-dir ../gofnt --out-dir ./dataset \
    --upload --repo slimguat/GofChianti --verbose 1
```

End users download from the repository's `releases/latest/download/` assets, so
the release that should be served to users must be marked **latest**.

### 3. Run the tests

```bash
pytest
```

The suite runs fully offline against the freshly built `dataset/`. See
[`maintainers/convert_dat_to_npz.py`](maintainers/convert_dat_to_npz.py) for the
converter's full docstring and the IDL output quirks it handles.
