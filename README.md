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

Regenerating the dataset from IDL `.dat` output is a maintainer task and lives in
[`maintainers/convert_dat_to_npz.py`](maintainers/convert_dat_to_npz.py). It is
**not** part of the installed package. See that file's docstring for details.
