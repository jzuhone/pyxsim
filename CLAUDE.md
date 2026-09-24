# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is pyXSIM

pyXSIM generates synthetic X-ray observations from astrophysical simulation data (grid codes like
FLASH/Enzo/Athena, particle codes like Gadget/AREPO, or hand-built NumPy datasets). It implements the
PHOX algorithm: simulation data -> a `PhotonList` of unabsorbed, un-instrumented photons -> a projected
`EventList` of on-sky events, optionally exported to instrument simulators (e.g. SOXS) for realistic
end products. It depends on `yt` for reading/traversing simulation data and on `soxs` for spectral
generation, response models, and instrument simulation.

## Build & Install

pyXSIM has a compiled Cython core (`pyxsim/lib/*.pyx` -> `sky_functions`, `spectra`, `interpolate`).
After editing any `.pyx` file, or on first setup, rebuild in place:

```
pip install -e .
```

`clean.sh` runs `git clean -fx` to wipe all build artifacts/caches — this is destructive to anything
untracked, use with care.

## Testing

Tests live in `pyxsim/tests/` (not top-level `tests/`, which only holds CI shell scripts). Run:

```
python -m pytest pyxsim/tests
python -m pytest pyxsim/tests/test_beta_model.py::test_beta_model_offaxis   # single test
```

Notes on the test setup:
- Many tests are answer tests that compare against stored reference data (`--answer_store` regenerates
  answers instead of checking them; `--check_dir` points at where spectrum check files live — see
  `pyxsim/tests/conftest.py`). Reference answer tarballs (`pyxsimNN.tar.gz`) and yt/SOXS test data are
  fetched by `tests/ci_install.sh` in CI and configured via `~/.config/yt/yt.toml` and
  `~/.config/soxs/soxs.cfg`; they are not present by default in a local checkout.
  `pyxsim/tests/utils.py` has shared fixtures/data-source builders (e.g. `BetaModelSource`,
  `ParticleBetaModelSource`) used across multiple test files.
- CI (`.github/workflows/`) runs the matrix across Python 3.11-3.14 and both numpy 1.x/2.x on
  Linux/macOS/Windows.

## Linting/Formatting

Ruff is configured in `pyproject.toml` (line length 110) and run via pre-commit
(`.pre-commit-config.yaml`: `ruff-format` + `ruff-check --fix`). Run manually with:

```
ruff format .
ruff check --fix .
```

## Architecture

### Pipeline stages

1. **Source models** (`pyxsim/source_models/`) turn simulation fields into photon-generating recipes.
   All inherit from `SourceModel` (`sources.py`), which defines the `process_data` /
   `_process_chunk` / `_process_data` / `setup_model` contract that `make_photons` drives per data
   chunk. Concrete families:
   - `thermal_sources/` — `ThermalSourceModel` (base.py) with `CIESourceModel` / `NEISourceModel`
     (collisional.py) and `PionSourceModel` / `IGMSourceModel` (photoionization.py, non-equilibrium /
     photoionized plasma via `soxs`'s Cloudy/pion tables).
   - `line_sources.py` — `LineSourceModel` for narrow/broadened emission lines.
   - `power_law_sources.py` — `PowerLawSourceModel` for non-thermal power-law spectra.
   - `xray_binaries.py` is a standalone generator (not a `SourceModel` subclass) that populates mock
     X-ray binary populations from star-formation-rate fields and produces photons directly via
     `make_xrb_particles` / `make_xrb_photons`.
2. **Photon generation** (`photon_list.py`) — `make_photons(...)` walks a yt `data_source` in chunks
   (optionally MPI-parallel via yt's `parallel_objects`/`communication_system`), asks the source model
   to produce photon energies/positions per chunk, and writes an HDF5 photon list
   (`{prefix}.h5` or `{prefix}.{rank}.h5` under MPI). `determine_fields` figures out whether the
   source type is particle- or grid-based and picks the right position/velocity/width fields;
   `find_object_bounds` handles periodic-boundary wrapping for the data region. `PhotonList` wraps a
   photon file for later re-projection.
3. **Projection** (`photon_list.py`: `_project_photons`, `project_photons`, `project_photons_allsky`) —
   projects photons along a line of sight (or as an all-sky map), applying Doppler shifting, foreground
   Galactic absorption (`spectral_models.py: AbsorptionModel`/`TBabsModel`/`WabsModel`), and sky-plane
   scattering (Cython `pyxsim.lib.sky_functions`), producing an `EventList` (`event_list.py`).
4. **Events** (`event_list.py`) — `EventList` holds on-sky event positions/energies and can be written
   out, merged (`utils.merge_files`), and handed off to SOXS for instrument simulation.

### Spectral models (`spectral_models.py`)

`ThermalSpectralModel` and its subclasses (`TableCIEModel`, `Atable1DSpectralModel` ->
`MekalSpectralModel`/`CloudyCIESpectralModel`, `PionSpectralModel`) wrap SOXS's precomputed spectral
tables and interpolate them (`SpectralInterpolator1D`/`2D`) onto the temperature/density/redshift grid
needed by the thermal source models. These are the numerical backends that `ThermalSourceModel`
subclasses call into; they're independent of yt.

### Ongoing refactor: decoupling from yt

Recent work on this branch (see git log: "removing yt-isms out of the core parts", "ripping yt out of
the guts of pyXSIM") is actively separating the yt-specific data-access glue (chunk iteration, field
resolution, parallelism via `yt.utilities.parallel_tools`) from the numerical "kernels" that compute
photon/event quantities from arrays, so the kernels can eventually run against non-yt data sources.
When touching `source_models/` or `photon_list.py`, check whether new code belongs in the yt-facing
layer or in a yt-agnostic kernel, and prefer keeping array-crunching logic free of yt types
(`YTArray`/`YTQuantity`/dataset objects) — `unyt` arrays and plain NumPy are the target.

### Misc

- `internal_absorption.py` computes column-density maps for intrinsic (as opposed to foreground)
  absorption.
- `light_cone.py`'s `XrayLightCone` builds light-cone stacks of photon lists across redshift using yt's
  `LightCone` machinery.
- `utils.py` has cross-cutting helpers: unit parsing (`parse_value`), abundance-table/metallicity field
  handling (`compute_elem_mass_fraction`, `create_metal_fields`), and `merge_files` for combining
  per-rank photon/event files after an MPI run.
- MPI parallelism throughout is via yt's `communication_system`/`parallel_objects`, not `mpi4py`
  directly — the `comm = communication_system.communicators[-1]` pattern recurs across modules.
