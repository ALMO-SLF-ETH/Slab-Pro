# Slab-Pro: Depth-Averaged Material Point Method for Geophysical Mass Flows

[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-blue.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2610.32907-b31b1b.svg)](https://arxiv.org/abs/2610.32907)

## Description
We provide a C++ implementation of the depth-averaged material point method (MPM) for simulating geophysical mass flows such as snow and rock avalanches.

The code supports:

- **Terrain:** a digital elevation model (DEM) in ESRI ASCII grid format, or an analytical surface for verification cases.
- **Release areas:** rectangles, ellipses, polygons from a shapefile, or a raster of spatially varying release heights.
- **Basal friction laws:** Coulomb, Voellmy and μ(I) rheology.
- **Particle sampling:** uniform (regular 2 × 2 per cell) or Poisson-disk sampling.
- **Outputs:** particle .ply data, grid .tif data, and maps of maximum flow height and velocity.

All simulation settings are read from a JSON configuration file.

## Table of contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [User guide](#user-guide)
5. [Examples](#examples)
6. [Contributors](#contributors)
7. [How to cite](#how-to-cite)
8. [License](#license)

---

## Requirements

### Dependencies

- C++ compiler with C++20 support
- CMake ≥ 3.14
- Eigen3 ≥ 3.3
- GDAL ≥ 3.0
- nlohmann_json ≥ 3.10.5

### Computational requirements

- **Operating system:** We tested the software on Ubuntu. On Windows, use the Windows Subsystem for Linux (WSL).
- **CPU/Memory:** The hardware requirements are very low. Any computer that can run a modern operating system should suffice.
- **Disk space:** One particle and one raster file are written per output time. To reduce disk usage, increase `output_dt` or decrease DEM resolution.

---

## Installation

### 1. Install the dependencies

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install build-essential cmake libeigen3-dev libgdal-dev nlohmann-json3-dev
```

### 2. Obtain the source code

```bash
git clone https://github.com/ALMO-SLF-ETH/Slab-Pro.git
cd Slab-Pro
```

### 3. Build

```bash
mkdir build # make build directory
cd build    # change into build directory
cmake ..    # create make files
make        # compile
```

This produces the executable `build/slab-pro`. Compile on the machine where you run the simulations.

---

## Quick start

From the `build` directory:

```bash
./slab-pro ../simConfig.json # run simulation
```

---

## User guide

### Workflow

1. Prepare the terrain (DEM file) and the release area.
2. Write a JSON configuration file (see [Configuration file reference](#configuration-file-reference)).
3. Run `./slab-pro path/to/simConfig.json`.
4. Inspect the results in `<project_dir>/<sim_name>/`.


### Configuration file reference

The configuration is a JSON file. All quantities are in SI units. Every key shown in the tables is required, except that only the sub-object matching the selected `type` is read in the `friction`, `terrain` and `release` sections. Sub-objects for unused types are ignored.

#### General

| Key | Type | Unit | Description |
|---|---|---|---|
| `project_dir` | string | – | Directory under which the simulation folder is created. Relative or absolute paths can be used. |
| `sim_name` | string | – | Name of the simulation. Results are written to `<project_dir>/<sim_name>/`. |

#### `physical_params`

| Key | Type | Unit | Description |
|---|---|---|---|
| `gravity` | number | m s⁻² | Gravitational acceleration. |
| `rho` | number | kg m⁻³ | Bulk density of the flowing material. |

#### `time_params`

| Key | Type | Unit | Description |
|---|---|---|---|
| `output_dt` | number | s | Interval between output snapshots ("dump step"). |
| `total_time` | number | s | Simulated physical time. |

#### `numerical_params`

| Key | Type | Unit | Description |
|---|---|---|---|
| `C_cfl` | number | – | Courant number for the adaptive time step. |
| `dt_max` | number | s | Upper bound for the time step. |
| `flip_ratio` | number | – | Blending between FLIP and PIC velocity updates, from 0 (pure PIC) to 1 (pure FLIP). |
| `sampling_method` | string | – | `"POISSON"` (Poisson-disk sampling) or `"UNIFORM"` (regular 2 × 2 particles per cell). |

#### `friction`

`type` selects the basal friction law: `"COULOMB"`, `"VOELLMY"` or `"MU_I"`.

| `type` | Sub-object | Key | Unit | Description |
|---|---|---|---|---|
| `COULOMB` | `coulomb` | `mu` | – | Friction coefficient μ = tan δ. |
| `VOELLMY` | `voellmy` | `mu` | – | Friction coefficient μ. |
| | | `xi` | m s⁻² | Turbulent friction coefficient ξ. |
| `MU_I` | `mu_I` | `mu_1` | – | Static friction coefficient μ₁. |
| | | `mu_2` | – | Limiting friction coefficient μ₂. |
| | | `I0` | – | Reference inertial number I₀. |
| | | `d` | m | Grain diameter. |
| | | `phi` | – | Solid volume fraction φ. |

#### `terrain`

| `type` | Sub-object | Key | Unit | Description |
|---|---|---|---|---|
| `DEM` | `dem` | `filepath` | – | Path to the DEM in ESRI ASCII grid format (`.asc`). The grid spacing is taken from the file's `cellsize`. |
| `ANALYTICAL` | `analytical` | `Lx` | m | Domain length in x. |
| | | `Ly` | m | Domain length in y. |
| | | `dx` | m | Grid spacing. |

The analytical terrain surface z(x, y) is defined in the source code (function `alt_z` in `Simulation::loadFromJSON` in `src/simulation.cpp`). Changing it requires recompiling.

#### `release`

| Key | Type | Unit | Description |
|---|---|---|---|
| `type` | string | – | `"POLYGON"`, `"RASTER"`, `"RECTANGLE"` or `"ELLIPSE"`. |
| `h0` | number | m | Uniform release height (measured vertically). |

| `type` | Sub-object | Key | Unit | Description |
|---|---|---|---|---|
| `POLYGON` | `polygon` | `filepath` | – | Shapefile (`.shp`) with a single polygon. |
| `RASTER` | `raster` | `filepath` | – | ESRI ASCII grid (`.asc`) of release heights in metres. |
| | | `h_min` | m | Cells with release height below `h_min` are not released. |
| `RECTANGLE` | `rectangle` | `center` | m | `[x, y]` of the rectangle center. |
| | | `Lx`, `Ly` | m | Side lengths. |
| | | `rotation_deg` | ° | Counter-clockwise rotation about the center. |
| `ELLIPSE` | `ellipse` | `center` | m | `[x, y]` of the ellipse center. |
| | | `rx`, `ry` | m | Semi-axes. |
| | | `rotation_deg` | ° | Counter-clockwise rotation about the center. |

### Coordinate systems and units

The DEM must use a projected, metric coordinate reference system with square cells, measured in meters. The polygon and raster release must be in the same CRS as the DEM.

### Outputs

Each run creates the following structure:

```
<project_dir>/
└── <sim_name>/
    ├── <config>.json   # copy of the configuration file used for this run
    ├── max_values.tif  # map of maximum flow height and velocity
    ├── particles/      # particle snapshots, one file per output time
    └── rasters/        # gridded fields, one file per output time
```

| File | Content | Units |
|---|---|---|
| `max_values.tif` | Maximum flow height and velocity over the whole simulation | m, ms⁻¹ |
| `particles/particles_frame_<frame>.ply` | Particle positions, area A, flow height h, velocity vector | m, m², m, ms⁻¹ |
| `rasters/raster_frame_<frame>.tif` | Grid flow height, velocity norm, mass, velocity vector | m, ms⁻¹, kg, ms⁻¹ |

---

## Examples

### Example 1: Brienz Rockslide with raster release

Run the simulation with `./slab-pro ../examples/brienz_raster.json` from the `build` directory.

Note: The resolution in this example is 1 m.

### Example 2: Brienz Rockslide with polygon release

Run the simulation with `./slab-pro ../examples/brienz_polygon.json` from the `build` directory.

---

## Contributors

This code branched off in 2021 from [Matter](https://github.com/larsblatny/matter), a 3D MPM implementation. Early development of the depth-averaged formulation is described in [Guillet et al. (2023)](https://doi.org/10.1029/2023JF007092). The current work is documented in [Kammholz et al. (2026)](https://doi.org/10.48550/arXiv.2610.32907) (see [How to cite](#how-to-cite)).

### Current contributors

| Name | Period | Contribution |
|---|---|---|
| Johan Gaume | 2021–present | Conceptualization, project lead, supervision, funding acquisition |
| Lars Blatny | 2021–present | Supervision, conceptualization, verification, plasticity, developer of [Matter](https://github.com/larsblatny/matter) |
| Johann Kammholz | 2024–present | Main model development |
| Yidong Zhao | 2026–present | Performance optimization, verification |
| Tiziano Di Pietro | 2026–present | JSON configuration input |

### Past contributors

| Name | Period | Contribution |
|---|---|---|
| Denis Steffen | 2021 | Initial implementation on simple topography ([Guillet et al., 2023](https://doi.org/10.1029/2023JF007092)) |
| Louis Guillet | 2021–2023 | Adaptation for snow avalanche release, curvilinear coordinate system, plasticity ([Guillet et al., 2023](https://doi.org/10.1029/2023JF007092)) |
| Thomas Pauze | 2022 | Flow dynamics, including verification and complex terrain, first tests in curvilinear coordinates |

---

## How to cite

If you use this code in your research, please cite:

> Kammholz, J., Blatny, L., Zhao, Y., Pauze, T., & Gaume, J. (2026). A Highly Efficient Depth-Averaged Material Point Method for Geophysical Mass Flows. *arXiv preprint* arXiv:2610.32907. https://doi.org/10.48550/arXiv.2610.32907

```bibtex
@misc{kammholz2026highly,
  title         = {A Highly Efficient Depth-Averaged Material Point Method for Geophysical Mass Flows},
  author        = {Kammholz, Johann and Blatny, Lars and Zhao, Yidong and Pauze, Thomas and Gaume, Johan},
  year          = {2026},
  eprint        = {2610.32907},
  archivePrefix = {arXiv},
  primaryClass  = {cs.CE},
  doi           = {10.48550/arXiv.2610.32907},
  url           = {https://arxiv.org/abs/2610.32907}
}
```

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. See the [LICENSE](LICENSE) for details.

Copyright (C) 2026 ETH Zurich, SLF Davos
