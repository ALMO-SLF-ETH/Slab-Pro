/*
 * Slab-Pro: Depth-Averaged Material Point Method for Geophysical Mass Flows
 * Copyright (C) 2026 ETH Zurich, SLF Davos
 *
 * Slab-Pro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Slab-Pro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Slab-Pro. If not, see <https://www.gnu.org/licenses/>.
 */

#include "simulation.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::ifstream;
using std::runtime_error;
using std::string;
using std::vector;

template <typename T, typename TV2> RasterData Simulation<T, TV2>::readASC(const std::filesystem::path& path) const {
    ifstream file(path);
    if (!file.is_open()) {
        throw runtime_error("Failed to open ASC file: " + path.string());
    }

    RasterData raster;
    string line;

    // Read header
    file >> line >> raster.ncols;
    file >> line >> raster.nrows;
    file >> line >> raster.xllcorner;
    file >> line >> raster.yllcorner;
    file >> line >> raster.cellsize;
    file >> line >> raster.nodata;

    // Compute x and y
    for (int i = 0; i < raster.ncols; ++i)
        raster.x.push_back(i * raster.cellsize);

    for (int j = 0; j < raster.nrows; ++j)
        raster.y.push_back(j * raster.cellsize);

    // Read elevation values
    raster.z.resize(raster.ncols * raster.nrows);
    for (int j = raster.nrows - 1; j >= 0; --j) { // Reverse row order to match Cartesian (y up ↑)
        for (int i = 0; i < raster.ncols; ++i) {
            file >> raster.z[ind(i, j, raster.nrows)];
        }
    }

    file.close();
    return raster;
}

template <typename T, typename TV2> vector<T> Simulation<T, TV2>::generateGrid(T length, T dx_param) const {
    // Function to generate a grid with a fixed step size dx covering at least the range [0,length]
    vector<T> grid_vec;
    int idx = 0;
    T x = 0.;
    grid_vec.push_back(x);

    // Keep generating points until the next point would be greater than "length"
    while (x < length) {
        idx++;
        x = static_cast<T>(idx) * dx_param;
        grid_vec.push_back(x);
    }

    return grid_vec;
}

template <typename T, typename TV2> void Simulation<T, TV2>::terrainGradients() {
    // Compute gradient at grid points
    T dz_x = 0.;
    T dz_y = 0.;

    for (int i = 0; i < static_cast<int>(terrain.Nx) - 1; ++i) { // gradient at the top and right boundary cannot be computed and will always be zero as there is no
                                                                 // neighboring grid point
        for (int j = 0; j < terrain.Ny - 1; ++j) {
            dz_x = (terrain.z[ind(i + 1, j, terrain.Ny)] + terrain.z[ind(i + 1, j + 1, terrain.Ny)] - terrain.z[ind(i, j, terrain.Ny)] - terrain.z[ind(i, j + 1, terrain.Ny)]) / 2.
                   * one_over_dx;

            dz_y = (terrain.z[ind(i, j + 1, terrain.Ny)] + terrain.z[ind(i + 1, j + 1, terrain.Ny)] - terrain.z[ind(i, j, terrain.Ny)] - terrain.z[ind(i + 1, j, terrain.Ny)]) / 2.
                   * one_over_dx;
            terrain.gradients[ind(i, j, terrain.Ny)] << dz_x, dz_y; // the gradient of each cell is stored in the bottom-left node
        }
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::terrainTheta() {
    // Compute slope angle at grid points in direction of the gradient
    for (int i = 0; i < static_cast<int>(terrain.Nx); ++i) {
        for (int j = 0; j < static_cast<int>(terrain.Ny); ++j) {
            int grid_idx = ind(i, j, terrain.Ny);
            TV2 grad = terrain.gradients[grid_idx];
            terrain.theta[grid_idx] = atan(grad.norm());
        }
    }
}

template <typename T, typename TV2> T Simulation<T, TV2>::particleSlopeAngle(const int p) const {
    // Extract slope angle for a given particle in direction of the gradient from terrain slope angles
    TV2 x_p = particles.x[p];
    int i_base = static_cast<int>(std::floor((x_p(0) - terrain.x[0]) * one_over_dx));
    int j_base = static_cast<int>(std::floor((x_p(1) - terrain.y[0]) * one_over_dx));
    T theta = terrain.theta[ind(i_base, j_base, terrain.Ny)];
    return theta;
}

template <typename T, typename TV2> TV2 Simulation<T, TV2>::particleGradients(const int p) const {
    TV2 x_p = particles.x[p];

    int i_base = static_cast<int>(std::floor((x_p(0) - terrain.x[0]) * one_over_dx)); // Lower-left corner of grid cell that contains particle
    int j_base = static_cast<int>(std::floor((x_p(1) - terrain.y[0]) * one_over_dx));

    TV2 out = terrain.gradients[ind(i_base, j_base, terrain.Ny)];
    return out;
}

template <typename T, typename TV2> T Simulation<T, TV2>::getReleaseHeightFromRaster(const TV2 position) const {
    // Query release height from the raster at a given position.
    // Returns height value or 0.0 if position is out of bounds or nodata

    // assume raster has been loaded when this function is used
    if (release_height_raster.Nx == 0 || release_height_raster.Ny == 0) {
        return 0.0;
    }

    // Position in raster coordinates
    T x_rel = position(0) - release_height_raster.x[0];
    T y_rel = position(1) - release_height_raster.y[0];

    // Check if position is within raster bounds
    if (x_rel < 0 || y_rel < 0 || x_rel >= release_height_raster.x.back() || y_rel >= release_height_raster.y.back()) {
        return 0.0; // Out of bounds
    }

    T dx_local = release_height_raster.x[1] - release_height_raster.x[0];
    T dy_local = release_height_raster.y[1] - release_height_raster.y[0];

    int i = static_cast<int>(floor(x_rel / dx_local));
    int j = static_cast<int>(floor(y_rel / dy_local));

    if (i < 0 || j < 0 || i >= release_height_raster.Nx - 1 || j >= release_height_raster.Ny - 1) {
        return 0.0;
    }

    T fx = (x_rel - i * dx_local) / dx_local;
    T fy = (y_rel - j * dy_local) / dy_local;

    // Get corner values
    T z00 = release_height_raster.z[ind(i, j, release_height_raster.Ny)];
    T z10 = release_height_raster.z[ind(i + 1, j, release_height_raster.Ny)];
    T z01 = release_height_raster.z[ind(i, j + 1, release_height_raster.Ny)];
    T z11 = release_height_raster.z[ind(i + 1, j + 1, release_height_raster.Ny)];

    // Check for nodata values using a simple heuristic: if any corner is nodata, assign height 0.0
    T nodata = release_height_raster.nodata;

    if (abs(z00 - nodata) < 1e-6 || abs(z10 - nodata) < 1e-6 || abs(z01 - nodata) < 1e-6 || abs(z11 - nodata) < 1e-6) {
        return 0.0;
    }

    // Bilinear interpolation
    T z0 = z00 * (1.0 - fx) + z10 * fx;
    T z1 = z01 * (1.0 - fx) + z11 * fx;
    T height = z0 * (1.0 - fy) + z1 * fy;

    // Return the interpolated height, or 0 if it's <= 0
    return height > 0.0 ? height : 0.0;
}

template class Simulation<double, Eigen::Matrix<double, 2, 1>>;