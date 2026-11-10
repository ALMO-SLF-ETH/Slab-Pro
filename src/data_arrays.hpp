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

#include "tools.hpp"
#include <vector>

class Particles {
  public:
    Particles(int Np = 400) {
        x.resize(Np, TV2::Zero()); // position of particles in 2D space
        z.resize(Np, 0.);          // z component for visualization
        v.resize(Np, TV2::Zero()); // velocity of particles in 2D
        vz.resize(Np, 0.);         // vertical velocity component from surface kinematics

        h.resize(Np, 0.);     // particle height
        A.resize(Np, 0.);     // particle area
        mass.resize(Np, 0.);  // particle mass
        theta.resize(Np, 0.); // particle slope angle

        order.resize(Np, 0); // sort
    }

    std::vector<TV2> x;
    std::vector<T> z;
    std::vector<TV2> v;
    std::vector<T> vz;
    std::vector<T> h;
    std::vector<T> A;
    std::vector<T> mass;
    std::vector<T> theta;
    std::vector<int> order;
};

class Grid {
  public:
    Grid(int const& Nx_in = 11, int const& Ny_in = 11) : Nx(Nx_in), Ny(Ny_in) {
        x.resize(Nx, 0.); // x component of position of nodes in 2D space
        y.resize(Ny, 0.); // y component of position of nodes in 2D space

        v.resize(Nx * Ny, TV2::Zero()); // velocity in 2D at nodes
        vz.resize(Nx * Ny, 0.);         // vertical velocity component at nodes
        a.resize(Nx * Ny, TV2::Zero()); // acceleration in 2D

        mass.resize(Nx * Ny, 0.);       // nodal masses
        F.resize(Nx * Ny, TV2::Zero()); // nodal forces

        z.resize(Nx * Ny, 0.);      // elevation
        h.resize(Nx * Ny, 0.);      // flow height
        h_max.resize(Nx * Ny, 0.);  // maximum flow height for each grid point over time (only used for terrain)
        v_max2.resize(Nx * Ny, 0.); // maximum tangent velocity squared for each grid point over time (only used for terrain)

        gradients.resize(Nx * Ny, TV2::Zero()); // gradients in x and y of the terrain on grid
        flip.resize(Nx * Ny, TV2::Zero());      // transfer scheme term
        theta.resize(Nx * Ny, 0.);              // slope angle
    }

    std::vector<T> x;
    std::vector<T> y;
    std::vector<TV2> v;
    std::vector<T> vz;
    std::vector<TV2> a;
    std::vector<T> mass;
    std::vector<TV2> F;
    std::vector<T> z;
    std::vector<T> h;
    std::vector<T> h_max;
    std::vector<T> v_max2;
    std::vector<TV2> gradients;
    std::vector<TV2> flip;
    std::vector<T> theta;

    T xllcorner = 0.;
    T yllcorner = 0.;

    int Nx;
    int Ny;
    T nodata;
};

struct RasterData {
    std::vector<T> x, y, z;
    int ncols, nrows;
    T xllcorner, yllcorner, cellsize, nodata;
};