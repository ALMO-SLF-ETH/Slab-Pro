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

#include <vector>
#include <cmath>
#include <cstdlib>

class Poisson_disk {
  private:
    int k = 30;               // Max number of attempts to find a new particle
    T r;                      // Minimum distance between particles
    T Lx;                     // Surface width
    T Ly;                     // Surface height
    T a;                      // Cell side length
    int Nx;                   // Number of cells in the x direction
    int Ny;                   // Number of cells in the y direction
    std::vector<TV2> samples; // Resulting particle positions
    std::vector<int> cells;   // Particle cell table (index of particles, -1 if empty)

  public:
    Poisson_disk(T dx, T L1, T L2, int ppc) {
        r = std::sqrt(15.0 / 24.0 * dx * dx / ppc);
        Lx = L1;
        Ly = L2;
        a = r / std::sqrt(2);
        Nx = static_cast<int>(Lx / a) + 1;
        Ny = static_cast<int>(Ly / a) + 1;

        cells.resize(Nx * Ny);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                cells[i * Ny + j] = -1; // Initialize cells to -1 (empty)
            }
        }
    }

    std::vector<TV2> poisson_disk_sampling() {
        TV2 pt = {static_cast<double>(std::rand()) / RAND_MAX * Lx, static_cast<double>(std::rand()) / RAND_MAX * Ly};
        samples.push_back(pt);
        std::vector<int> cell_coords = get_cell_coords(pt);
        cells[cell_coords[0] * Ny + cell_coords[1]] = 0;

        std::vector<int> active = {0};
        int nsamples = 0;

        while (!active.empty()) {
            int idx = std::rand() % static_cast<int>(active.size());
            TV2 refpt = samples[active[idx]];

            pt = get_point(refpt);
            if (pt[0] == -1 && pt[1] == -1) {
                active.erase(active.begin() + idx);
            } else {
                samples.push_back(pt);
                nsamples += 1;
                active.push_back(nsamples);

                cell_coords = get_cell_coords(pt);
                cells[cell_coords[0] * Ny + cell_coords[1]] = nsamples;
            }
        }
        return samples;
    }

  private:
    std::vector<int> get_cell_coords(TV2 pt) {
        std::vector<int> res(2);
        res[0] = static_cast<int>(std::floor(pt[0] / a));
        res[1] = static_cast<int>(std::floor(pt[1] / a));
        return res;
    }

    std::vector<int> get_neighbours(std::vector<int> coords) {
        int dxdy[21][2] = {{-1, -2}, {0, -2}, {1, -2}, {-2, -1}, {-1, -1}, {0, -1}, {1, -1}, {2, -1}, {-2, 0}, {-1, 0}, {1, 0},
                           {2, 0},   {-2, 1}, {-1, 1}, {0, 1},   {1, 1},   {2, 1},  {-1, 2}, {0, 2},  {1, 2},  {0, 0}};
        std::vector<int> neighbours;

        for (int i = 0; i < 21; i++) {
            std::vector<int> neighbour_coords(2);
            neighbour_coords[0] = coords[0] + dxdy[i][0];
            neighbour_coords[1] = coords[1] + dxdy[i][1];
            if (neighbour_coords[0] >= 0 && neighbour_coords[0] < Nx && neighbour_coords[1] >= 0 && neighbour_coords[1] < Ny) {
                int neighbour_cell = cells[neighbour_coords[0] * Ny + neighbour_coords[1]];
                if (neighbour_cell != -1) {
                    neighbours.push_back(neighbour_cell);
                }
            }
        }
        return neighbours;
    }

    bool point_valid(TV2 pt) {
        std::vector<int> cell_coords = get_cell_coords(pt);
        std::vector<int> neighbours = get_neighbours(cell_coords);

        for (int i = 0; i < static_cast<int>(neighbours.size()); i++) {
            TV2 nearby_pt = samples[neighbours[i]];
            T distance2 = (nearby_pt[0] - pt[0]) * (nearby_pt[0] - pt[0]) + (nearby_pt[1] - pt[1]) * (nearby_pt[1] - pt[1]);
            if (distance2 < r * r) {
                return false;
            }
        }
        return true;
    }

    TV2 get_point(TV2 refpt) {
        for (int i = 0; i < k; i++) {
            T rho = static_cast<double>(std::rand()) / RAND_MAX * r + r;
            T theta = static_cast<double>(std::rand()) / RAND_MAX * 2 * M_PI;
            TV2 pt = {refpt[0] + rho * std::cos(theta), refpt[1] + rho * std::sin(theta)};
            if (pt[0] >= 0 && pt[0] <= Lx && pt[1] >= 0 && pt[1] <= Ly) {
                if (point_valid(pt)) {
                    return pt;
                }
            }
        }
        return {-1, -1}; // If no point is found
    }
};
