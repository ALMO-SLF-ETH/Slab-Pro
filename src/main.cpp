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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <json_file>" << endl;
        return 1;
    }

    fs::path json_path = argv[1];

    Simulation<T, TV2> sim = Simulation<T, TV2>();

    try {
        sim.loadFromJSON(json_path);

        // ---- Initialize simulation ----
        sim.initialize();
        printf("Initialization complete.\n");
        sim.initializeTerrain();
        printf("Terrain initialized.\n");
        sim.initializeParticles();
        printf("Particles initialized.\n");

        // ---- Run simulation ----
        sim.simulate();
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}