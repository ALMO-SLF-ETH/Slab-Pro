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

namespace fs = std::filesystem;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::string;

fs::path createSimulationDirectories(const string& simName, const fs::path& baseDirectory) {
    // Create folder structure for a new simulation run
    fs::path rootDir = baseDirectory / simName;

    try {
        // Create main simulation directory
        if (fs::create_directories(rootDir)) {
            cout << "Directory " << rootDir << " has been created." << endl;
        } else {
            cerr << "Directory " << rootDir << " already exists or could not be created." << endl;
        }

        // Create subdirectories
        fs::create_directories(rootDir / "rasters");
        fs::create_directories(rootDir / "particles");

    } catch (const fs::filesystem_error& e) {
        cerr << "Trying to create directory: " << rootDir << endl;
        cerr << "Filesystem error: " << e.what() << endl;
    } catch (const exception& e) {
        cerr << "General error: " << e.what() << endl;
    }

    return rootDir;
}