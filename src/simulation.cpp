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

#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <iomanip>
#include <vector>
#include <cmath>
#include <filesystem>
#include <limits>
#include <chrono>
#include <iostream>

#include "polygon_tools.hpp"
#include "poisson_disk_sampling.hpp"

using std::ceil;
using std::copy;
using std::cout;
using std::endl;
using std::fill;
using std::left;
using std::max;
using std::max_element;
using std::min;
using std::min_element;
using std::minmax_element;
using std::numeric_limits;
using std::right;
using std::setw;
using std::string;
using std::vector;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;

template <typename T, typename TV2> void Simulation<T, TV2>::initialize() {
    // Initialize simulation.
    sim_dir_path = createSimulationDirectories(sim_name, project_dir_path);
    particle_path = sim_dir_path / "particles";
    raster_path = sim_dir_path / "rasters";
    nFiles = static_cast<int>(ceil(total_time / output_dt));
    padding = static_cast<int>(floor(log10(nFiles))) + 1;

    // Copy JSON file to sim folder
    try {
        // Determine filename (keep original name)
        std::filesystem::path json_filename = std::filesystem::path(json_file_path).filename();
        std::filesystem::path dest = sim_dir_path / json_filename;

        std::filesystem::copy_file(json_file_path, dest, std::filesystem::copy_options::overwrite_existing);
        std::cout << "Copied JSON to " << dest << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to copy JSON file: " << e.what() << "\n";
    }

    time = 0.0;
    frame = 0;

    runtime_updateDt = 0;
    runtime_sort = 0.;
    runtime_remesh = 0;
    runtime_P2G = 0;
    runtime_computeGridVelocity = 0;
    runtime_G2P = 0;
}

// Helper to convert friction type string to enum
template <typename T, typename TV2> FrictionType Simulation<T, TV2>::parseFrictionType(const std::string& str) {
    if (str == "COULOMB")
        return FrictionType::Coulomb;
    if (str == "VOELLMY")
        return FrictionType::Voellmy;
    if (str == "MU_I")
        return FrictionType::mu_I;
    throw std::runtime_error("Unknown friction type: " + str);
}

// Helper to convert terrain type string to enum
template <typename T, typename TV2> TerrainType Simulation<T, TV2>::parseTerrainType(const std::string& str) {
    if (str == "DEM")
        return TerrainType::DEM;
    if (str == "ANALYTICAL")
        return TerrainType::Analytical;
    throw std::runtime_error("Unknown terrain type: " + str);
}

// Helper to convert release shape string to enum
template <typename T, typename TV2> ReleaseShapeType Simulation<T, TV2>::parseReleaseShapeType(const std::string& str) {
    if (str == "RECTANGLE")
        return ReleaseShapeType::Rectangle;
    if (str == "ELLIPSE")
        return ReleaseShapeType::Ellipse;
    if (str == "POLYGON")
        return ReleaseShapeType::Polygon;
    if (str == "RASTER")
        return ReleaseShapeType::Raster;
    throw std::runtime_error("Unknown release shape type: " + str);
}

template <typename T, typename TV2> void Simulation<T, TV2>::loadFromJSON(const std::filesystem::path& json_path) {
    json_file_path = json_path;
    if (!std::filesystem::exists(json_path))
        throw std::runtime_error("JSON file not found: " + json_path.string());

    nlohmann::json j;
    std::ifstream json_file(json_path);
    json_file >> j;

    // ---- General parameters ----
    std::filesystem::path project_dir = j.at("project_dir").get<std::string>();
    project_dir_path = project_dir.string();
    sim_name = j.at("sim_name").get<std::string>();

    // ---- Physical parameters ----
    const auto& phys = j.at("physical_params");
    gravity = phys.at("gravity").get<double>();
    rho = phys.at("rho").get<double>();

    // ---- Time parameters ----
    const auto& time_params = j.at("time_params");
    output_dt = time_params.at("output_dt").get<double>();
    total_time = time_params.at("total_time").get<double>();

    // ---- Numerical parameters ----
    const auto& num_params = j.at("numerical_params");
    C_cfl = num_params.at("C_cfl").get<double>();
    dt_max = num_params.at("dt_max").get<double>();
    flip_ratio = num_params.at("flip_ratio").get<double>();
    std::string method = num_params.at("sampling_method").get<std::string>();
    if (method == "UNIFORM")
        sampling_method = SamplingMethod::Uniform;
    else if (method == "POISSON")
        sampling_method = SamplingMethod::Poisson;
    else
        throw std::runtime_error("Unknown sampling method: " + method);

    // ---- Friction ----
    const auto& fr = j.at("friction");
    friction.type = parseFrictionType(fr.at("type").get<std::string>());
    if (friction.type == FrictionType::Coulomb) {
        friction.coulomb.mu = fr.at("coulomb").at("mu").get<double>();
    } else if (friction.type == FrictionType::Voellmy) {
        friction.voellmy.mu = fr.at("voellmy").at("mu").get<double>();
        friction.voellmy.xi = fr.at("voellmy").at("xi").get<double>();
    } else if (friction.type == FrictionType::mu_I) {
        friction.muI.mu_1 = fr.at("mu_I").at("mu_1").get<double>();
        friction.muI.mu_2 = fr.at("mu_I").at("mu_2").get<double>();
        friction.muI.I0 = fr.at("mu_I").at("I0").get<double>();
        friction.muI.d = fr.at("mu_I").at("d").get<double>();
        friction.muI.phi = fr.at("mu_I").at("phi").get<double>();
    }

    // ---- Terrain ----
    const auto& terrain_params = j.at("terrain");
    terrain_type = parseTerrainType(terrain_params.at("type").get<std::string>());
    if (terrain_type == TerrainType::DEM) {
        dem_filepath = terrain_params.at("dem").at("filepath").get<std::string>();
    } else if (terrain_type == TerrainType::Analytical) {
        terrain_Lx = terrain_params.at("analytical").at("Lx").get<double>();
        terrain_Ly = terrain_params.at("analytical").at("Ly").get<double>();
        dx = terrain_params.at("analytical").at("dx").get<double>();

        // Hardcoded analytical terrain function (only used for verification cases)
        alt_z = [](const T& x, const T& y) -> T {
            double theta_rad = 5.0 * M_PI / 180.0;
            return -x * std::tan(theta_rad) - y * std::tan(theta_rad);
        };
    }

    // ---- Release ----
    const auto& release = j.at("release");
    release_shape_type = parseReleaseShapeType(release.at("type").get<std::string>());
    h0 = release.at("h0").get<double>();
    if (release_shape_type == ReleaseShapeType::Rectangle) {
        rect.center = TV2(release.at("rectangle").at("center")[0].get<double>(), release.at("rectangle").at("center")[1].get<double>());
        rect.Lx = release.at("rectangle").at("Lx").get<double>();
        rect.Ly = release.at("rectangle").at("Ly").get<double>();
        rect.rotation = release.at("rectangle").at("rotation_deg").get<double>() * M_PI / 180.0;
    } else if (release_shape_type == ReleaseShapeType::Ellipse) {
        ell.center = TV2(release.at("ellipse").at("center")[0].get<double>(), release.at("ellipse").at("center")[1].get<double>());
        ell.rx = release.at("ellipse").at("rx").get<double>();
        ell.ry = release.at("ellipse").at("ry").get<double>();
        ell.rotation = release.at("ellipse").at("rotation_deg").get<double>() * M_PI / 180.0;
    } else if (release_shape_type == ReleaseShapeType::Polygon) {
        poly.filepath = release.at("polygon").at("filepath").get<std::string>();
    } else if (release_shape_type == ReleaseShapeType::Raster) {
        rast.filepath = release.at("raster").at("filepath").get<std::string>();
        rast.h_min = release.at("raster").at("h_min").get<double>();
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::initializeTerrain() {
    // Create grid called terrain which contains the whole domain.
    switch (terrain_type) {
        case TerrainType::DEM: {
            RasterData raster = readASC(dem_filepath); // Load terrain from ASC file

            terrain = Grid(raster.ncols, raster.nrows); // Set up grid
            terrain.x = raster.x;
            terrain.y = raster.y;
            terrain.z = raster.z;
            terrain.xllcorner = raster.xllcorner;
            terrain.yllcorner = raster.yllcorner;

            terrain.Nx = static_cast<int>(terrain.x.size());
            terrain.Ny = static_cast<int>(terrain.y.size());

            dx = raster.cellsize; // Update dx based on ASC file
            one_over_dx = 1. / dx;
            break;
        }
        case TerrainType::Analytical: {
            one_over_dx = 1. / dx; // dx is defined in input script

            auto terrain_x = generateGrid(terrain_Lx, dx); // Generate terrain x and y coordinates using step size dx
            auto terrain_y = generateGrid(terrain_Ly, dx);

            terrain.Nx = static_cast<int>(terrain_x.size());
            terrain.Ny = static_cast<int>(terrain_y.size());

            terrain = Grid(terrain.Nx, terrain.Ny);
            terrain.x = terrain_x; // xllcorner, yllcorner stays initialized to zero for analytical terrain
            terrain.y = terrain_y;

            for (int i = 0; i < terrain.Nx; ++i) {
                for (int j = 0; j < terrain.Ny; ++j) {
                    terrain.z[ind(i, j, terrain.Ny)] = alt_z(terrain.x[i], terrain.y[j]); // terrain height
                }
            }
            break;
        }
    }

    terrainGradients(); // gradient at each grid point
    terrainTheta();     // terrain angle in direction of the gradient

    fill(terrain.h_max.begin(), terrain.h_max.end(), 0.);   // Initialize h_max over the full terrain for max flow height tracking
    fill(terrain.v_max2.begin(), terrain.v_max2.end(), 0.); // Initialize max tangent velocity squared
}

template <typename T, typename TV2> void Simulation<T, TV2>::initializeParticles() {
    // Initialize particles with position, velocity, density, etc.
    vector<TV2> samples, particle_positions;
    T x_min = 0, x_max = 0, y_min = 0, y_max = 0; // Bounding box
    T release_volume = 0;

    switch (release_shape_type) { // Determine bounding box
        case ReleaseShapeType::Rectangle: {
            T diag = sqrt(rect.Lx * rect.Lx + rect.Ly * rect.Ly) / 2.;
            x_min = rect.center(0) - diag;
            x_max = rect.center(0) + diag;
            y_min = rect.center(1) - diag;
            y_max = rect.center(1) + diag;
            release_area = rect.Lx * rect.Ly;
            break;
        }

        case ReleaseShapeType::Ellipse: {
            T biggest_radius = max(ell.rx, ell.ry);
            x_min = ell.center(0) - biggest_radius;
            x_max = ell.center(0) + biggest_radius;
            y_min = ell.center(1) - biggest_radius;
            y_max = ell.center(1) + biggest_radius;
            release_area = M_PI * ell.rx * ell.ry;
            break;
        }

        case ReleaseShapeType::Polygon: {
            GDALAllRegister();

            x_min = y_min = numeric_limits<T>::max();
            x_max = y_max = numeric_limits<T>::lowest();

            poly.vertices = loadPolygonFromShapefile(poly.filepath);

            for (auto& pt : poly.vertices) { // Shift polygon into local coordinates (starting at 0,0)
                pt(0) -= terrain.xllcorner;
                pt(1) -= terrain.yllcorner;
            }

            for (const auto& pt : poly.vertices) {
                x_min = std::min(x_min, static_cast<T>(pt(0)));
                x_max = std::max(x_max, static_cast<T>(pt(0)));
                y_min = std::min(y_min, static_cast<T>(pt(1)));
                y_max = std::max(y_max, static_cast<T>(pt(1)));
            }
            release_area = computePolygonArea(poly.vertices);
            break;
        }

        case ReleaseShapeType::Raster: {
            // load raster now and compute bounds/area
            RasterData raster = readASC(rast.filepath);
            release_height_raster = Grid(raster.ncols, raster.nrows);
            release_height_raster.x.resize(raster.ncols);
            release_height_raster.y.resize(raster.nrows);

            for (int i = 0; i < raster.ncols; ++i)
                release_height_raster.x[i] = raster.xllcorner - terrain.xllcorner + i * raster.cellsize;

            for (int j = 0; j < raster.nrows; ++j)
                release_height_raster.y[j] = raster.yllcorner - terrain.yllcorner + j * raster.cellsize;
            release_height_raster.z = raster.z;
            release_height_raster.Nx = raster.ncols;
            release_height_raster.Ny = raster.nrows;
            release_height_raster.xllcorner = raster.xllcorner;
            release_height_raster.yllcorner = raster.yllcorner;
            release_height_raster.nodata = raster.nodata;

            x_min = raster.xllcorner - terrain.xllcorner;
            y_min = raster.yllcorner - terrain.yllcorner;

            x_max = x_min + raster.cellsize * (raster.ncols - 1);
            y_max = y_min + raster.cellsize * (raster.nrows - 1);

            // compute release area and release volume from raster
            int count = 0;
            T cell_area = raster.cellsize * raster.cellsize;
            for (int j = 0; j < raster.nrows; ++j) {
                for (int i = 0; i < raster.ncols; ++i) {
                    const T h_cell = raster.z[ind(i, j, raster.nrows)];
                    if (h_cell >= rast.h_min && h_cell > 0.) {
                        ++count;
                        release_volume += h_cell * cell_area;
                    }
                }
            }
            release_area = count * cell_area;
            break;
        }
    }

    if (release_shape_type != ReleaseShapeType::Raster) {
        release_volume = release_area * h0; // uniform release height
    }

    x_min = floor(x_min * one_over_dx) * dx;
    x_max = ceil(x_max * one_over_dx) * dx;
    y_min = floor(y_min * one_over_dx) * dx;
    y_max = ceil(y_max * one_over_dx) * dx;

    switch (sampling_method) {
        case SamplingMethod::Uniform: {
            int Nx_cells = static_cast<int>(round((x_max - x_min) * one_over_dx)); // number of grid cells occupied by particles
            int Ny_cells = static_cast<int>(round((y_max - y_min) * one_over_dx));

            vector<T> disp_x(ppc); // displacement to distribute particles within each cell
            vector<T> disp_y(ppc);

            disp_x[0] = 0.25;
            disp_x[1] = 0.75;
            disp_x[2] = 0.25;
            disp_x[3] = 0.75;
            disp_y[0] = 0.25;
            disp_y[1] = 0.25;
            disp_y[2] = 0.75;
            disp_y[3] = 0.75;

            for (int i = 0; i < Nx_cells; ++i) {
                for (int j = 0; j < Ny_cells; ++j) {
                    for (int d = 0; d < ppc; ++d) {
                        T px = ((double)i + disp_x[d]) * dx;
                        T py = ((double)j + disp_y[d]) * dx;
                        samples.push_back({px, py});
                    }
                }
            }
            break;
        }

        case SamplingMethod::Poisson: {
            Poisson_disk sampler(dx, x_max - x_min, y_max - y_min, ppc);
            samples = sampler.poisson_disk_sampling();
            break;
        }
    }

    // filter samples according to release shape
    for (const TV2& p : samples) {
        TV2 shifted = p + TV2(x_min, y_min);
        bool keep = false;
        switch (release_shape_type) {
            case ReleaseShapeType::Rectangle: {
                TV2 local = rotatePoint(shifted - rect.center, -rect.rotation);
                keep = abs(local(0)) <= rect.Lx / 2. && abs(local(1)) <= rect.Ly / 2.;
                break;
            }
            case ReleaseShapeType::Ellipse: {
                TV2 local = rotatePoint(shifted - ell.center, -ell.rotation);
                keep = pow(local(0) / ell.rx, 2.) + pow(local(1) / ell.ry, 2.) <= 1.;
                break;
            }
            case ReleaseShapeType::Polygon: {
                keep = isPointInPolygon(shifted(0), shifted(1), poly.vertices);
                break;
            }
            case ReleaseShapeType::Raster: {
                // raster defines both location and height, keep if height >= h_min and > 0
                const T h_rel = getReleaseHeightFromRaster(shifted);
                keep = (h_rel >= rast.h_min && h_rel > 0.);
                break;
            }
        }
        if (keep) {
            particle_positions.push_back(shifted);
        }
    }

    Np = static_cast<int>(particle_positions.size());
    particles = Particles(Np);

    // Initial particle height = release height at the particle position
    vector<T> h_init(Np, h0);
    T sum_h = 0.0;
    for (int p = 0; p < Np; ++p) {
        if (release_shape_type == ReleaseShapeType::Raster)
            h_init[p] = getReleaseHeightFromRaster(particle_positions[p]);
        sum_h += h_init[p];
    }

    // Initial particle area: both samplers give a spatially uniform particle density, so every particle represents the same planform area.
    // Scale it such that sum_p h_p A_p = V_0 exactly. For a uniform release height this is A_p = release_area / Np.
    const T A_init = release_volume / sum_h;

    for (int p = 0; p < Np; ++p) {
        particles.x[p] = particle_positions[p];

        particles.order[p] = p;

        particles.h[p] = h_init[p];

        particles.A[p] = A_init;

        particles.mass[p] = rho * particles.h[p] * particles.A[p];
        particles.theta[p] = particleSlopeAngle(p);

        const TV2& x_p = particles.x[p];

        int i_base = static_cast<int>((x_p(0) - terrain.x[0]) * one_over_dx) - 1;
        int j_base = static_cast<int>((x_p(1) - terrain.y[0]) * one_over_dx) - 1;

        for (int i = i_base; i < i_base + 4; ++i) {
            for (int j = j_base; j < j_base + 4; ++j) {
                T weight = wip(x_p(0), x_p(1), terrain.x[i], terrain.y[j], one_over_dx);
                particles.z[p] += terrain.z[ind(i, j, terrain.Ny)] * weight;
            }
        }
    }

    // print total volume
    T total_volume = 0.0;
    for (int p = 0; p < Np; ++p) {
        total_volume += particles.h[p] * particles.A[p];
    }
    cout << "Total volume: " << total_volume << endl;

    cout << "Number of particles: " << Np << endl;
}

template <typename T, typename TV2> TV2 Simulation<T, TV2>::rotatePoint(const TV2& p, T angle) const {
    T c = cos(angle);
    T s = sin(angle);
    return TV2(c * p(0) - s * p(1), s * p(0) + c * p(1));
}

template <typename T, typename TV2> void Simulation<T, TV2>::simulate() {
    // Manage overall execution of simulation.
    // saveTerrain();
    cout << "Starting simulation..." << endl;
    int time_step = 0;
    frame = 0;
    auto start_time = steady_clock::now(); // Record start time
    int width = 8;                         // Number of digits to align cout

    remesh();
    P2G();
    computeGridVelocity(); // dt is initialized to zero, so this is safe; it is needed for normalizing grid.h for time step 0 output
    saveParticleData();
    saveGridData();
    cout << "Time step: " << setw(width) << time_step << "    Time: " << time << " / " << total_time << " s" << endl;

    while (time < total_time) {
        time_step++;
        advanceStep(time_step);
        time += dt;
        if (output_dt * (frame + 1) - time < 1e-13) {
            frame++;
            saveParticleData();
            saveGridData();
            cout << "Time step: " << setw(width) << time_step << "    Time: " << left << setw(10) << time << right << " / " << total_time << " s" << endl;
        }
    }

    auto end_time = steady_clock::now(); // Record end time
    cout << "Simulation took " << duration_cast<milliseconds>(end_time - start_time).count() << " ms." << endl;
    cout << "updateDt took " << runtime_updateDt << " s." << endl;
    cout << "sort took " << runtime_sort << " s." << endl;
    cout << "remesh took " << runtime_remesh << " s." << endl;
    cout << "P2G took " << runtime_P2G << " s." << endl;
    cout << "computeGridVelocity took " << runtime_computeGridVelocity << " s." << endl;
    cout << "G2P took " << runtime_G2P << " s." << endl;

    saveMaxData();
}

template <typename T, typename TV2> void Simulation<T, TV2>::advanceStep(int time_step) {
    // Advance simulation by one time step.
    timer t_updateDt;
    t_updateDt.start();
    updateDt();
    t_updateDt.stop();
    runtime_updateDt += t_updateDt.get_timing();

    timer t_remesh;
    t_remesh.start();
    remesh();
    t_remesh.stop();
    runtime_remesh += t_remesh.get_timing();

    timer t_sort;
    t_sort.start();
    if (time_step % 50 == 0) {
        sort_particles_by_grid_block();
    }
    t_sort.stop();
    runtime_sort += t_sort.get_timing();

    timer t_P2G;
    t_P2G.start();
    P2G();
    t_P2G.stop();
    runtime_P2G += t_P2G.get_timing();

    timer t_computeGridVelocity;
    t_computeGridVelocity.start();
    computeGridVelocity();
    t_computeGridVelocity.stop();
    runtime_computeGridVelocity += t_computeGridVelocity.get_timing();

    timer t_G2P;
    t_G2P.start();
    G2P();
    t_G2P.stop();
    runtime_G2P += t_G2P.get_timing();
}

template <typename T, typename TV2> void Simulation<T, TV2>::updateDt() {
    // Update time step based on current simulation conditions.
    auto max_velocity_it = max_element(particles.v.begin(), particles.v.end(), [](const TV2& v1, const TV2& v2) {
        return v1.squaredNorm() < v2.squaredNorm();
    });                                    // Find particle with maximum velocity
    T max_speed = max_velocity_it->norm(); // Get maximum velocity from iterator

    if (max_speed > 0.0) {
        dt = min(C_cfl * dx / max_speed, dt_max); // Minimum of CFL condition and maximum time step
    } else {
        dt = dt_max; // All particles at rest
    }

    T next_frame_time = output_dt * (frame + 1);
    if (time + dt > next_frame_time - 1e-12) { // If next time step almost reaches or overshoots next frame time, snap to it
        dt = next_frame_time - time;
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::remesh() {
    // Remesh grid based on particle positions.
    auto x_cmp = [](const TV2& a, const TV2& b) { // Comparison function
        return a(0) < b(0);
    };
    auto y_cmp = [](const TV2& a, const TV2& b) {
        return a(1) < b(1);
    };

    // Iterator pointing to min/max of particle positions
    auto [min_x_it, max_x_it] = minmax_element(particles.x.begin(), particles.x.end(), x_cmp);
    auto [min_y_it, max_y_it] = minmax_element(particles.x.begin(), particles.x.end(), y_cmp);

    T min_x = (*min_x_it)(0); // min/max particle positions in x and y
    T max_x = (*max_x_it)(0);
    T min_y = (*min_y_it)(1);
    T max_y = (*max_y_it)(1);

    const T epsilon = 1e-8 * dx;
    int i_base_left = static_cast<int>((min_x - terrain.x[0]) * one_over_dx - 1.0 - epsilon);  // same as in P2G/G2P with additional epsilon for floating point safety
    int i_base_right = static_cast<int>((max_x - terrain.x[0]) * one_over_dx + 2.0 + epsilon); // +1.5 from lower limit because quadratic is 3x3
    int j_base_bottom = static_cast<int>((min_y - terrain.y[0]) * one_over_dx - 1.0 - epsilon);
    int j_base_top = static_cast<int>((max_y - terrain.y[0]) * one_over_dx + 2.0 + epsilon);

    grid_i_offset = i_base_left; // Store offset of grid relative to terrain for saving temporal maximum values
    grid_j_offset = j_base_bottom;

    grid.Nx = i_base_right - i_base_left + 1; // Grid dimensions in x and y
    grid.Ny = j_base_top - j_base_bottom + 1;
    int totalSize = grid.Nx * grid.Ny;

    grid.x.resize(grid.Nx); // Resize grid to new grid size
    grid.y.resize(grid.Ny);
    grid.z.resize(totalSize);
    grid.h.resize(totalSize);
    fill(grid.h.begin(), grid.h.end(), 0.);
    grid.mass.resize(totalSize);
    fill(grid.mass.begin(), grid.mass.end(), 0.);
    grid.F.resize(totalSize);
    fill(grid.F.begin(), grid.F.end(), TV2::Zero());
    grid.flip.resize(totalSize);
    fill(grid.flip.begin(), grid.flip.end(), TV2::Zero());
    grid.a.resize(totalSize);
    fill(grid.a.begin(), grid.a.end(), TV2::Zero());
    grid.v.resize(totalSize);
    fill(grid.v.begin(), grid.v.end(), TV2::Zero());
    grid.vz.resize(totalSize);
    fill(grid.vz.begin(), grid.vz.end(), 0.);

    // Populate grid x and y values from terrain
    copy(terrain.x.begin() + i_base_left, terrain.x.begin() + i_base_right + 1, grid.x.begin());
    copy(terrain.y.begin() + j_base_bottom, terrain.y.begin() + j_base_top + 1, grid.y.begin());

    for (int i = i_base_left; i <= i_base_right; ++i) { // Populate grid z values from terrain
        for (int j = j_base_bottom; j <= j_base_top; ++j) {
            grid.z[ind(i - i_base_left, j - j_base_bottom, grid.Ny)] = terrain.z[ind(i, j, terrain.Ny)];
        }
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::sort_particles_by_grid_block() {
    constexpr int SORT_BLOCK_SIZE = 4;

    const int nbx = (grid.Nx + SORT_BLOCK_SIZE - 1) / SORT_BLOCK_SIZE;
    const int nby = (grid.Ny + SORT_BLOCK_SIZE - 1) / SORT_BLOCK_SIZE;
    const int nblocks = nbx * nby;

    std::vector<int> counts(nblocks, 0);
    std::vector<int> block_id(Np, 0);

    // Particles loop
    for (int p = 0; p < Np; ++p) {
        const TV2& x_p = particles.x[p];

        int i_base = static_cast<int>((x_p(0) - grid.x[0]) * one_over_dx - 0.5);
        int j_base = static_cast<int>((x_p(1) - grid.y[0]) * one_over_dx - 0.5);

        // Clamp for safety
        if (i_base < 0)
            i_base = 0;
        if (j_base < 0)
            j_base = 0;
        if (i_base >= grid.Nx)
            i_base = grid.Nx - 1;
        if (j_base >= grid.Ny)
            j_base = grid.Ny - 1;

        const int bi = i_base / SORT_BLOCK_SIZE;
        const int bj = j_base / SORT_BLOCK_SIZE;
        const int b = bi * nby + bj;

        block_id[p] = b;
        counts[b]++;
    }

    // Prefix sum
    std::vector<int> offsets(nblocks, 0);
    {
        int sum = 0;
        for (int b = 0; b < nblocks; ++b) {
            offsets[b] = sum;
            sum += counts[b];
        }
    }

    // Permutation
    std::vector<int> next = offsets;
    for (int p = 0; p < Np; ++p) {
        const int b = block_id[p];
        particles.order[next[b]++] = p;
    }

    apply_particle_permutation(particles.order);
}

template <typename T, typename TV2> void Simulation<T, TV2>::apply_particle_permutation(const std::vector<int>& order) {
    auto permute_vec2 = [&](std::vector<TV2>& a) {
        std::vector<TV2> tmp(a.size());
        for (int new_idx = 0; new_idx < Np; ++new_idx) {
            tmp[new_idx] = a[order[new_idx]];
        }
        a.swap(tmp);
    };

    auto permute_scalar = [&](std::vector<T>& a) {
        std::vector<T> tmp(a.size());
        for (int new_idx = 0; new_idx < Np; ++new_idx) {
            tmp[new_idx] = a[order[new_idx]];
        }
        a.swap(tmp);
    };

    permute_vec2(particles.x);
    permute_scalar(particles.z);
    permute_vec2(particles.v);
    permute_scalar(particles.vz);
    permute_scalar(particles.h);
    permute_scalar(particles.A);
    permute_scalar(particles.mass);
    permute_scalar(particles.theta);
}

template <typename T, typename TV2> void Simulation<T, TV2>::P2G() {
    // Interpolate particle quantities to grid.

    for (int p = 0; p < Np; ++p) {
        const TV2 x_p = particles.x[p]; // Extract particle properties
        const T height_p = particles.h[p];
        const T A_p = particles.A[p];
        const TV2 v_p = particles.v[p];
        const T m_p = particles.mass[p];

        const TV2 vm_p = v_p * m_p;
        const T vzm_p = particles.vz[p] * m_p;
        const T hm_p = height_p * m_p;
        const T h_A = height_p * A_p;

        int i_base = static_cast<int>((x_p(0) - grid.x[0]) * one_over_dx - 0.5); // Index of particle position in grid array
        int j_base = static_cast<int>((x_p(1) - grid.y[0]) * one_over_dx - 0.5);

        const TV2 part_grad = particleGradients(p); // Slope angle at particles
        particles.theta[p] = particleSlopeAngle(p);
        const T cos_part_theta = cos(particles.theta[p]);

        const T coeff = -m_p * gravity * cos_part_theta * cos_part_theta;
        const T m_b_px = coeff * part_grad(0); // coeff has minus sign because the force goes in opposite direction of the slope
        const T m_b_py = coeff * part_grad(1);

        TM2 sigma;
        sigma(0, 0) = -0.5 * rho * gravity * cos_part_theta * cos_part_theta * height_p; // compressive stress is negative
        sigma(1, 1) = sigma(0, 0);
        sigma(0, 1) = 0.0;
        sigma(1, 0) = 0.0;

        // Precompute weights
        T wx[3], wy[3], dwx[3], dwy[3];
        for (int di = 0; di < 3; ++di) {
            int i = i_base + di;
            T dx_ = (x_p(0) - grid.x[i]) * one_over_dx;
            wx[di] = N(dx_);
            dwx[di] = dNdu(dx_) * one_over_dx;
        }
        for (int dj = 0; dj < 3; ++dj) {
            int j = j_base + dj;
            T dy_ = (x_p(1) - grid.y[j]) * one_over_dx;
            wy[dj] = N(dy_);
            dwy[dj] = dNdu(dy_) * one_over_dx;
        }

        for (int di = 0; di < 3; ++di) {
            int i = i_base + di;
            for (int dj = 0; dj < 3; ++dj) {
                int j = j_base + dj;

                const T weight = wx[di] * wy[dj];

                if (weight > 3e-3) {
                    const T grad_weight_x = dwx[di] * wy[dj];
                    const T grad_weight_y = dwy[dj] * wx[di];

                    int grid_idx = ind(i, j, grid.Ny);

                    grid.mass[grid_idx] += m_p * weight;
                    grid.v[grid_idx] += vm_p * weight;   // momentum
                    grid.vz[grid_idx] += vzm_p * weight; // vertical momentum
                    grid.h[grid_idx] += hm_p * weight;   // mass-weighted height, for export of grid data and friction laws only

                    // Internal force has minus sign from integration by parts in weak form
                    grid.F[grid_idx](0) += m_b_px * weight - h_A * (sigma(0, 0) * grad_weight_x + sigma(0, 1) * grad_weight_y);
                    grid.F[grid_idx](1) += m_b_py * weight - h_A * (sigma(1, 1) * grad_weight_y + sigma(1, 0) * grad_weight_x);
                }
            }
        }
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::computeGridVelocity() {
    // Compute acceleration and velocity on grid.
    for (int i = 0; i < grid.Nx; ++i) {
        for (int j = 0; j < grid.Ny; ++j) {
            int grid_idx = ind(i, j, grid.Ny);
            const T m_i = grid.mass[grid_idx];

            if (m_i > 1e-6) {
                grid.v[grid_idx] /= m_i;
                grid.vz[grid_idx] /= m_i;
                grid.h[grid_idx] /= m_i;

                TV2 v_old = grid.v[grid_idx];

                grid.a[grid_idx] = grid.F[grid_idx] / m_i;   // Note: Does not contain friction!
                TV2 v_trial = v_old + dt * grid.a[grid_idx]; // Trial velocity without friction

                const T v_sq = v_trial.squaredNorm();

                int terrain_idx = ind(grid_i_offset + i, grid_j_offset + j, terrain.Ny); // Add grid offset to get corresponding terrain indices

                if (v_sq > 1e-16) { // Moving
                    const TV2 grad = terrain.gradients[terrain_idx];
                    const T theta = terrain.theta[terrain_idx];

                    const T cos_theta = cos(theta);

                    const T v_z = v_trial.dot(grad);
                    const T v_3D_sq = v_sq + v_z * v_z;
                    const T v_norm_3D = sqrt(v_3D_sq);

                    T a_fric = 0.; // Slope-parallel (not horizontal)

                    switch (friction.type) {
                        case FrictionType::Coulomb: {
                            a_fric = friction.coulomb.mu * gravity * cos_theta;
                            break;
                        }
                        case FrictionType::Voellmy: {
                            a_fric = friction.voellmy.mu * gravity * cos_theta + gravity * v_3D_sq / (friction.voellmy.xi * grid.h[grid_idx] * cos_theta);
                            break;
                        }
                        case FrictionType::mu_I: {
                            const T I = 5. * friction.muI.d * v_norm_3D / (2. * grid.h[grid_idx] * cos_theta * cos_theta * sqrt(gravity * grid.h[grid_idx] * friction.muI.phi));
                            const T mu_I = friction.muI.mu_1 + (friction.muI.mu_2 - friction.muI.mu_1) / (1. + friction.muI.I0 / I);
                            a_fric = mu_I * gravity * cos_theta;
                            break;
                        }
                    }

                    // Maximum removable speed due to friction. Slope-parallel (not horizontal), therefore division by 3D velocity norm in grid.v instead of division by 2D norm
                    const T dv = dt * a_fric;

                    if (dv < v_norm_3D) { // Reduce velocity by dv in the direction of v_trial
                        grid.v[grid_idx] = v_trial * (1. - dv / v_norm_3D);
                    } else { // If friction would reverse velocity, set velocity to zero
                        grid.v[grid_idx].setZero();
                    }
                } else {                        // Static
                    grid.v[grid_idx].setZero(); // Don't zero grid.vz because the error is only a lag of one time step
                }

                grid.flip[grid_idx] = grid.v[grid_idx] - v_old; // dt * acceleration, including friction

                terrain.h_max[terrain_idx] = std::max(terrain.h_max[terrain_idx], grid.h[grid_idx]);

                const T v_2d_sq = grid.v[grid_idx].squaredNorm();
                terrain.v_max2[terrain_idx] = std::max(terrain.v_max2[terrain_idx], v_2d_sq);
            } else { // If this is not zeroed, non-normalized terms might remain
                grid.v[grid_idx].setZero();
                grid.vz[grid_idx] = 0.;
                grid.h[grid_idx] = 0.;
            }
        }
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::G2P() {
    // Interpolate grid velocity to particles, update particle positions.

    for (int p = 0; p < Np; ++p) {
        TV2& x_p = particles.x[p];
        const TV2 part_grad = particleGradients(p);

        int i_base = static_cast<int>((x_p(0) - grid.x[0]) * one_over_dx - 0.5);
        int j_base = static_cast<int>((x_p(1) - grid.y[0]) * one_over_dx - 0.5);

        // Particle-wise local variables
        TV2 pic_p = TV2::Zero();
        TV2 flip_p = TV2::Zero();
        T z_p = 0.;

        // Precompute weights
        T wx[3], wy[3], dwx[3], dwy[3];
        for (int di = 0; di < 3; ++di) {
            int i = i_base + di;
            T dx_ = (x_p(0) - grid.x[i]) * one_over_dx;
            wx[di] = N(dx_);
            dwx[di] = dNdu(dx_) * one_over_dx;
        }
        for (int dj = 0; dj < 3; ++dj) {
            int j = j_base + dj;
            T dy_ = (x_p(1) - grid.y[j]) * one_over_dx;
            wy[dj] = N(dy_);
            dwy[dj] = dNdu(dy_) * one_over_dx;
        }

        T strain_trace = 0.;

        for (int di = 0; di < 3; ++di) {
            int i = i_base + di;
            for (int dj = 0; dj < 3; ++dj) {
                int j = j_base + dj;

                T weight = wx[di] * wy[dj];
                int grid_idx = ind(i, j, grid.Ny);

                T grad_weight_x = dwx[di] * wy[dj];
                T grad_weight_y = dwy[dj] * wx[di];

                pic_p += grid.v[grid_idx] * weight;
                flip_p += grid.flip[grid_idx] * weight;
                z_p += grid.z[grid_idx] * weight;

                strain_trace += grad_weight_x * grid.v[grid_idx](0) + grad_weight_y * grid.v[grid_idx](1);
            }
        }

        strain_trace *= dt;

        particles.z[p] = z_p;

        particles.v[p] = flip_ratio * (particles.v[p] + flip_p) + (1. - flip_ratio) * pic_p;
        particles.vz[p] = particles.v[p].dot(part_grad);

        x_p += dt * pic_p;

        particles.h[p] /= 1. + strain_trace;                         // h update
        particles.A[p] = particles.mass[p] / (particles.h[p] * rho); // A update
    }
}

template class Simulation<double, Eigen::Matrix<double, 2, 1>>;