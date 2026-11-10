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

#pragma once

#include "data_arrays.hpp"
#include "tools.hpp"
#include "timer.hpp"

template <typename T, typename TV2>

class Simulation {
  public:
    // Constructor with safe defaults
    Simulation()
        : sim_name("default_sim"),
          nFiles(0),
          padding(0),
          rho(static_cast<T>(2500.0)),      // from rhoEdit
          h0(static_cast<T>(30.0)),         // from h0Edit
          gravity(static_cast<T>(9.8045)),  // from gravityEdit
          total_time(static_cast<T>(80.0)), // from totalTimeEdit
          C_cfl(static_cast<T>(0.5)),       // from cflEdit
          dt(static_cast<T>(0.0)),          // initial dt has to be zero for correctly saving time step 0 data
          dt_max(static_cast<T>(0.06)),     // from dtMaxEdit
          flip_ratio(static_cast<T>(0.9)),  // from flipRatioEdit
          output_dt(static_cast<T>(5.0)),   // from outputDtEdit
          frame(0),
          friction({
              FrictionType::Coulomb,      // default type
              {0.531709431759},           // coulomb.mu
              {0.2, 200.0},               // voellmy.mu, voellmy.xi
              {0.2, 0.9, 0.1, 0.005, 0.6} // muI: mu_1, mu_2, I0, d, phi
          }),
          terrain_type(TerrainType::Analytical), // from terrainTypeCombo
          terrain_Lx(static_cast<T>(5000.0)),
          terrain_Ly(static_cast<T>(5000.0)),
          dx(static_cast<T>(8.0)),
          one_over_dx(static_cast<T>(1.0) / 8.0),
          alt_z([](const T&, const T&) {
              return T(0);
          }),
          release_shape_type(ReleaseShapeType::Rectangle),
          rect{
              TV2(static_cast<T>(1200.0), static_cast<T>(1400.0)), // center
              static_cast<T>(100.0),                               // Lx
              static_cast<T>(80.0),                                // Ly
              static_cast<T>(0.0)                                  // rotation
          },
          ell{
              TV2(static_cast<T>(1200.0), static_cast<T>(1400.0)), // center
              static_cast<T>(100.0),                               // Rx
              static_cast<T>(80.0),                                // Ry
              static_cast<T>(0.0)                                  // rotation
          },
          release_height_raster_filepath(""),
          sampling_method(SamplingMethod::Uniform),
          ppc(4),
          Np(0)

    {}
    ~Simulation() = default;

    void loadFromJSON(const std::filesystem::path& folder);

    // setters for simulation configuration
    void setSimName(const std::string& name) {
        sim_name = name;
    }
    void setProjectDir(const std::filesystem::path& path) {
        project_dir_path = path;
    }
    void setParticlePath(const std::filesystem::path& path) {
        particle_path = path;
    }
    void setRasterPath(const std::filesystem::path& path) {
        raster_path = path;
    }
    void setDemFilePath(const std::filesystem::path& path) {
        dem_filepath = path;
    }
    void setPolygonPath(const std::filesystem::path& path) {
        poly.filepath = path;
    }
    void setNFiles(int n) {
        nFiles = n;
    }
    void setPadding(int p) {
        padding = p;
    }
    void setRho(T r) {
        rho = r;
    }
    void setH0(T h) {
        h0 = h;
    }
    void setGravity(T g) {
        gravity = g;
    }
    void setTotalTime(T t) {
        total_time = t;
    }
    void setCFL(T cfl) {
        C_cfl = cfl;
    }
    void setDt(T delta_t) {
        dt = delta_t;
    }
    void setDtMax(T dtm) {
        dt_max = dtm;
    }
    void setFlipRatio(T fr) {
        flip_ratio = fr;
    }
    void setOutputDt(T odt) {
        output_dt = odt;
    }
    void setFrame(int f) {
        frame = f;
    }
    void setFrictionType(const FrictionType& fric) {
        friction.type = fric;
    }
    void setFrictionCoulombMu(double mu) {
        friction.coulomb.mu = mu;
    }
    void setFrictionVoellmyMu(double mu) {
        friction.voellmy.mu = mu;
    }
    void setFrictionVoellmyXi(double xi) {
        friction.voellmy.xi = xi;
    }
    void setFrictionMuIMu1(double mu1) {
        friction.muI.mu_1 = mu1;
    }
    void setFrictionMuIMu2(double mu2) {
        friction.muI.mu_2 = mu2;
    }
    void setFrictionMuII0(double I0) {
        friction.muI.I0 = I0;
    }
    void setFrictionMuID(double d) {
        friction.muI.d = d;
    }
    void setFrictionMuIPhi(double p) {
        friction.muI.phi = p;
    }
    void setTerrainType(TerrainType tt) {
        terrain_type = tt;
    }
    void setTerrainLx(T l1) {
        terrain_Lx = l1;
    }
    void setTerrainLy(T l2) {
        terrain_Ly = l2;
    }
    void setDx(T delta_x) {
        dx = delta_x;
        one_over_dx = T(1) / dx;
    }
    void setOneOverDx(T oodx) {
        one_over_dx = oodx;
        dx = T(1) / one_over_dx;
    }
    void setAltZFunction(const std::function<T(const T&, const T&)>& func) {
        alt_z = func;
    }
    void setReleaseShapeType(ReleaseShapeType rst) {
        release_shape_type = rst;
    }
    void setRectangle(const Rectangle& r) {
        rect = r;
    }
    void setEllipse(const Ellipse& e) {
        ell = e;
    }
    void setPolygon(const Polygon& p) {
        poly = p;
    }
    void setReleaseArea(T area) {
        release_area = area;
    }
    void setReleaseHeightRasterFilepath(const std::string& filepath) {
        release_height_raster_filepath = filepath;
    }
    void setSamplingMethod(SamplingMethod sm) {
        sampling_method = sm;
    }
    void setNp(int np) {
        Np = np;
    }
    void setParticles(const Particles& parts) {
        particles = parts;
    }
    void setGrid(const Grid& g) {
        grid = g;
    }
    void setTerrain(const Grid& t) {
        terrain = t;
    }

    std::filesystem::path getProjectDir() const {
        return project_dir_path;
    }
    std::string getSimName() const {
        return sim_name;
    }
    std::filesystem::path getDemFilePath() const {
        return dem_filepath;
    }
    std::filesystem::path getParticlePath() const {
        return particle_path;
    }
    std::filesystem::path getRasterPath() const {
        return raster_path;
    }

    FrictionType getFrictionType() const {
        return friction.type;
    }

    // Main functions
    // Simulation lifecycle
    void initialize();
    void initializeTerrain();
    void initializeParticles();
    void simulate();

    // Particle/grid computations
    TV2 rotatePoint(const TV2& p, T angle) const;
    TV2 particleGradients(const int p) const;
    T particleSlopeAngle(const int p) const;
    void terrainGradients();
    void terrainTheta();
    std::vector<T> generateGrid(T length, T dx) const;
    T getReleaseHeightFromRaster(const TV2 position) const;

    // I/O
    void saveTerrain() const;
    void saveParticleData() const;
    void saveGridData() const;
    void saveMaxData() const;
    RasterData readASC(const std::filesystem::path& path) const;

    // Prevent unsafe copy/move
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) = default;
    Simulation& operator=(Simulation&&) = default;

    // Sort
    void sort_particles_by_grid_block();
    void apply_particle_permutation(const std::vector<int>& order);

  private:
    // Helper to convert friction type string to enum
    FrictionType parseFrictionType(const std::string& str);

    // Helper to convert terrain type string to enum
    TerrainType parseTerrainType(const std::string& str);

    // Helper to convert release shape string to enum
    ReleaseShapeType parseReleaseShapeType(const std::string& str);

    // Internal step functions
    void advanceStep(int time_step);
    void updateDt();
    void remesh();
    void P2G();
    void computeGridVelocity();
    void G2P();

    // Simulation info
    std::string sim_name;
    std::filesystem::path json_file_path;
    std::filesystem::path sim_dir_path;
    std::filesystem::path project_dir_path;
    std::filesystem::path particle_path;
    std::filesystem::path raster_path;
    std::filesystem::path dem_filepath;
    int nFiles;
    int padding;

    // Physical parameters
    T rho;
    T h0;
    T gravity;
    T time;
    T total_time;
    T C_cfl;
    T dt;
    T dt_max;
    T flip_ratio;
    T output_dt;
    int frame;

    // Friction
    Friction friction;

    // Terrain
    TerrainType terrain_type;
    T terrain_Lx;
    T terrain_Ly;
    T dx;
    T one_over_dx;
    std::function<T(const T&, const T&)> alt_z;

    // Release shape
    ReleaseShapeType release_shape_type;
    Rectangle rect;
    Ellipse ell;
    Polygon poly;
    Raster rast;
    T release_area;
    std::string release_height_raster_filepath;
    Grid release_height_raster;

    // Sampling
    SamplingMethod sampling_method;
    const int ppc;

    // Particle/grid data
    int Np;
    Particles particles;
    Grid grid;
    Grid terrain;

    // Timer
    T runtime_updateDt = 0.;
    T runtime_sort = 0.;
    T runtime_remesh = 0.;
    T runtime_P2G = 0.;
    T runtime_computeGridVelocity = 0.;
    T runtime_G2P = 0.;

    // Offset of grid relative to terrain
    int grid_i_offset = 0;
    int grid_j_offset = 0;
};