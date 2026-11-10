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
#pragma GCC diagnostic push // Save current state of warning settings
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Warray-bounds"
#include "tinyply.h"
#pragma GCC diagnostic pop // Restore previous state of warning settings
#include "gdal_priv.h"
#include <cmath>

namespace fs = std::filesystem;
using std::ios;
using std::ofstream;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::vector;

template <typename T, typename TV2> void Simulation<T, TV2>::saveParticleData() const {
    ostringstream ss;
    ss << "particles_frame_" << setw(padding) << setfill('0') << frame << ".ply";
    fs::path filename = particle_path / ss.str();

    ofstream out(filename, ios::out | ios::binary);
    tinyply::PlyFile file;

    vector<float> x_f(2 * Np);
    vector<float> v_f(3 * Np);
    vector<float> z_f(Np);
    vector<float> h_f(Np);
    vector<float> A_f(Np);

    for (int i = 0; i < Np; ++i) {
        x_f[2 * i + 0] = static_cast<float>(particles.x[i][0]);
        x_f[2 * i + 1] = static_cast<float>(particles.x[i][1]);

        v_f[3 * i + 0] = static_cast<float>(particles.v[i][0]);
        v_f[3 * i + 1] = static_cast<float>(particles.v[i][1]);
        v_f[3 * i + 2] = static_cast<float>(particles.vz[i]);

        z_f[i] = static_cast<float>(particles.z[i]);

        h_f[i] = static_cast<float>(particles.h[i]);

        A_f[i] = static_cast<float>(particles.A[i]);
    }

    file.add_properties_to_element("vertex", {"x", "y"}, tinyply::Type::FLOAT32, Np, reinterpret_cast<uint8_t*>(x_f.data()), tinyply::Type::INVALID, 0);

    file.add_properties_to_element("vertex", {"vx", "vy", "vz"}, tinyply::Type::FLOAT32, Np, reinterpret_cast<uint8_t*>(v_f.data()), tinyply::Type::INVALID, 0);

    file.add_properties_to_element("vertex", {"z"}, tinyply::Type::FLOAT32, Np, reinterpret_cast<uint8_t*>(z_f.data()), tinyply::Type::INVALID, 0);

    file.add_properties_to_element("vertex", {"h"}, tinyply::Type::FLOAT32, Np, reinterpret_cast<uint8_t*>(h_f.data()), tinyply::Type::INVALID, 0);

    file.add_properties_to_element("vertex", {"A"}, tinyply::Type::FLOAT32, Np, reinterpret_cast<uint8_t*>(A_f.data()), tinyply::Type::INVALID, 0);

    file.write(out, true);
}

template <typename T, typename TV2> void Simulation<T, TV2>::saveTerrain() const {
    ofstream outFile(project_dir_path / sim_name / "out_terrain.csv");
    outFile << "x"
            << ","
            << "y"
            << ","
            << "z"
            << ","
            << "theta in °"
            << "\n";

    for (int i = 0; i < terrain.Nx; ++i) {
        for (int j = 0; j < terrain.Ny; ++j) {
            outFile << terrain.x[i] << "," << terrain.y[j] << "," << terrain.z[ind(i, j, terrain.Ny)] << "," << terrain.theta[ind(i, j, terrain.Ny)] * 180. / M_PI << "\n";
        }
    }
}

template <typename T, typename TV2> void Simulation<T, TV2>::saveGridData() const {
    std::ostringstream ss;
    ss << "rasters_frame_" << std::setw(padding) << std::setfill('0') << frame << ".tif";
    fs::path filename = raster_path / ss.str();

    GDALAllRegister();
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver == nullptr) {
        throw std::runtime_error("GeoTIFF driver not available.");
    }

    char** papszOptions = nullptr;
    papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "ZSTD");
    papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "3");
    GDALDataset* poDstDS = poDriver->Create(filename.string().c_str(), grid.Nx, grid.Ny, 6, GDT_Float32, papszOptions);
    CSLDestroy(papszOptions);

    if (poDstDS == nullptr) {
        throw std::runtime_error("Could not create GeoTIFF: " + filename.string());
    }

    // GeoTransform
    T adfGeoTransform[6] = {
        terrain.xllcorner + grid.x.front(), // x-coordinate of top-left corner
        dx,                                 // pixel width
        0,
        terrain.yllcorner + grid.y.back(), // y-coordinate of top-left corner
        0,
        -dx // pixel height (negative for north-up image)
    };

    poDstDS->SetGeoTransform(adfGeoTransform);

    // Buffers for raster bands
    std::vector<float> h_f(grid.Nx * grid.Ny);
    std::vector<float> v_f(grid.Nx * grid.Ny);
    std::vector<float> m_f(grid.Nx * grid.Ny);
    std::vector<float> v_x(grid.Nx * grid.Ny);
    std::vector<float> v_y(grid.Nx * grid.Ny);
    std::vector<float> v_z(grid.Nx * grid.Ny);

    constexpr float nodata = -9999.0f;

    // Copy from grid -> GDAL ordering (row-major, top to bottom)
    for (int j = 0; j < grid.Ny; ++j) {
        for (int i = 0; i < grid.Nx; ++i) {
            int grid_idx = ind(i, j, grid.Ny);
            int dst = (grid.Ny - 1 - j) * grid.Nx + i;

            float hval = static_cast<float>(grid.h[grid_idx]);
            if (hval == 0.0f) {
                h_f[dst] = nodata;
                v_f[dst] = nodata;
                m_f[dst] = nodata;
                v_x[dst] = nodata;
                v_y[dst] = nodata;
                v_z[dst] = nodata;
            } else {
                h_f[dst] = hval;
                const T v_2d_sq = grid.v[grid_idx].squaredNorm();
                v_f[dst] = static_cast<float>(std::sqrt(v_2d_sq));
                m_f[dst] = static_cast<float>(grid.mass[grid_idx]);
                v_x[dst] = static_cast<float>(grid.v[grid_idx][0]);
                v_y[dst] = static_cast<float>(grid.v[grid_idx][1]);
                v_z[dst] = static_cast<float>(grid.vz[grid_idx]);
            }
        }
    }

    // Write band 1 (h)
    GDALRasterBand* poBand = poDstDS->GetRasterBand(1);
    poBand->SetDescription("h");
    poBand->SetNoDataValue(nodata);
    CPLErr err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, h_f.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band h to GeoTIFF");
    }

    // Write band 2 (horizontal velocity norm)
    poBand = poDstDS->GetRasterBand(2);
    poBand->SetDescription("v_norm");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, v_f.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band v_norm to GeoTIFF");
    }

    // Write band 3 (mass)
    poBand = poDstDS->GetRasterBand(3);
    poBand->SetDescription("m");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, m_f.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band m to GeoTIFF");
    }

    // Write band 4 (v_x)
    poBand = poDstDS->GetRasterBand(4);
    poBand->SetDescription("v_x");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, v_x.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band v_x to GeoTIFF");
    }

    // Write band 5 (v_y)
    poBand = poDstDS->GetRasterBand(5);
    poBand->SetDescription("v_y");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, v_y.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band v_y to GeoTIFF");
    }

    // Write band 6 (v_z)
    poBand = poDstDS->GetRasterBand(6);
    poBand->SetDescription("v_z");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, grid.Nx, grid.Ny, v_z.data(), grid.Nx, grid.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band v_z to GeoTIFF");
    }

    GDALClose(poDstDS);
}

template <typename T, typename TV2> void Simulation<T, TV2>::saveMaxData() const {
    fs::path filename = sim_dir_path / "max_values.tif";

    GDALAllRegister();
    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver == nullptr) {
        throw std::runtime_error("GeoTIFF driver not available.");
    }

    char** papszOptions = nullptr;
    papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "ZSTD");
    papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "3");
    GDALDataset* poDstDS = poDriver->Create(filename.string().c_str(), terrain.Nx, terrain.Ny, 2, GDT_Float32, papszOptions);
    CSLDestroy(papszOptions);

    if (poDstDS == nullptr) {
        throw std::runtime_error("Could not create GeoTIFF: " + filename.string());
    }

    // GeoTransform for full terrain extent
    T adfGeoTransform[6] = {terrain.xllcorner, dx, 0, terrain.yllcorner + terrain.y.back(), 0, -dx};

    poDstDS->SetGeoTransform(adfGeoTransform);

    std::vector<float> h_f(terrain.Nx * terrain.Ny);
    std::vector<float> v_f(terrain.Nx * terrain.Ny);
    constexpr float nodata = -9999.0f;

    for (int j = 0; j < terrain.Ny; ++j) {
        for (int i = 0; i < terrain.Nx; ++i) {
            int idx = ind(i, j, terrain.Ny);
            int dst = (terrain.Ny - 1 - j) * terrain.Nx + i;
            float hval = static_cast<float>(terrain.h_max[idx]);
            float vval = static_cast<float>(std::sqrt(terrain.v_max2[idx]));
            if (hval == 0.0f) {
                h_f[dst] = nodata;
                v_f[dst] = nodata;
            } else {
                h_f[dst] = hval;
                v_f[dst] = vval;
            }
        }
    }

    GDALRasterBand* poBand = poDstDS->GetRasterBand(1);
    poBand->SetDescription("h_max");
    poBand->SetNoDataValue(nodata);
    CPLErr err = poBand->RasterIO(GF_Write, 0, 0, terrain.Nx, terrain.Ny, h_f.data(), terrain.Nx, terrain.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band h_max to GeoTIFF");
    }

    poBand = poDstDS->GetRasterBand(2);
    poBand->SetDescription("v_norm_max");
    poBand->SetNoDataValue(nodata);
    err = poBand->RasterIO(GF_Write, 0, 0, terrain.Nx, terrain.Ny, v_f.data(), terrain.Nx, terrain.Ny, GDT_Float32, 0, 0);
    if (err != CE_None) {
        GDALClose(poDstDS);
        throw std::runtime_error("Failed writing raster band v_norm_max to GeoTIFF");
    }

    GDALClose(poDstDS);
}

template class Simulation<double, Eigen::Matrix<double, 2, 1>>;