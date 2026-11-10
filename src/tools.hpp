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

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <filesystem>

/// DOUBLE TYPE DEFINITION
typedef double T;

//////////////////////// EIGEN objects
typedef Eigen::Matrix<T, 2, 1> TV2;
typedef Eigen::Matrix<T, 2, 2> TM2;

// Strongly typed enums
enum class TerrainType { Analytical, DEM };
enum class ReleaseShapeType { Rectangle, Ellipse, Polygon, Raster };
enum class SamplingMethod { Uniform, Poisson };
enum class FrictionType { Coulomb, Voellmy, mu_I };
struct Rectangle {
    TV2 center;
    T Lx;
    T Ly;
    T rotation;
};

struct Ellipse {
    TV2 center;
    T rx;
    T ry;
    T rotation;
};

struct Polygon {
    std::filesystem::path filepath;
    std::vector<TV2> vertices;
};

struct Raster {
    std::filesystem::path filepath;
    T h_min;
};

struct CoulombFriction {
    double mu;
};

struct VoellmyFriction {
    double mu;
    double xi;
};

struct MuIFriction {
    double mu_1;
    double mu_2;
    double I0;
    double d;
    double phi;
};

struct Friction {
    FrictionType type;
    CoulombFriction coulomb;
    VoellmyFriction voellmy;
    MuIFriction muI;
};

std::filesystem::path createSimulationDirectories(const std::string& simName, const std::filesystem::path& baseDirectory);

///////////////////// TOOLS ////////////////////////

inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// // cubic shape function
// inline T N(const T& u) {
//     T uabs = std::abs(u);
//     if (uabs < 1.0) {
//         return 0.5 * uabs * uabs * uabs - u * u + 0.6666666666666666666666666666666666666666666666666;
//     } else if (uabs < 2.0) {
//         return 0.1666666666666666666666666666666666666666666 * (2.0 - uabs) * (2.0 - uabs) * (2.0 - uabs);
//     } else {
//         return 0.0;
//     }
// }

// // cubic shape function derivative
// inline T dNdu(const T& u) {
//     T uabs = std::abs(u);
//     if (uabs < 1.0) {
//         return u * (1.5 * uabs - 2.0);
//     } else if (uabs < 2.0) {
//         return -0.5 * sgn(u) * (2.0 - uabs) * (2.0 - uabs);
//     } else {
//         return 0.0;
//     }
// }

// quadratic shape function
inline T N(const T& u) {
    T uabs = std::abs(u);
    if (uabs < 0.5) {
        return 0.75 - uabs * uabs;
    } else if (uabs < 1.5) {
        return 0.5 * (1.5 - uabs) * (1.5 - uabs);
    } else {
        return 0.0;
    }
}

// quadratic shape function derivative
inline T dNdu(const T& u) {
    T uabs = std::abs(u);
    if (uabs < 0.5) {
        return (-2.0 * u);
    } else if (uabs < 1.5) {
        return (u - 1.5 * sgn(u));
    } else {
        return 0.0;
    }
}

// linear shape function
// inline T N(const T& u) {
//     T uabs = std::abs(u);
//     if (uabs < 1.0) {
//         return 1.0 - uabs;
//     } else {
//         return 0.0;
//     }
// }

// linear shape function derivative
// inline T dNdu(const T& u) {
//     T uabs = std::abs(u);
//     if (uabs < 1.0) {
//         return -sgn(u);
//     } else {
//         return 0.0;
//     }
// }

inline T wip(const T& xp, const T& yp, const T& xi, const T& yi, const T& one_over_h) {
    return N((xp - xi) * one_over_h) * N((yp - yi) * one_over_h);
}

inline TV2 grad_wip(const T& xp, const T& yp, const T& xi, const T& yi, const T& one_over_h) {
    TV2 out;
    out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h, dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * one_over_h;
    return out;
}

inline int ind(const int& i, const int& j, const int& Ny) {
    return i * Ny + j;
}