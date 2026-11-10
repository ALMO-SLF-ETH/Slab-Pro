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

#include "polygon_tools.hpp"
#include <ogrsf_frmts.h>
#include <stdexcept>

std::vector<TV2> loadPolygonFromShapefile(const std::string& shapefilePath) {
    std::vector<TV2> vertices;

    GDALAllRegister();
    GDALDataset* poDS = (GDALDataset*)GDALOpenEx(shapefilePath.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (!poDS)
        throw std::runtime_error("Failed to open shapefile.");

    OGRLayer* poLayer = poDS->GetLayer(0);
    if (!poLayer) {
        GDALClose(poDS);
        throw std::runtime_error("Failed to get layer.");
    }

    OGRFeature* poFeature = poLayer->GetNextFeature();
    if (!poFeature) {
        GDALClose(poDS);
        throw std::runtime_error("No features found.");
    }

    OGRGeometry* poGeometry = poFeature->GetGeometryRef();
    if (!poGeometry || wkbFlatten(poGeometry->getGeometryType()) != wkbPolygon) {
        OGRFeature::DestroyFeature(poFeature);
        GDALClose(poDS);
        throw std::runtime_error("Not a polygon.");
    }

    OGRPolygon* poPolygon = poGeometry->toPolygon();
    OGRLinearRing* ring = poPolygon->getExteriorRing();
    if (!ring) {
        OGRFeature::DestroyFeature(poFeature);
        GDALClose(poDS);
        throw std::runtime_error("No exterior ring.");
    }

    for (int i = 0; i < ring->getNumPoints(); ++i) {
        vertices.emplace_back(TV2(ring->getX(i), ring->getY(i)));
    }

    OGRFeature::DestroyFeature(poFeature);
    GDALClose(poDS);
    return vertices;
}

T computePolygonArea(const std::vector<TV2>& vertices) {
    if (vertices.size() < 3)
        return static_cast<T>(0); // Not a polygon

    OGRLinearRing ring;
    for (const auto& pt : vertices) {
        ring.addPoint(pt(0), pt(1));
    }

    // Close ring if not closed
    if (!(vertices.front() == vertices.back())) {
        ring.addPoint(vertices.front()(0), vertices.front()(1));
    }

    OGRPolygon ogrPoly;
    ogrPoly.addRing(&ring);

    return static_cast<T>(ogrPoly.get_Area());
}

bool isPointInPolygon(T x, T y, const std::vector<TV2>& polygon) {
    int n = static_cast<int>(polygon.size());
    bool inside = false;

    for (int i = 0, j = n - 1; i < n; j = i++) {
        const T xi = polygon[i].x(), yi = polygon[i].y();
        const T xj = polygon[j].x(), yj = polygon[j].y();

        bool intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-12) + xi);
        if (intersect)
            inside = !inside;
    }

    return inside;
}