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

#include <chrono>

class timer {
  public:
    timer() = default;

    inline void start() {
        start_time = clock::now();
    }

    inline void stop() {
        stop_time = clock::now();
    }

    double get_timing() const {
        return std::chrono::duration<double>(stop_time - start_time).count();
    }

  private:
    using clock = std::chrono::steady_clock;
    clock::time_point start_time, stop_time;
};