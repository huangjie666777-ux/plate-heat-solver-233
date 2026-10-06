#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

namespace plate_heat233 {

struct Result {
    bool success = false;
    std::string error;

    std::map<int, double> temperature;              // node id -> T [K]
    std::vector<std::array<double, 2>> heat_flux;   // per triangle: -k grad T [W/m^2]
    double max_residual = 0.0;                      // max |residual| on free nodes
    double total_source_power = 0.0;                // integral of Q over domain [W]
    std::map<int, double> boundary_power;           // boundary tag -> outward power [W]
    double power_imbalance = 0.0;                   // |sources - outward boundary sum| [W]
};

}  // namespace plate_heat233

