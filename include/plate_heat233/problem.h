#pragma once

#include <map>
#include <variant>

#include "plate_heat233/mesh.h"

namespace plate_heat233 {

struct Material {
    double k = 0.0;  // thermal conductivity [W/(m K)], must be positive finite
    double Q = 0.0;  // volumetric heat source [W/m^3], must be finite
};

struct DirichletBC {
    double T = 0.0;  // prescribed temperature [K]
};

struct FluxBC {
    double q = 0.0;  // outward heat flux [W/m^2], positive out of the domain
};

struct ConvectionBC {
    double h = 0.0;   // heat transfer coefficient [W/(m^2 K)], positive finite
    double Ta = 0.0;  // ambient temperature [K]
};

using BoundaryCondition = std::variant<DirichletBC, FluxBC, ConvectionBC>;

struct Problem {
    Mesh mesh;
    std::map<int, Material> materials;             // region tag -> material
    std::map<int, BoundaryCondition> boundaries;   // boundary tag -> condition
    double tolerance = 1e-8;  // max allowed absolute residual on free nodes
};

}  // namespace plate_heat233

