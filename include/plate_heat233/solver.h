#pragma once

#include "plate_heat233/problem.h"
#include "plate_heat233/result.h"

namespace plate_heat233 {

// Solve -div(k grad T) = Q on the unit-thickness plate with P1 triangles.
// Never throws for physical/modelling errors; reports them in Result.error.
Result solve(const Problem& problem);

}  // namespace plate_heat233

