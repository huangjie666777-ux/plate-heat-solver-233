#pragma once

#include <string>

#include "plate_heat233/mesh.h"
#include "plate_heat233/result.h"

namespace plate_heat233 {

// Write an ASCII VTK (legacy) unstructured grid with the temperature field
// (point data) and the per-triangle heat flux vector (cell data).
// Returns false and sets error on failure.
bool write_vtk(const std::string& path, const Mesh& mesh, const Result& result,
               std::string& error);

}  // namespace plate_heat233

