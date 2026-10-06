#pragma once

#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace plate_heat233 {

class Error : public std::runtime_error {
public:
  explicit Error(const std::string& msg) : std::runtime_error(msg) {}
};

struct Material {
  double k = 0.0;  // thermal conductivity [W/(m*K)], must be positive finite
  double Q = 0.0;  // volumetric heat source [W/m^3] (per unit thickness), finite
};

struct DirichletBC { double T = 0.0; };      // prescribed temperature [K]
struct NeumannBC { double q = 0.0; };        // outward heat flux [W/m^2], positive outward
struct ConvectionBC { double h = 0.0; double Ta = 0.0; };  // q_out = h*(T - Ta), h > 0

using BoundaryCondition = std::variant<DirichletBC, NeumannBC, ConvectionBC>;

struct Node {
  long id = 0;
  double x = 0.0;
  double y = 0.0;
};

struct Triangle {
  long id = 0;
  int tag = 0;  // region tag (first Gmsh physical tag)
  std::array<long, 3> nodes{};
};

struct BoundaryLine {
  long id = 0;
  int tag = 0;  // boundary tag (first Gmsh physical tag)
  std::array<long, 2> nodes{};
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<Triangle> triangles;
  std::vector<BoundaryLine> lines;
};

struct SolveOptions {
  std::map<int, Material> materials;           // region tag -> material
  std::map<int, BoundaryCondition> boundaries; // boundary tag -> condition (missing tag = adiabatic)
  double residualTolerance = 1e-8;             // max absolute residual on free nodes
};

struct TriangleHeatFlux {
  long elementId = 0;
  int tag = 0;
  double qx = 0.0;  // -k * dT/dx [W/m^2]
  double qy = 0.0;  // -k * dT/dy [W/m^2]
};

struct PowerBalance {
  double source = 0.0;        // integral of Q over the domain [W/m]
  double dirichletOut = 0.0;  // power leaving through Dirichlet boundaries (from reactions)
  double neumannOut = 0.0;    // power leaving through prescribed-flux boundaries
  double convectionOut = 0.0; // power leaving through convection boundaries
  double imbalance = 0.0;     // source - (dirichletOut + neumannOut + convectionOut)
};

struct SolveResult {
  std::map<long, double> temperature;      // original node id -> T [K]
  std::vector<TriangleHeatFlux> heatFlux;  // per-triangle -k grad T
  double maxResidual = 0.0;                // max |residual| over free nodes
  PowerBalance power;
};

// Load a Gmsh 2.2 ASCII 2D mesh (z = 0). Element types: 1 (2-node line),
// 2 (3-node triangle). The first tag of each element is the physical tag.
Mesh loadGmsh(const std::string& path);

// Validate mesh topology and limits. Throws Error on any violation.
void validateMesh(const Mesh& mesh);

// Assemble and solve -div(k grad T) = Q with P1 triangles.
// Throws Error on invalid input, non-unique temperature, or solver failure.
SolveResult solve(const Mesh& mesh, const SolveOptions& options);

// Export legacy ASCII VTK with node temperature and per-triangle heat flux.
void writeVtk(const std::string& path, const Mesh& mesh, const SolveResult& result);

}  // namespace plate_heat233
