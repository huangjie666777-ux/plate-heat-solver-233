// Dual-material heat spreader example for plate_heat233.
//
// Plate [0,1] x [0,0.5] m, unit thickness, split at x = 0.5:
//   region 1 (left,  k = 50 W/(m K), no source)
//   region 2 (right, k = 5  W/(m K), Q = 2e4 W/m^3)
// Boundary tags: 1 = left (T = 320 K), 2 = right (convection h = 25, Ta = 300 K),
// 3 = top, 4 = bottom (both left unspecified -> adiabatic).

#include "plate_heat233/plate_heat233.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace {

void writeDualMaterialMesh(const std::string& path, int nx, int ny) {
  std::ofstream out(path);
  if (!out) throw plate_heat233::Error("cannot write mesh file " + path);
  const double lx = 1.0, ly = 0.5;
  const int nodesX = nx + 1, nodesY = ny + 1;

  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  out << "$Nodes\n" << nodesX * nodesY << "\n";
  auto nodeId = [&](int i, int j) { return j * nodesX + i + 1; };
  for (int j = 0; j < nodesY; ++j)
    for (int i = 0; i < nodesX; ++i)
      out << nodeId(i, j) << " " << lx * i / nx << " " << ly * j / ny << " 0\n";

  const int nTri = 2 * nx * ny;
  const int nLine = 2 * nx + 2 * ny;
  out << "$Elements\n" << nTri + nLine << "\n";
  long eid = 1;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const int region = (i + 0.5) / nx * lx < 0.5 ? 1 : 2;
      const long n00 = nodeId(i, j), n10 = nodeId(i + 1, j);
      const long n01 = nodeId(i, j + 1), n11 = nodeId(i + 1, j + 1);
      out << eid++ << " 2 1 " << region << " " << n00 << " " << n10 << " " << n11 << "\n";
      out << eid++ << " 2 1 " << region << " " << n00 << " " << n11 << " " << n01 << "\n";
    }
  }
  for (int j = 0; j < ny; ++j)  // left edge, tag 1
    out << eid++ << " 1 1 1 " << nodeId(0, j) << " " << nodeId(0, j + 1) << "\n";
  for (int j = 0; j < ny; ++j)  // right edge, tag 2
    out << eid++ << " 1 1 2 " << nodeId(nx, j) << " " << nodeId(nx, j + 1) << "\n";
  for (int i = 0; i < nx; ++i)  // bottom edge, tag 4 (adiabatic)
    out << eid++ << " 1 1 4 " << nodeId(i, 0) << " " << nodeId(i + 1, 0) << "\n";
  for (int i = 0; i < nx; ++i)  // top edge, tag 3 (adiabatic)
    out << eid++ << " 1 1 3 " << nodeId(i, ny) << " " << nodeId(i + 1, ny) << "\n";
  out << "$EndElements\n";
}

}  // namespace

int main() {
  using namespace plate_heat233;
  try {
    const std::string meshPath = "dual_material.msh";
    writeDualMaterialMesh(meshPath, 40, 20);
    Mesh mesh = loadGmsh(meshPath);

    SolveOptions opts;
    opts.materials[1] = Material{50.0, 0.0};
    opts.materials[2] = Material{5.0, 2.0e4};
    opts.boundaries[1] = DirichletBC{320.0};
    opts.boundaries[2] = ConvectionBC{25.0, 300.0};
    opts.residualTolerance = 1e-8;

    SolveResult res = solve(mesh, opts);

    double tMin = 1e300, tMax = -1e300;
    for (const auto& [id, t] : res.temperature) {
      tMin = std::min(tMin, t);
      tMax = std::max(tMax, t);
    }
    std::cout << "nodes:            " << mesh.nodes.size() << "\n";
    std::cout << "triangles:        " << mesh.triangles.size() << "\n";
    std::cout << "T min / max [K]:  " << tMin << " / " << tMax << "\n";
    std::cout << "max residual:     " << res.maxResidual << "\n";
    std::cout << "source power:     " << res.power.source << " W/m\n";
    std::cout << "dirichlet out:    " << res.power.dirichletOut << " W/m\n";
    std::cout << "neumann out:      " << res.power.neumannOut << " W/m\n";
    std::cout << "convection out:   " << res.power.convectionOut << " W/m\n";
    std::cout << "power imbalance:  " << res.power.imbalance << " W/m\n";

    writeVtk("dual_material.vtk", mesh, res);
    std::cout << "wrote dual_material.vtk\n";
  } catch (const Error& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
