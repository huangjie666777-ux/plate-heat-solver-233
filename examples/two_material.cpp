// Two-material plate example for plate_heat233.
//
// A 0.1 m x 0.1 m plate, split at x = 0.05 m:
//   region 1 (left):  k = 200 W/(m K), Q = 0
//   region 2 (right): k = 50  W/(m K), Q = 1e5 W/m^3
//   boundary 10 (left edge):  Dirichlet T = 300 K
//   boundary 11 (right edge): convection h = 500 W/(m^2 K), Ta = 290 K
//   top/bottom edges:         adiabatic (no boundary condition assigned)
// The example generates a Gmsh 2.2 ASCII mesh (non-contiguous node ids,
// mixed triangle windings), solves, prints the power balance and exports
// an ASCII VTK file.

#include <fstream>
#include <iostream>
#include <string>

#include "plate_heat233/plate_heat233.h"

using namespace plate_heat233;

namespace {

void write_mesh(const std::string& path, int nx, int ny) {
    const double lx = 0.1, ly = 0.1;
    const int nx1 = nx / 2;  // elements of region 1 along x
    auto node_id = [&](int i, int j) { return 1000 + j * (nx + 1) + i; };

    std::ofstream out(path);
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    out << "$PhysicalNames\n5\n"
        << "2 1 \"copper_half\"\n"
        << "2 2 \"steel_half_heated\"\n"
        << "1 10 \"left_fixed\"\n"
        << "1 11 \"right_convection\"\n"
        << "1 12 \"top_bottom\"\n"
        << "$EndPhysicalNames\n";
    out << "$Nodes\n" << (nx + 1) * (ny + 1) << "\n";
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            out << node_id(i, j) << " " << lx * i / nx << " " << ly * j / ny << " 0\n";
    out << "$EndNodes\n";

    const int ntri = 2 * nx * ny;
    const int nline = 2 * nx + 2 * ny;
    out << "$Elements\n" << ntri + nline << "\n";
    int eid = 1;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int region = (i < nx1) ? 1 : 2;
            int n00 = node_id(i, j), n10 = node_id(i + 1, j);
            int n01 = node_id(i, j + 1), n11 = node_id(i + 1, j + 1);
            // Alternate winding to prove orientation independence.
            if ((i + j) % 2 == 0) {
                out << eid++ << " 2 1 " << region << " " << n00 << " " << n10 << " " << n11 << "\n";
                out << eid++ << " 2 1 " << region << " " << n00 << " " << n11 << " " << n01 << "\n";
            } else {
                out << eid++ << " 2 1 " << region << " " << n00 << " " << n11 << " " << n10 << "\n";
                out << eid++ << " 2 1 " << region << " " << n00 << " " << n01 << " " << n11 << "\n";
            }
        }
    }
    for (int j = 0; j < ny; ++j)  // left edge, tag 10
        out << eid++ << " 1 1 10 " << node_id(0, j) << " " << node_id(0, j + 1) << "\n";
    for (int j = 0; j < ny; ++j)  // right edge, tag 11
        out << eid++ << " 1 1 11 " << node_id(nx, j) << " " << node_id(nx, j + 1) << "\n";
    for (int i = 0; i < nx; ++i) {  // bottom and top edges, tag 12 (adiabatic)
        out << eid++ << " 1 1 12 " << node_id(i, 0) << " " << node_id(i + 1, 0) << "\n";
        out << eid++ << " 1 1 12 " << node_id(i, ny) << " " << node_id(i + 1, ny) << "\n";
    }
    out << "$EndElements\n";
}

}  // namespace

int main() {
    const std::string mesh_path = "two_material.msh";
    const std::string vtk_path = "two_material.vtk";

    write_mesh(mesh_path, 40, 40);

    Problem problem;
    try {
        problem.mesh = Mesh::read_gmsh(mesh_path);
    } catch (const MeshError& e) {
        std::cerr << "mesh error: " << e.what() << "\n";
        return 1;
    }
    problem.materials[1] = Material{200.0, 0.0};
    problem.materials[2] = Material{50.0, 1.0e5};
    problem.boundaries[10] = DirichletBC{300.0};
    problem.boundaries[11] = ConvectionBC{500.0, 290.0};
    // Tag 12 intentionally left unspecified: adiabatic.
    problem.tolerance = 1e-8;

    Result result = solve(problem);
    if (!result.success) {
        std::cerr << "solve failed: " << result.error << "\n";
        return 1;
    }

    double tmin = 1e300, tmax = -1e300;
    for (const auto& kv : result.temperature) {
        tmin = std::min(tmin, kv.second);
        tmax = std::max(tmax, kv.second);
    }
    std::cout << "nodes:            " << problem.mesh.nodes().size() << "\n";
    std::cout << "triangles:        " << problem.mesh.triangles().size() << "\n";
    std::cout << "T min / max [K]:  " << tmin << " / " << tmax << "\n";
    std::cout << "max residual:     " << result.max_residual << "\n";
    std::cout << "source power [W]: " << result.total_source_power << "\n";
    for (const auto& kv : result.boundary_power)
        std::cout << "boundary " << kv.first << " outward power [W]: " << kv.second << "\n";
    std::cout << "power imbalance [W]: " << result.power_imbalance << "\n";

    std::string error;
    if (!write_vtk(vtk_path, problem.mesh, result, error)) {
        std::cerr << "VTK export failed: " << error << "\n";
        return 1;
    }
    std::cout << "wrote " << vtk_path << "\n";
    return 0;
}

