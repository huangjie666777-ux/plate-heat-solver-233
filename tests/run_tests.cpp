// Self-tests for plate_heat233. Returns non-zero if any check fails.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "plate_heat233/plate_heat233.h"

using namespace plate_heat233;

namespace {

int failures = 0;

void check(bool ok, const std::string& name) {
    std::cout << (ok ? "[PASS] " : "[FAIL] ") << name << "\n";
    if (!ok) ++failures;
}

std::string replace_first(std::string s, const std::string& from, const std::string& to) {
    std::size_t pos = s.find(from);
    if (pos == std::string::npos) {
        std::cerr << "test setup error: pattern not found: " << from << "\n";
        std::exit(2);
    }
    return s.replace(pos, from.size(), to);
}

void write_file(const std::string& path, const std::string& content) {
    std::ofstream(path) << content;
}

// Simple 2-triangle unit square mesh. Region tag 1, boundary tags 10..13.
const char* kSquareMesh =
    "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    "$Nodes\n4\n"
    "1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 0\n"
    "$EndNodes\n"
    "$Elements\n6\n"
    "1 2 1 1 1 2 3\n"
    "2 2 1 1 1 3 4\n"
    "3 1 1 10 1 2\n"
    "4 1 1 11 2 3\n"
    "5 1 1 12 3 4\n"
    "6 1 1 13 4 1\n"
    "$EndElements\n";

bool mesh_fails(const std::string& content, const std::string& needle) {
    write_file("tests/tmp.msh", content);
    try {
        Mesh::read_gmsh("tests/tmp.msh");
    } catch (const MeshError& e) {
        return std::string(e.what()).find(needle) != std::string::npos;
    }
    return false;
}

}  // namespace

int main() {
    // ---- Mesh validation ----
    check(mesh_fails(replace_first(kSquareMesh, "\n2 1 0 0\n", "\n1 1 0 0\n"),
                     "duplicate node id"),
          "duplicate node id rejected");

    check(mesh_fails(replace_first(kSquareMesh, "1 2 1 1 1 2 3", "1 2 1 1 1 2 9"),
                     "unknown node"),
          "unknown node reference rejected");

    // Node 3 moved onto the x axis: triangle (1,2,3) becomes collinear.
    check(mesh_fails(replace_first(kSquareMesh, "\n3 1 1 0\n", "\n3 2 0 0\n"),
                     "degenerate triangle"),
          "degenerate triangle rejected");

    // Third triangle sharing edge 1-2 -> non-manifold.
    {
        std::string m = replace_first(kSquareMesh, "$Elements\n6\n", "$Elements\n7\n");
        m = replace_first(m, "6 1 1 13 4 1\n", "6 1 1 13 4 1\n7 2 1 1 2 1 3\n");
        check(mesh_fails(m, "non-manifold"), "non-manifold edge rejected");
    }

    // Boundary line on the interior diagonal 1-3.
    {
        std::string m = replace_first(kSquareMesh, "$Elements\n6\n", "$Elements\n7\n");
        m = replace_first(m, "6 1 1 13 4 1\n", "6 1 1 13 4 1\n7 1 1 14 1 3\n");
        check(mesh_fails(m, "not on the outer boundary"),
              "interior boundary line rejected");
    }

    // Non-zero z coordinate.
    check(mesh_fails(replace_first(kSquareMesh, "\n2 1 0 0\n", "\n2 1 0 0.5\n"),
                     "non-zero z"),
          "non-zero z coordinate rejected");

    // ---- Physics on the unit square ----
    write_file("tests/tmp.msh", kSquareMesh);
    Mesh mesh = Mesh::read_gmsh("tests/tmp.msh");

    // Missing material.
    {
        Problem p;
        p.mesh = mesh;
        p.boundaries[10] = DirichletBC{300.0};
        check(!solve(p).success, "missing material fails");
    }

    // Conflicting Dirichlet at shared node (tags 10 and 11 share node 2).
    {
        Problem p;
        p.mesh = mesh;
        p.materials[1] = Material{1.0, 0.0};
        p.boundaries[10] = DirichletBC{300.0};
        p.boundaries[11] = DirichletBC{400.0};
        Result r = solve(p);
        check(!r.success && r.error.find("conflicting") != std::string::npos,
              "conflicting Dirichlet fails");
    }

    // No Dirichlet / convection anywhere -> not unique.
    {
        Problem p;
        p.mesh = mesh;
        p.materials[1] = Material{1.0, 0.0};
        p.boundaries[10] = FluxBC{0.0};
        check(!solve(p).success, "pure Neumann fails (temperature not unique)");
    }

    // 1D conduction: T(x=0)=300, convection at x=1 (h=10, Ta=290), k=1, Q=0.
    // Exact linear profile; flux = 300 - T1 = 10 (T1 - 290) -> T1 = 290.9090...
    {
        Problem p;
        p.mesh = mesh;
        p.materials[1] = Material{1.0, 0.0};
        p.boundaries[13] = DirichletBC{300.0};         // x = 0
        p.boundaries[11] = ConvectionBC{10.0, 290.0};  // x = 1
        p.tolerance = 1e-10;
        Result r = solve(p);
        double exact = (300.0 + 10.0 * 290.0) / 11.0;
        check(r.success && std::abs(r.temperature.at(2) - exact) < 1e-9 &&
                  std::abs(r.temperature.at(3) - exact) < 1e-9,
              "1D conduction + convection matches analytic solution");
        double q = 300.0 - exact;  // outward power per unit depth
        check(std::abs(r.boundary_power.at(11) - q) < 1e-9 &&
                  std::abs(r.boundary_power.at(13) + q) < 1e-9 &&
                  r.power_imbalance < 1e-9,
              "power balance and Dirichlet reaction correct");
    }

    // Uniform source Q=2, k=1, T(x=0)=0, adiabatic elsewhere. On this
    // 2-triangle mesh the consistent P1 solution is T2 = 8/9, T3 = 10/9
    // (the continuum limit T(x) = Q x (2 - x) / (2k) gives T(1) = 1).
    {
        Problem p;
        p.mesh = mesh;
        p.materials[1] = Material{1.0, 2.0};
        p.boundaries[13] = DirichletBC{0.0};
        p.tolerance = 1e-10;
        Result r = solve(p);
        check(r.success && std::abs(r.temperature.at(2) - 8.0 / 9.0) < 1e-9 &&
                  std::abs(r.temperature.at(3) - 10.0 / 9.0) < 1e-9,
              "uniform source conduction matches reference FEM solution");
        check(std::abs(r.total_source_power - 2.0) < 1e-12 &&
                  std::abs(r.boundary_power.at(13) - 2.0) < 1e-9,
              "source power and reaction balance");
    }

    // Outward flux q=5 at x=1, T(x=0)=100, k=1 -> T(1) = 95.
    {
        Problem p;
        p.mesh = mesh;
        p.materials[1] = Material{1.0, 0.0};
        p.boundaries[13] = DirichletBC{100.0};
        p.boundaries[11] = FluxBC{5.0};
        p.tolerance = 1e-10;
        Result r = solve(p);
        check(r.success && std::abs(r.temperature.at(2) - 95.0) < 1e-9,
              "outward flux boundary matches analytic solution");
    }

    std::remove("tests/tmp.msh");
    std::cout << (failures == 0 ? "all tests passed" : "FAILURES PRESENT") << "\n";
    return failures == 0 ? 0 : 1;
}
