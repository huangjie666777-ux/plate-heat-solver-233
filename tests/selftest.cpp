// Self-tests for plate_heat233: physics, parser, and error handling.

#include "plate_heat233/plate_heat233.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace plate_heat233;

static int g_failures = 0;

#define CHECK(cond, msg)                                       \
  do {                                                         \
    if (!(cond)) {                                             \
      ++g_failures;                                            \
      std::cerr << "FAIL: " << msg << " (line " << __LINE__ << ")\n"; \
    }                                                          \
  } while (0)

// Unit square, two triangles. Nodes: 1(0,0) 2(1,0) 3(1,1) 4(0,1).
// Boundary tags: 1=bottom, 2=right, 3=top, 4=left.
static Mesh unitSquare(int regionTag = 10) {
  Mesh m;
  m.nodes = {{1, 0, 0}, {2, 1, 0}, {3, 1, 1}, {4, 0, 1}};
  m.triangles = {{101, regionTag, {1, 2, 3}}, {102, regionTag, {1, 3, 4}}};
  m.lines = {{201, 1, {1, 2}}, {202, 2, {2, 3}}, {203, 3, {3, 4}}, {204, 4, {4, 1}}};
  return m;
}

static void testDirichletLinear() {
  SolveOptions opts;
  opts.materials[10] = Material{1.0, 0.0};
  opts.boundaries[4] = DirichletBC{0.0};
  opts.boundaries[2] = DirichletBC{10.0};
  SolveResult r = solve(unitSquare(), opts);
  CHECK(std::abs(r.temperature.at(1) - 0.0) < 1e-12, "T node1");
  CHECK(std::abs(r.temperature.at(2) - 10.0) < 1e-12, "T node2");
  CHECK(std::abs(r.temperature.at(3) - 10.0) < 1e-12, "T node3");
  CHECK(std::abs(r.temperature.at(4) - 0.0) < 1e-12, "T node4");
  CHECK(r.maxResidual <= 1e-8, "residual within tolerance");
  for (const auto& f : r.heatFlux) {
    CHECK(std::abs(f.qx + 10.0) < 1e-9, "flux qx = -k dT/dx");
    CHECK(std::abs(f.qy) < 1e-9, "flux qy = 0");
  }
  CHECK(std::abs(r.power.dirichletOut) < 1e-9, "net dirichlet power ~ 0");
  CHECK(std::abs(r.power.imbalance) < 1e-9, "power balance");
}

static void testConvectionSource() {
  SolveOptions opts;
  opts.materials[10] = Material{2.0, 100.0};
  for (int tag = 1; tag <= 4; ++tag) opts.boundaries[tag] = ConvectionBC{10.0, 20.0};
  SolveResult r = solve(unitSquare(), opts);
  CHECK(r.maxResidual <= 1e-8, "convection case residual");
  CHECK(std::abs(r.power.source - 100.0) < 1e-9, "source power = Q*A");
  CHECK(std::abs(r.power.convectionOut - 100.0) < 1e-6, "convection carries source");
  CHECK(std::abs(r.power.imbalance) < 1e-6, "convection power balance");
  for (const auto& [id, temp] : r.temperature)
    CHECK(temp > 21.0 && temp < 25.0, "temperature in physical range");
}

static void testNeumannDirichlet() {
  // Left T=0, right outward flux q=5 with k=2 -> dT/dx = -q_out/k... check balance only.
  SolveOptions opts;
  opts.materials[10] = Material{2.0, 0.0};
  opts.boundaries[4] = DirichletBC{0.0};
  opts.boundaries[2] = NeumannBC{5.0};  // outward flux +5 on right edge
  SolveResult r = solve(unitSquare(), opts);
  // Outward flux on right = -k dT/dx = 5 -> dT/dx = -2.5 -> T(right) = -2.5
  CHECK(std::abs(r.temperature.at(2) + 2.5) < 1e-9, "neumann sign convention");
  CHECK(std::abs(r.power.neumannOut - 5.0) < 1e-12, "neumann power");
  CHECK(std::abs(r.power.imbalance) < 1e-9, "neumann power balance");
}

static void testParser() {
  const std::string path = "selftest_mesh.msh";
  {
    std::ofstream out(path);
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    out << "$PhysicalNames\n2\n2 10 \"plate\"\n1 2 \"right\"\n$EndPhysicalNames\n";
    out << "$Nodes\n4\n";
    out << "10 0 0 0\n20 1 0 0\n30 1 1 0\n40 0 1 0\n$EndNodes\n";
    out << "$Elements\n6\n";
    out << "100 2 2 10 7 10 20 30\n";      // extra elementary tag, normal winding
    out << "101 2 1 10 10 40 30\n";         // reversed winding
    out << "200 1 1 1 10 20\n";
    out << "201 1 1 2 20 30\n";
    out << "202 1 1 3 30 40\n";
    out << "203 1 1 4 40 10\n";
    out << "$EndElements\n";
  }
  Mesh m = loadGmsh(path);
  CHECK(m.nodes.size() == 4 && m.triangles.size() == 2 && m.lines.size() == 4, "parser counts");
  CHECK(m.triangles[0].tag == 10, "first physical tag used");
  SolveOptions opts;
  opts.materials[10] = Material{1.0, 0.0};
  opts.boundaries[4] = DirichletBC{0.0};
  opts.boundaries[2] = DirichletBC{10.0};
  SolveResult r = solve(m, opts);
  CHECK(std::abs(r.temperature.at(20) - 10.0) < 1e-12, "parser solve T");
  CHECK(std::abs(r.heatFlux[0].qx + 10.0) < 1e-9, "parser solve flux");
  writeVtk("selftest_mesh.vtk", m, r);
  std::remove(path.c_str());
  std::remove("selftest_mesh.vtk");
}

template <typename F>
static void expectThrow(F&& f, const std::string& what) {
  try {
    f();
    CHECK(false, "expected Error: " + what);
  } catch (const Error&) {
  }
}

static void testErrors() {
  {
    Mesh m = unitSquare();
    m.nodes.push_back({1, 0.5, 0.5});  // duplicate id
    expectThrow([&] { validateMesh(m); }, "duplicate node id");
  }
  {
    Mesh m = unitSquare();
    m.triangles[0].nodes[2] = 99;  // unknown node
    expectThrow([&] { validateMesh(m); }, "unknown node reference");
  }
  {
    Mesh m = unitSquare();
    m.triangles.push_back({103, 10, {1, 2, 5}});
    m.nodes.push_back({5, 2, 0});  // collinear with 1-2
    expectThrow([&] { validateMesh(m); }, "degenerate triangle");
  }
  {
    Mesh m = unitSquare();
    m.nodes.push_back({5, 2, 1});
    m.triangles.push_back({103, 10, {1, 3, 5}});  // edge 1-3 now shared by 3 triangles
    expectThrow([&] { validateMesh(m); }, "non-manifold edge");
  }
  {
    Mesh m = unitSquare();
    m.lines.push_back({205, 5, {1, 3}});  // diagonal is interior
    expectThrow([&] { validateMesh(m); }, "boundary line not on outer boundary");
  }
  {
    Mesh m = unitSquare();
    m.triangles.push_back({103, 10, {1, 2, 3}});  // duplicate element id? no, new id but dup element ids below
    m.triangles.back().id = 101;
    expectThrow([&] { validateMesh(m); }, "duplicate element id");
  }
  {
    SolveOptions opts;  // no material for tag 10
    opts.boundaries[4] = DirichletBC{0.0};
    expectThrow([&] { solve(unitSquare(), opts); }, "missing material");
  }
  {
    SolveOptions opts;
    opts.materials[10] = Material{-1.0, 0.0};
    opts.boundaries[4] = DirichletBC{0.0};
    expectThrow([&] { solve(unitSquare(), opts); }, "negative conductivity");
  }
  {
    SolveOptions opts;
    opts.materials[10] = Material{1.0, 0.0};
    opts.boundaries[1] = NeumannBC{5.0};  // only flux -> singular
    expectThrow([&] { solve(unitSquare(), opts); }, "temperature not unique");
  }
  {
    SolveOptions opts;
    opts.materials[10] = Material{1.0, 0.0};
    opts.boundaries[1] = DirichletBC{0.0};
    opts.boundaries[4] = DirichletBC{3.0};  // share node 1 -> conflict
    expectThrow([&] { solve(unitSquare(), opts); }, "conflicting Dirichlet");
  }
  {
    SolveOptions opts;
    opts.materials[10] = Material{1.0, 0.0};
    opts.boundaries[4] = ConvectionBC{0.0, 20.0};  // h must be positive
    expectThrow([&] { solve(unitSquare(), opts); }, "non-positive h");
  }
  {
    // Two disconnected components, one without any anchoring BC.
    Mesh m = unitSquare();
    const long off = 10;
    for (auto& nd : std::vector<Node>{{1, 0, 0}, {2, 1, 0}, {3, 1, 1}, {4, 0, 1}})
      m.nodes.push_back({nd.id + off, nd.x + 3.0, nd.y});
    m.triangles.push_back({111, 10, {11, 12, 13}});
    m.triangles.push_back({112, 10, {11, 13, 14}});
    SolveOptions opts;
    opts.materials[10] = Material{1.0, 0.0};
    opts.boundaries[4] = DirichletBC{0.0};  // anchors only the first component
    expectThrow([&] { solve(m, opts); }, "unanchored component");
  }
}

int main() {
  testDirichletLinear();
  testConvectionSource();
  testNeumannDirichlet();
  testParser();
  testErrors();
  if (g_failures == 0) {
    std::cout << "all self-tests passed\n";
    return 0;
  }
  std::cerr << g_failures << " check(s) failed\n";
  return 1;
}
