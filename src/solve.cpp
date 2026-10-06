#include "plate_heat233/plate_heat233.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <map>
#include <numeric>
#include <vector>

namespace plate_heat233 {

namespace {

struct DSU {
  std::vector<int> parent;
  explicit DSU(int n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
  int find(int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); }
  void unite(int a, int b) { parent[find(a)] = find(b); }
};

void validateOptions(const Mesh& mesh, const SolveOptions& opts) {
  if (!(opts.residualTolerance > 0.0) || !std::isfinite(opts.residualTolerance))
    throw Error("residual tolerance must be positive and finite");
  for (const auto& [tag, m] : opts.materials) {
    if (!(m.k > 0.0) || !std::isfinite(m.k))
      throw Error("region tag " + std::to_string(tag) + ": conductivity k must be positive and finite");
    if (!std::isfinite(m.Q))
      throw Error("region tag " + std::to_string(tag) + ": heat source Q must be finite");
  }
  for (const auto& [tag, bc] : opts.boundaries) {
    if (const auto* d = std::get_if<DirichletBC>(&bc)) {
      if (!std::isfinite(d->T))
        throw Error("boundary tag " + std::to_string(tag) + ": Dirichlet T must be finite");
    } else if (const auto* nbc = std::get_if<NeumannBC>(&bc)) {
      if (!std::isfinite(nbc->q))
        throw Error("boundary tag " + std::to_string(tag) + ": flux q must be finite");
    } else {
      const auto& c = std::get<ConvectionBC>(bc);
      if (!(c.h > 0.0) || !std::isfinite(c.h))
        throw Error("boundary tag " + std::to_string(tag) + ": convection h must be positive and finite");
      if (!std::isfinite(c.Ta))
        throw Error("boundary tag " + std::to_string(tag) + ": ambient temperature Ta must be finite");
    }
  }
  for (const auto& t : mesh.triangles) {
    if (!opts.materials.count(t.tag))
      throw Error("missing material for region tag " + std::to_string(t.tag));
  }
}

}  // namespace

SolveResult solve(const Mesh& mesh, const SolveOptions& opts) {
  validateMesh(mesh);
  validateOptions(mesh, opts);

  const int n = static_cast<int>(mesh.nodes.size());
  std::map<long, int> indexOf;
  for (int i = 0; i < n; ++i) indexOf[mesh.nodes[i].id] = i;

  auto findBc = [&](int tag) -> const BoundaryCondition* {
    auto it = opts.boundaries.find(tag);
    return it == opts.boundaries.end() ? nullptr : &it->second;
  };

  // Dirichlet values per node, with conflict detection.
  std::vector<char> isDir(n, 0);
  std::vector<double> dirT(n, 0.0);
  for (const auto& ln : mesh.lines) {
    const BoundaryCondition* bc = findBc(ln.tag);
    const auto* d = bc ? std::get_if<DirichletBC>(bc) : nullptr;
    if (!d) continue;
    for (long nid : ln.nodes) {
      const int i = indexOf[nid];
      if (isDir[i] && dirT[i] != d->T)
        throw Error("conflicting Dirichlet temperatures at node " + std::to_string(nid));
      isDir[i] = 1;
      dirT[i] = d->T;
    }
  }

  // Assembly of the full (unconstrained) system.
  using Triplet = Eigen::Triplet<double>;
  std::vector<Triplet> triplets;
  triplets.reserve(mesh.triangles.size() * 9 + mesh.lines.size() * 4);
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
  double sourcePower = 0.0;

  for (const auto& tri : mesh.triangles) {
    const Material& mat = opts.materials.at(tri.tag);
    int id[3];
    double x[3], y[3];
    for (int i = 0; i < 3; ++i) {
      id[i] = indexOf[tri.nodes[i]];
      x[i] = mesh.nodes[id[i]].x;
      y[i] = mesh.nodes[id[i]].y;
    }
    const double detJ = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
    const double area = 0.5 * std::abs(detJ);
    const double bvec[3] = {y[1] - y[2], y[2] - y[0], y[0] - y[1]};
    const double cvec[3] = {x[2] - x[1], x[0] - x[2], x[1] - x[0]};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double ke = mat.k * (bvec[i] * bvec[j] + cvec[i] * cvec[j]) / (4.0 * area);
        triplets.emplace_back(id[i], id[j], ke);
      }
      rhs[id[i]] += mat.Q * area / 3.0;
    }
    sourcePower += mat.Q * area;
  }

  double neumannPower = 0.0;
  for (const auto& ln : mesh.lines) {
    const BoundaryCondition* bc = findBc(ln.tag);
    if (!bc) continue;  // unspecified boundary: adiabatic
    const int i0 = indexOf[ln.nodes[0]];
    const int i1 = indexOf[ln.nodes[1]];
    const double len = std::hypot(mesh.nodes[i1].x - mesh.nodes[i0].x,
                                  mesh.nodes[i1].y - mesh.nodes[i0].y);
    if (const auto* nbc = std::get_if<NeumannBC>(bc)) {
      // Weak form: int k gradT.gradv = int Qv - oint q_out v, so an
      // outward-positive flux contributes -q to the load vector.
      const double f = -nbc->q * len / 2.0;
      rhs[i0] += f;
      rhs[i1] += f;
      neumannPower += nbc->q * len;
    } else if (const auto* c = std::get_if<ConvectionBC>(bc)) {
      const double s = c->h * len / 6.0;
      triplets.emplace_back(i0, i0, 2.0 * s);
      triplets.emplace_back(i1, i1, 2.0 * s);
      triplets.emplace_back(i0, i1, s);
      triplets.emplace_back(i1, i0, s);
      const double f = c->h * c->Ta * len / 2.0;
      rhs[i0] += f;
      rhs[i1] += f;
    }
  }

  // Connectivity: components of triangles joined through shared edges.
  const int nTri = static_cast<int>(mesh.triangles.size());
  DSU dsu(nTri);
  std::map<std::pair<int, int>, int> edgeToTri;
  std::vector<char> nodeInTri(n, 0);
  for (int t = 0; t < nTri; ++t) {
    int v[3];
    for (int i = 0; i < 3; ++i) {
      v[i] = indexOf[mesh.triangles[t].nodes[i]];
      nodeInTri[v[i]] = 1;
    }
    for (int e = 0; e < 3; ++e) {
      int a = v[e], b = v[(e + 1) % 3];
      if (a > b) std::swap(a, b);
      auto [it, inserted] = edgeToTri.emplace(std::make_pair(a, b), t);
      if (!inserted) dsu.unite(it->second, t);
    }
  }
  std::map<int, char> compHasDir, compHasConv;
  for (const auto& ln : mesh.lines) {
    const BoundaryCondition* bc = findBc(ln.tag);
    if (!bc) continue;
    int a = indexOf[ln.nodes[0]], b = indexOf[ln.nodes[1]];
    if (a > b) std::swap(a, b);
    const int root = dsu.find(edgeToTri.at({a, b}));
    if (std::holds_alternative<DirichletBC>(*bc)) compHasDir[root] = 1;
    if (std::holds_alternative<ConvectionBC>(*bc)) compHasConv[root] = 1;
  }
  for (int t = 0; t < nTri; ++t) {
    const int root = dsu.find(t);
    if (!compHasDir[root] && !compHasConv[root])
      throw Error("temperature not unique: connected component containing triangle " +
                  std::to_string(mesh.triangles[t].id) +
                  " has no Dirichlet or convection boundary");
  }
  for (int i = 0; i < n; ++i) {
    if (!nodeInTri[i] && !isDir[i])
      throw Error("node " + std::to_string(mesh.nodes[i].id) +
                  " belongs to no triangle and has no Dirichlet condition");
  }

  // Symmetric elimination of Dirichlet dofs: A_ff x_f = b_f - A_fc T_c.
  std::vector<int> freeOf(n, -1);
  int nFree = 0;
  for (int i = 0; i < n; ++i)
    if (!isDir[i]) freeOf[i] = nFree++;

  std::vector<Triplet> freeTriplets;
  freeTriplets.reserve(triplets.size());
  Eigen::VectorXd bf = Eigen::VectorXd::Zero(nFree);
  for (int i = 0; i < n; ++i)
    if (freeOf[i] >= 0) bf[freeOf[i]] = rhs[i];
  for (const auto& tr : triplets) {
    const int r = tr.row(), c = tr.col();
    if (freeOf[r] >= 0 && freeOf[c] >= 0)
      freeTriplets.emplace_back(freeOf[r], freeOf[c], tr.value());
    else if (freeOf[r] >= 0)
      bf[freeOf[r]] -= tr.value() * dirT[c];  // symmetric twin covers (c, r)
  }

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i)
    if (isDir[i]) x[i] = dirT[i];

  double maxResidual = 0.0;
  if (nFree > 0) {
    Eigen::SparseMatrix<double> A(nFree, nFree);
    A.setFromTriplets(freeTriplets.begin(), freeTriplets.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    if (solver.info() != Eigen::Success)
      throw Error("sparse factorization failed (singular or indefinite system)");
    Eigen::VectorXd xf = solver.solve(bf);
    if (solver.info() != Eigen::Success || !xf.allFinite())
      throw Error("linear solve failed");
    const Eigen::VectorXd res = bf - A * xf;
    maxResidual = res.cwiseAbs().maxCoeff();
    if (!(maxResidual <= opts.residualTolerance))
      throw Error("max free-node residual " + std::to_string(maxResidual) +
                  " exceeds tolerance " + std::to_string(opts.residualTolerance));
    for (int i = 0; i < n; ++i)
      if (freeOf[i] >= 0) x[i] = xf[freeOf[i]];
  }

  // Reactions from the original unconstrained system: R = A_full x - b.
  Eigen::SparseMatrix<double> AFull(n, n);
  AFull.setFromTriplets(triplets.begin(), triplets.end());
  const Eigen::VectorXd reaction = AFull * x - rhs;
  double dirichletPower = 0.0;  // outward power through Dirichlet boundaries
  for (int i = 0; i < n; ++i)
    if (isDir[i]) dirichletPower -= reaction[i];

  double convectionPower = 0.0;
  for (const auto& ln : mesh.lines) {
    const BoundaryCondition* bc = findBc(ln.tag);
    const auto* c = bc ? std::get_if<ConvectionBC>(bc) : nullptr;
    if (!c) continue;
    const int i0 = indexOf[ln.nodes[0]];
    const int i1 = indexOf[ln.nodes[1]];
    const double len = std::hypot(mesh.nodes[i1].x - mesh.nodes[i0].x,
                                  mesh.nodes[i1].y - mesh.nodes[i0].y);
    convectionPower += c->h * len * (0.5 * (x[i0] + x[i1]) - c->Ta);
  }

  SolveResult result;
  for (int i = 0; i < n; ++i) result.temperature[mesh.nodes[i].id] = x[i];

  result.heatFlux.reserve(mesh.triangles.size());
  for (const auto& tri : mesh.triangles) {
    const Material& mat = opts.materials.at(tri.tag);
    double xv[3], yv[3], tv[3];
    for (int i = 0; i < 3; ++i) {
      const int idx = indexOf[tri.nodes[i]];
      xv[i] = mesh.nodes[idx].x;
      yv[i] = mesh.nodes[idx].y;
      tv[i] = x[idx];
    }
    const double detJ = (xv[1] - xv[0]) * (yv[2] - yv[0]) - (xv[2] - xv[0]) * (yv[1] - yv[0]);
    const double bvec[3] = {yv[1] - yv[2], yv[2] - yv[0], yv[0] - yv[1]};
    const double cvec[3] = {xv[2] - xv[1], xv[0] - xv[2], xv[1] - xv[0]};
    double dTx = 0.0, dTy = 0.0;
    for (int i = 0; i < 3; ++i) {
      dTx += tv[i] * bvec[i];
      dTy += tv[i] * cvec[i];
    }
    dTx /= detJ;
    dTy /= detJ;
    result.heatFlux.push_back({tri.id, tri.tag, -mat.k * dTx, -mat.k * dTy});
  }

  result.maxResidual = maxResidual;
  result.power.source = sourcePower;
  result.power.dirichletOut = dirichletPower;
  result.power.neumannOut = neumannPower;
  result.power.convectionOut = convectionPower;
  result.power.imbalance = sourcePower - (dirichletPower + neumannPower + convectionPower);
  return result;
}

}  // namespace plate_heat233
