#include "plate_heat233/solver.h"

#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <set>
#include <vector>

namespace plate_heat233 {
namespace {

Result fail(std::string msg) {
    Result r;
    r.success = false;
    r.error = std::move(msg);
    return r;
}

bool finite(double v) { return std::isfinite(v); }

}  // namespace

Result solve(const Problem& problem) {
    const Mesh& mesh = problem.mesh;
    const std::size_t n = mesh.nodes().size();
    const auto& tris = mesh.triangles();
    const auto& lines = mesh.boundary_lines();

    // ---- Validate materials and boundary conditions ----
    for (const auto& t : tris) {
        auto it = problem.materials.find(t.region_tag);
        if (it == problem.materials.end())
            return fail("missing material for region tag " + std::to_string(t.region_tag));
        const Material& m = it->second;
        if (!(m.k > 0.0) || !finite(m.k))
            return fail("region tag " + std::to_string(t.region_tag) +
                        " has non-positive or non-finite conductivity");
        if (!finite(m.Q))
            return fail("region tag " + std::to_string(t.region_tag) +
                        " has non-finite heat source");
    }
    for (const auto& kv : problem.boundaries) {
        const BoundaryCondition& bc = kv.second;
        if (const auto* c = std::get_if<ConvectionBC>(&bc)) {
            if (!(c->h > 0.0) || !finite(c->h) || !finite(c->Ta))
                return fail("boundary tag " + std::to_string(kv.first) +
                            " has invalid convection data");
        } else if (const auto* d = std::get_if<DirichletBC>(&bc)) {
            if (!finite(d->T))
                return fail("boundary tag " + std::to_string(kv.first) +
                            " has non-finite prescribed temperature");
        } else {
            const auto& f = std::get<FluxBC>(bc);
            if (!finite(f.q))
                return fail("boundary tag " + std::to_string(kv.first) +
                            " has non-finite heat flux");
        }
    }

    // ---- Geometry helpers ----
    auto coords = [&](std::size_t i) {
        const Node& nd = mesh.nodes()[i];
        return std::array<double, 2>{nd.x, nd.y};
    };

    // ---- Assemble global system K T = b (unit thickness) ----
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(tris.size() * 9 + lines.size() * 4);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

    std::vector<std::array<int, 3>> tri_nodes(tris.size());
    std::vector<double> tri_area(tris.size());
    std::vector<std::array<double, 3>> tri_b(tris.size()), tri_c(tris.size());

    for (std::size_t e = 0; e < tris.size(); ++e) {
        const Triangle& t = tris[e];
        const Material& mat = problem.materials.at(t.region_tag);
        std::array<int, 3> idx;
        for (int j = 0; j < 3; ++j) idx[j] = static_cast<int>(mesh.index_of(t.nodes[j]));
        auto p0 = coords(idx[0]), p1 = coords(idx[1]), p2 = coords(idx[2]);
        double det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        double area = 0.5 * std::abs(det);
        double b[3] = {p1[1] - p2[1], p2[1] - p0[1], p0[1] - p1[1]};
        double c[3] = {p2[0] - p1[0], p0[0] - p2[0], p1[0] - p0[0]};
        double scale = mat.k / (4.0 * area);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                triplets.emplace_back(idx[i], idx[j],
                                      scale * (b[i] * b[j] + c[i] * c[j]));
        double qe = mat.Q * area / 3.0;
        for (int i = 0; i < 3; ++i) rhs[idx[i]] += qe;
        tri_nodes[e] = idx;
        tri_area[e] = area;
        for (int i = 0; i < 3; ++i) {
            tri_b[e][i] = b[i] / (2.0 * area);
            tri_c[e][i] = c[i] / (2.0 * area);
        }
    }

    // Boundary edges: natural conditions contribute to K / b.
    std::vector<char> node_dirichlet(n, 0);
    std::vector<double> node_dirichlet_value(n, 0.0);
    std::vector<char> node_has_convection(n, 0);

    for (const auto& l : lines) {
        auto bcit = problem.boundaries.find(l.boundary_tag);
        if (bcit == problem.boundaries.end()) continue;  // adiabatic by default
        std::size_t i0 = mesh.index_of(l.nodes[0]);
        std::size_t i1 = mesh.index_of(l.nodes[1]);
        auto p0 = coords(i0), p1 = coords(i1);
        double len = std::hypot(p1[0] - p0[0], p1[1] - p0[1]);
        const BoundaryCondition& bc = bcit->second;
        if (const auto* d = std::get_if<DirichletBC>(&bc)) {
            for (std::size_t i : {i0, i1}) {
                if (node_dirichlet[i] && node_dirichlet_value[i] != d->T)
                    return fail("conflicting prescribed temperatures at node " +
                                std::to_string(mesh.nodes()[i].id));
                node_dirichlet[i] = 1;
                node_dirichlet_value[i] = d->T;
            }
        } else if (const auto* f = std::get_if<FluxBC>(&bc)) {
            // Outward flux positive: -integral(q v) on the right-hand side.
            rhs[i0] -= 0.5 * f->q * len;
            rhs[i1] -= 0.5 * f->q * len;
        } else {
            const auto* cv = std::get_if<ConvectionBC>(&bc);
            double hL = cv->h * len;
            triplets.emplace_back(i0, i0, hL / 3.0);
            triplets.emplace_back(i1, i1, hL / 3.0);
            triplets.emplace_back(i0, i1, hL / 6.0);
            triplets.emplace_back(i1, i0, hL / 6.0);
            rhs[i0] += 0.5 * hL * cv->Ta;
            rhs[i1] += 0.5 * hL * cv->Ta;
            node_has_convection[i0] = 1;
            node_has_convection[i1] = 1;
        }
    }

    // ---- Connected components of the triangle mesh ----
    // Each component needs at least one Dirichlet node or convection edge.
    {
        std::vector<std::vector<int>> node_tris(n);
        for (std::size_t e = 0; e < tris.size(); ++e)
            for (int j = 0; j < 3; ++j) node_tris[tri_nodes[e][j]].push_back(static_cast<int>(e));
        std::vector<int> comp(tris.size(), -1);
        int ncomp = 0;
        for (std::size_t e = 0; e < tris.size(); ++e) {
            if (comp[e] != -1) continue;
            std::vector<char> fixed(n, 0);
            std::queue<std::size_t> q;
            q.push(e);
            comp[e] = ncomp;
            bool anchored = false;
            while (!q.empty()) {
                std::size_t cur = q.front();
                q.pop();
                for (int j = 0; j < 3; ++j) {
                    int nd = tri_nodes[cur][j];
                    if (node_dirichlet[nd] || node_has_convection[nd]) anchored = true;
                    if (fixed[nd]) continue;
                    fixed[nd] = 1;
                    for (int nb : node_tris[nd])
                        if (comp[nb] == -1) {
                            comp[nb] = ncomp;
                            q.push(nb);
                        }
                }
            }
            if (!anchored)
                return fail("connected component " + std::to_string(ncomp) +
                            " has no Dirichlet or convection boundary: temperature is not unique");
            ++ncomp;
        }
    }

    // ---- Build the original sparse system (for residual / reactions) ----
    Eigen::SparseMatrix<double> K(n, n);
    K.setFromTriplets(triplets.begin(), triplets.end());
    const Eigen::VectorXd rhs_original = rhs;  // kept for residual / reactions

    // ---- Dirichlet elimination with right-hand-side adjustment ----
    std::vector<int> free_of(n, -1);
    std::vector<int> free_nodes;
    for (std::size_t i = 0; i < n; ++i)
        if (!node_dirichlet[i]) {
            free_of[i] = static_cast<int>(free_nodes.size());
            free_nodes.push_back(static_cast<int>(i));
        }

    const std::size_t nf = free_nodes.size();
    Eigen::VectorXd T = Eigen::VectorXd::Zero(n);
    for (std::size_t i = 0; i < n; ++i)
        if (node_dirichlet[i]) T[i] = node_dirichlet_value[i];

    if (nf > 0) {
        std::vector<Triplet> ftrip;
        ftrip.reserve(triplets.size());
        for (int k = 0; k < K.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
                std::size_t row = it.row(), col = it.col();
                if (node_dirichlet[row]) continue;
                if (node_dirichlet[col])
                    rhs[row] -= it.value() * node_dirichlet_value[col];
                else
                    ftrip.emplace_back(free_of[row], free_of[col], it.value());
            }
        Eigen::SparseMatrix<double> A(nf, nf);
        A.setFromTriplets(ftrip.begin(), ftrip.end());
        Eigen::VectorXd b(nf);
        for (std::size_t i = 0; i < nf; ++i) b[i] = rhs[free_nodes[i]];

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success)
            return fail("factorization of the reduced system failed");
        Eigen::VectorXd x = solver.solve(b);
        if (solver.info() != Eigen::Success || !x.allFinite())
            return fail("linear solve failed");
        for (std::size_t i = 0; i < nf; ++i) T[free_nodes[i]] = x[i];
    }

    // ---- Residual on free nodes of the original system ----
    Eigen::VectorXd residual = K * T - rhs_original;
    double max_res = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        if (!node_dirichlet[i]) max_res = std::max(max_res, std::abs(residual[i]));
    if (!(max_res <= problem.tolerance))
        return fail("residual tolerance not met: max free-node residual " +
                    std::to_string(max_res) + " > " + std::to_string(problem.tolerance));

    // ---- Post-processing ----
    Result result;
    result.success = true;
    result.max_residual = max_res;
    for (std::size_t i = 0; i < n; ++i) result.temperature[mesh.nodes()[i].id] = T[i];

    // Per-triangle heat flux -k grad T.
    result.heat_flux.resize(tris.size());
    for (std::size_t e = 0; e < tris.size(); ++e) {
        const Material& mat = problem.materials.at(tris[e].region_tag);
        double gx = 0.0, gy = 0.0;
        for (int i = 0; i < 3; ++i) {
            double Ti = T[tri_nodes[e][i]];
            gx += tri_b[e][i] * Ti;
            gy += tri_c[e][i] * Ti;
        }
        result.heat_flux[e] = {-mat.k * gx, -mat.k * gy};
        result.total_source_power += mat.Q * tri_area[e];
    }

    // Boundary powers (outward positive). Dirichlet power comes from the
    // reactions of the original (non-eliminated) system, summed once per
    // node so corner nodes shared by two lines are not double counted.
    std::map<int, std::set<std::size_t>> dirichlet_nodes;
    for (const auto& l : lines) {
        std::size_t i0 = mesh.index_of(l.nodes[0]);
        std::size_t i1 = mesh.index_of(l.nodes[1]);
        auto p0 = coords(i0), p1 = coords(i1);
        double len = std::hypot(p1[0] - p0[0], p1[1] - p0[1]);
        auto bcit = problem.boundaries.find(l.boundary_tag);
        if (bcit == problem.boundaries.end()) {
            result.boundary_power[l.boundary_tag] += 0.0;  // adiabatic
        } else if (std::holds_alternative<DirichletBC>(bcit->second)) {
            dirichlet_nodes[l.boundary_tag].insert(i0);
            dirichlet_nodes[l.boundary_tag].insert(i1);
        } else if (const auto* f = std::get_if<FluxBC>(&bcit->second)) {
            result.boundary_power[l.boundary_tag] += f->q * len;
        } else {
            const auto* cv = std::get_if<ConvectionBC>(&bcit->second);
            double Tm = 0.5 * (T[i0] + T[i1]);
            result.boundary_power[l.boundary_tag] += cv->h * (Tm - cv->Ta) * len;
        }
    }
    for (const auto& kv : dirichlet_nodes) {
        // Reaction r = K T - b at a constrained node is the heat the
        // boundary supplies to the domain; outward power is -r.
        double power = 0.0;
        for (std::size_t i : kv.second) power -= residual[i];
        result.boundary_power[kv.first] += power;
    }

    double outward = 0.0;
    for (const auto& kv : result.boundary_power) outward += kv.second;
    result.power_imbalance = std::abs(result.total_source_power - outward);
    return result;
}

}  // namespace plate_heat233
