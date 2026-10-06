#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace plate_heat233 {

struct MeshError : std::runtime_error {
    explicit MeshError(const std::string& msg) : std::runtime_error(msg) {}
};

struct Node {
    int id = 0;
    double x = 0.0;
    double y = 0.0;
};

struct Triangle {
    int id = 0;
    int region_tag = 0;          // first physical tag
    std::array<int, 3> nodes{};  // node ids
};

struct BoundaryLine {
    int id = 0;
    int boundary_tag = 0;        // first physical tag
    std::array<int, 2> nodes{};  // node ids
};

class Mesh {
public:
    static constexpr std::size_t kMaxNodes = 5000;
    static constexpr std::size_t kMaxTriangles = 10000;

    // Parse a Gmsh 2.2 ASCII file. Throws MeshError on any inconsistency.
    static Mesh read_gmsh(const std::string& path);

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Triangle>& triangles() const { return triangles_; }
    const std::vector<BoundaryLine>& boundary_lines() const { return lines_; }

    // Dense index (0..n-1) of a node id; throws MeshError if unknown.
    std::size_t index_of(int node_id) const { return id_to_index_.at(node_id); }

private:
    std::vector<Node> nodes_;
    std::vector<Triangle> triangles_;
    std::vector<BoundaryLine> lines_;
    std::map<int, std::size_t> id_to_index_;

    void validate();
};

}  // namespace plate_heat233

