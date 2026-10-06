#include "plate_heat233/mesh.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace plate_heat233 {
namespace {

[[noreturn]] void fail(const std::string& msg) { throw MeshError(msg); }

std::string require_line(std::istream& in, const std::string& what) {
    std::string line;
    if (!std::getline(in, line)) fail("unexpected end of file while reading " + what);
    return line;
}

}  // namespace

Mesh Mesh::read_gmsh(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail("cannot open mesh file: " + path);

    Mesh mesh;
    std::set<int> node_ids, elem_ids;

    while (true) {
        std::string marker;
        if (!(in >> marker)) break;
        if (marker == "$MeshFormat") {
            double version;
            int file_type, data_size;
            if (!(in >> version >> file_type >> data_size)) fail("bad $MeshFormat section");
            if (version < 2.0 || version >= 3.0) fail("only Gmsh 2.x ASCII format is supported");
            if (file_type != 0) fail("only ASCII Gmsh files are supported");
            std::string end;
            if (!(in >> end) || end != "$EndMeshFormat") fail("missing $EndMeshFormat");
        } else if (marker == "$Nodes") {
            std::size_t count;
            if (!(in >> count)) fail("bad $Nodes header");
            if (count > kMaxNodes) fail("node count exceeds limit of 5000");
            for (std::size_t i = 0; i < count; ++i) {
                Node n;
                double z;
                if (!(in >> n.id >> n.x >> n.y >> z)) fail("bad node record");
                if (z != 0.0) fail("node " + std::to_string(n.id) + " has non-zero z coordinate");
                if (!node_ids.insert(n.id).second)
                    fail("duplicate node id " + std::to_string(n.id));
                mesh.id_to_index_[n.id] = mesh.nodes_.size();
                mesh.nodes_.push_back(n);
            }
            std::string end;
            if (!(in >> end) || end != "$EndNodes") fail("missing $EndNodes");
        } else if (marker == "$Elements") {
            std::size_t count;
            if (!(in >> count)) fail("bad $Elements header");
            std::size_t ntri = 0;
            for (std::size_t i = 0; i < count; ++i) {
                int id, type, ntags;
                if (!(in >> id >> type >> ntags)) fail("bad element header");
                std::vector<int> tags(ntags);
                for (int j = 0; j < ntags; ++j)
                    if (!(in >> tags[j])) fail("bad element tags");
                int nverts = (type == 1) ? 2 : (type == 2) ? 3 : 0;
                if (nverts == 0) {
                    fail("unsupported element type " + std::to_string(type) +
                         " (only 2-node lines and 3-node triangles)");
                }
                std::vector<int> conn(nverts);
                for (int j = 0; j < nverts; ++j)
                    if (!(in >> conn[j])) fail("bad element connectivity");
                if (!elem_ids.insert(id).second)
                    fail("duplicate element id " + std::to_string(id));
                int phys = ntags >= 1 ? tags[0] : 0;
                if (phys <= 0)
                    fail("element " + std::to_string(id) + " has no positive physical tag");
                if (type == 2) {
                    if (++ntri > kMaxTriangles) fail("triangle count exceeds limit of 10000");
                    mesh.triangles_.push_back(Triangle{id, phys, {conn[0], conn[1], conn[2]}});
                } else {
                    mesh.lines_.push_back(BoundaryLine{id, phys, {conn[0], conn[1]}});
                }
            }
            std::string end;
            if (!(in >> end) || end != "$EndElements") fail("missing $EndElements");
        } else if (marker == "$PhysicalNames") {
            std::string line;
            std::getline(in, line);  // rest of header line
            std::size_t count = 0;
            std::istringstream(line) >> count;
            for (std::size_t i = 0; i < count; ++i) require_line(in, "physical name");
            require_line(in, "$EndPhysicalNames");
        } else {
            // Skip unknown section generically: read lines until matching $End*.
            std::string end_marker = "$End" + marker.substr(1);
            std::string line;
            while (std::getline(in, line))
                if (line.find(end_marker) != std::string::npos) break;
        }
    }

    if (mesh.nodes_.empty()) fail("mesh contains no nodes");
    if (mesh.triangles_.empty()) fail("mesh contains no triangles");
    mesh.validate();
    return mesh;
}

void Mesh::validate() {
    // Unknown node references and degenerate triangles.
    std::map<std::pair<int, int>, int> edge_count;  // sorted dense-index edge -> count
    std::map<std::pair<int, int>, std::vector<int>> edge_tris;

    for (const auto& t : triangles_) {
        std::array<int, 3> idx;
        for (int j = 0; j < 3; ++j) {
            auto it = id_to_index_.find(t.nodes[j]);
            if (it == id_to_index_.end())
                fail("triangle " + std::to_string(t.id) + " references unknown node " +
                     std::to_string(t.nodes[j]));
            idx[j] = static_cast<int>(it->second);
        }
        if (idx[0] == idx[1] || idx[1] == idx[2] || idx[2] == idx[0])
            fail("degenerate triangle " + std::to_string(t.id) + " (repeated node)");
        const Node& a = nodes_[idx[0]];
        const Node& b = nodes_[idx[1]];
        const Node& c = nodes_[idx[2]];
        double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        if (cross == 0.0)
            fail("degenerate triangle " + std::to_string(t.id) + " (zero area)");
        for (int e = 0; e < 3; ++e) {
            int u = idx[e], v = idx[(e + 1) % 3];
            auto key = std::make_pair(std::min(u, v), std::max(u, v));
            edge_count[key]++;
            edge_tris[key].push_back(t.id);
        }
    }

    // Non-manifold edges.
    for (const auto& kv : edge_count) {
        if (kv.second > 2)
            fail("non-manifold edge between nodes " +
                 std::to_string(nodes_[kv.first.first].id) + " and " +
                 std::to_string(nodes_[kv.first.second].id));
    }

    // Boundary lines must reference known nodes and lie on the outer boundary
    // (an edge belonging to exactly one triangle).
    for (const auto& l : lines_) {
        auto it0 = id_to_index_.find(l.nodes[0]);
        auto it1 = id_to_index_.find(l.nodes[1]);
        if (it0 == id_to_index_.end() || it1 == id_to_index_.end())
            fail("boundary line " + std::to_string(l.id) + " references unknown node");
        if (l.nodes[0] == l.nodes[1])
            fail("degenerate boundary line " + std::to_string(l.id));
        auto key = std::make_pair(static_cast<int>(it0->second), static_cast<int>(it1->second)); if (key.first > key.second) std::swap(key.first, key.second);
        auto it = edge_count.find(key);
        if (it == edge_count.end() || it->second != 1)
            fail("boundary line " + std::to_string(l.id) +
                 " is not on the outer boundary of the mesh");
    }
}

}  // namespace plate_heat233
