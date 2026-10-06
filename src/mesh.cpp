#include "plate_heat233/plate_heat233.hpp"

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace plate_heat233 {

namespace {

constexpr std::size_t kMaxNodes = 5000;
constexpr std::size_t kMaxTriangles = 10000;

std::string stripCR(std::string s) {
  if (!s.empty() && s.back() == '\r') s.pop_back();
  return s;
}

std::map<long, int> nodeIndexMap(const Mesh& mesh) {
  std::map<long, int> idx;
  for (std::size_t i = 0; i < mesh.nodes.size(); ++i) idx[mesh.nodes[i].id] = static_cast<int>(i);
  return idx;
}

}  // namespace

void validateMesh(const Mesh& mesh) {
  if (mesh.nodes.empty()) throw Error("mesh has no nodes");
  if (mesh.triangles.empty()) throw Error("mesh has no triangles");
  if (mesh.nodes.size() > kMaxNodes)
    throw Error("node count " + std::to_string(mesh.nodes.size()) + " exceeds limit 5000");
  if (mesh.triangles.size() > kMaxTriangles)
    throw Error("triangle count " + std::to_string(mesh.triangles.size()) + " exceeds limit 10000");

  std::set<long> nodeIds;
  for (const auto& n : mesh.nodes) {
    if (!std::isfinite(n.x) || !std::isfinite(n.y))
      throw Error("node " + std::to_string(n.id) + " has non-finite coordinates");
    if (!nodeIds.insert(n.id).second)
      throw Error("duplicate node id " + std::to_string(n.id));
  }

  std::set<long> elemIds;
  auto checkElemId = [&](long id) {
    if (!elemIds.insert(id).second)
      throw Error("duplicate element id " + std::to_string(id));
  };
  auto checkNodeRef = [&](long nid, long elemId) {
    if (!nodeIds.count(nid))
      throw Error("element " + std::to_string(elemId) + " references unknown node " + std::to_string(nid));
  };

  for (const auto& t : mesh.triangles) {
    checkElemId(t.id);
    for (long nid : t.nodes) checkNodeRef(nid, t.id);
    if (t.nodes[0] == t.nodes[1] || t.nodes[1] == t.nodes[2] || t.nodes[2] == t.nodes[0])
      throw Error("triangle " + std::to_string(t.id) + " repeats a node");
  }
  for (const auto& ln : mesh.lines) {
    checkElemId(ln.id);
    for (long nid : ln.nodes) checkNodeRef(nid, ln.id);
    if (ln.nodes[0] == ln.nodes[1])
      throw Error("boundary line " + std::to_string(ln.id) + " repeats a node");
  }

  const auto idx = nodeIndexMap(mesh);

  // Degenerate triangle check (scale-aware).
  for (const auto& t : mesh.triangles) {
    const Node& a = mesh.nodes[idx.at(t.nodes[0])];
    const Node& b = mesh.nodes[idx.at(t.nodes[1])];
    const Node& c = mesh.nodes[idx.at(t.nodes[2])];
    const double cross = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
    const double scale = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) +
                         (c.x - b.x) * (c.x - b.x) + (c.y - b.y) * (c.y - b.y) +
                         (a.x - c.x) * (a.x - c.x) + (a.y - c.y) * (a.y - c.y);
    if (std::abs(cross) <= 1e-12 * scale)
      throw Error("degenerate triangle " + std::to_string(t.id));
  }

  // Edge incidence: non-manifold detection and outer-boundary check.
  std::map<std::pair<int, int>, int> edgeCount;
  auto edgeKey = [](int u, int v) { return u < v ? std::make_pair(u, v) : std::make_pair(v, u); };
  for (const auto& t : mesh.triangles) {
    const int v[3] = {idx.at(t.nodes[0]), idx.at(t.nodes[1]), idx.at(t.nodes[2])};
    for (int e = 0; e < 3; ++e) {
      const int c = ++edgeCount[edgeKey(v[e], v[(e + 1) % 3])];
      if (c > 2) throw Error("non-manifold edge shared by more than two triangles");
    }
  }
  std::set<std::pair<int, int>> lineEdges;
  for (const auto& ln : mesh.lines) {
    const auto key = edgeKey(idx.at(ln.nodes[0]), idx.at(ln.nodes[1]));
    auto it = edgeCount.find(key);
    if (it == edgeCount.end() || it->second != 1)
      throw Error("boundary line " + std::to_string(ln.id) + " is not on the outer boundary");
    if (!lineEdges.insert(key).second)
      throw Error("duplicate boundary line on the same edge (element " + std::to_string(ln.id) + ")");
  }
}

Mesh loadGmsh(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw Error("cannot open mesh file: " + path);

  Mesh mesh;
  bool formatSeen = false, nodesSeen = false, elementsSeen = false;
  std::string line;
  long lineNo = 0;
  auto fail = [&](const std::string& msg) -> Error {
    return Error(path + ":" + std::to_string(lineNo) + ": " + msg);
  };

  while (std::getline(in, line)) {
    ++lineNo;
    line = stripCR(line);
    if (line.empty() || line[0] != '$') continue;

    if (line == "$MeshFormat") {
      std::string fmt;
      if (!std::getline(in, fmt)) throw fail("unexpected EOF in $MeshFormat");
      ++lineNo;
      std::istringstream fs(stripCR(fmt));
      std::string version;
      int fileType = -1;
      int dataSize = -1;
      if (!(fs >> version >> fileType >> dataSize)) throw fail("malformed $MeshFormat");
      if (version.rfind("2.2", 0) != 0) throw fail("unsupported Gmsh version " + version + " (need 2.2)");
      if (fileType != 0) throw fail("binary Gmsh files are not supported");
      formatSeen = true;
    } else if (line == "$Nodes") {
      std::string cnt;
      if (!std::getline(in, cnt)) throw fail("unexpected EOF in $Nodes");
      ++lineNo;
      long count = 0;
      { std::istringstream cs(stripCR(cnt)); if (!(cs >> count) || count < 0) throw fail("malformed node count"); }
      for (long i = 0; i < count; ++i) {
        std::string rec;
        if (!std::getline(in, rec)) throw fail("unexpected EOF in $Nodes");
        ++lineNo;
        std::istringstream rs(stripCR(rec));
        Node n;
        double z = 0.0;
        if (!(rs >> n.id >> n.x >> n.y >> z)) throw fail("malformed node record");
        if (std::abs(z) > 1e-12) throw fail("node " + std::to_string(n.id) + " has nonzero z");
        mesh.nodes.push_back(n);
      }
      nodesSeen = true;
    } else if (line == "$Elements") {
      std::string cnt;
      if (!std::getline(in, cnt)) throw fail("unexpected EOF in $Elements");
      ++lineNo;
      long count = 0;
      { std::istringstream cs(stripCR(cnt)); if (!(cs >> count) || count < 0) throw fail("malformed element count"); }
      for (long i = 0; i < count; ++i) {
        std::string rec;
        if (!std::getline(in, rec)) throw fail("unexpected EOF in $Elements");
        ++lineNo;
        std::istringstream rs(stripCR(rec));
        long id = 0;
        int type = 0;
        int ntags = 0;
        if (!(rs >> id >> type >> ntags)) throw fail("malformed element header");
        if (ntags < 1) throw fail("element " + std::to_string(id) + " has no physical tag");
        int physTag = 0;
        for (int t = 0; t < ntags; ++t) {
          int tag = 0;
          if (!(rs >> tag)) throw fail("malformed element tags");
          if (t == 0) physTag = tag;
        }
        if (type == 1) {
          BoundaryLine ln;
          ln.id = id;
          ln.tag = physTag;
          if (!(rs >> ln.nodes[0] >> ln.nodes[1])) throw fail("malformed line element");
          mesh.lines.push_back(ln);
        } else if (type == 2) {
          Triangle tri;
          tri.id = id;
          tri.tag = physTag;
          if (!(rs >> tri.nodes[0] >> tri.nodes[1] >> tri.nodes[2])) throw fail("malformed triangle element");
          mesh.triangles.push_back(tri);
        } else {
          throw fail("unsupported element type " + std::to_string(type) +
                     " (only 1=line and 2=triangle)");
        }
      }
      elementsSeen = true;
    } else if (line.rfind("$End", 0) == 0) {
      continue;  // closing tag of a known section
    } else {
      // Skip unknown section (e.g. $PhysicalNames, $NodeData).
      const std::string endTag = "$End" + line.substr(1);
      while (std::getline(in, line)) {
        ++lineNo;
        if (stripCR(line) == endTag) break;
      }
    }
  }

  if (!formatSeen) throw Error(path + ": missing $MeshFormat section");
  if (!nodesSeen) throw Error(path + ": missing $Nodes section");
  if (!elementsSeen) throw Error(path + ": missing $Elements section");

  validateMesh(mesh);
  return mesh;
}

}  // namespace plate_heat233
