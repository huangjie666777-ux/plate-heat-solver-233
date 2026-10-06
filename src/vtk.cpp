#include "plate_heat233/plate_heat233.hpp"

#include <fstream>
#include <map>

namespace plate_heat233 {

void writeVtk(const std::string& path, const Mesh& mesh, const SolveResult& result) {
  std::ofstream out(path);
  if (!out) throw Error("cannot open VTK output file: " + path);
  out.precision(17);

  const std::size_t nNodes = mesh.nodes.size();
  const std::size_t nTri = mesh.triangles.size();
  const std::size_t nLines = mesh.lines.size();

  std::map<long, int> indexOf;
  for (std::size_t i = 0; i < nNodes; ++i) indexOf[mesh.nodes[i].id] = static_cast<int>(i);

  out << "# vtk DataFile Version 3.0\n";
  out << "plate_heat233 steady-state heat conduction\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << nNodes << " double\n";
  for (const auto& n : mesh.nodes) out << n.x << " " << n.y << " 0\n";

  out << "CELLS " << (nTri + nLines) << " " << (nTri * 4 + nLines * 3) << "\n";
  for (const auto& t : mesh.triangles)
    out << "3 " << indexOf[t.nodes[0]] << " " << indexOf[t.nodes[1]] << " " << indexOf[t.nodes[2]] << "\n";
  for (const auto& ln : mesh.lines)
    out << "2 " << indexOf[ln.nodes[0]] << " " << indexOf[ln.nodes[1]] << "\n";

  out << "CELL_TYPES " << (nTri + nLines) << "\n";
  for (std::size_t i = 0; i < nTri; ++i) out << "5\n";   // VTK_TRIANGLE
  for (std::size_t i = 0; i < nLines; ++i) out << "3\n";  // VTK_LINE

  out << "POINT_DATA " << nNodes << "\n";
  out << "SCALARS TEMPERATURE double 1\n";
  out << "LOOKUP_TABLE default\n";
  for (const auto& n : mesh.nodes) {
    auto it = result.temperature.find(n.id);
    if (it == result.temperature.end())
      throw Error("VTK export: missing temperature for node " + std::to_string(n.id));
    out << it->second << "\n";
  }

  out << "CELL_DATA " << (nTri + nLines) << "\n";
  out << "VECTORS HEAT_FLUX double\n";
  std::map<long, TriangleHeatFlux> fluxById;
  for (const auto& f : result.heatFlux) fluxById[f.elementId] = f;
  for (const auto& t : mesh.triangles) {
    auto it = fluxById.find(t.id);
    if (it == fluxById.end())
      throw Error("VTK export: missing heat flux for triangle " + std::to_string(t.id));
    out << it->second.qx << " " << it->second.qy << " 0\n";
  }
  for (std::size_t i = 0; i < nLines; ++i) out << "0 0 0\n";

  out << "SCALARS REGION_TAG int 1\n";
  out << "LOOKUP_TABLE default\n";
  for (const auto& t : mesh.triangles) out << t.tag << "\n";
  for (const auto& ln : mesh.lines) out << ln.tag << "\n";

  if (!out) throw Error("failed writing VTK file: " + path);
}

}  // namespace plate_heat233
