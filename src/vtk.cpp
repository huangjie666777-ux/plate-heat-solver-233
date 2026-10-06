#include "plate_heat233/vtk.h"

#include <fstream>

namespace plate_heat233 {

bool write_vtk(const std::string& path, const Mesh& mesh, const Result& result,
               std::string& error) {
    std::ofstream out(path);
    if (!out) {
        error = "cannot open VTK file for writing: " + path;
        return false;
    }
    const auto& nodes = mesh.nodes();
    const auto& tris = mesh.triangles();
    const auto& lines = mesh.boundary_lines();

    out << "# vtk DataFile Version 3.0\nplate_heat233 temperature field\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << nodes.size() << " double\n";
    for (const auto& nd : nodes) out << nd.x << " " << nd.y << " 0\n";

    std::size_t ncells = tris.size() + lines.size();
    out << "CELLS " << ncells << " " << tris.size() * 4 + lines.size() * 3 << "\n";
    for (const auto& t : tris)
        out << "3 " << mesh.index_of(t.nodes[0]) << " " << mesh.index_of(t.nodes[1]) << " "
            << mesh.index_of(t.nodes[2]) << "\n";
    for (const auto& l : lines)
        out << "2 " << mesh.index_of(l.nodes[0]) << " " << mesh.index_of(l.nodes[1]) << "\n";

    out << "CELL_TYPES " << ncells << "\n";
    for (std::size_t i = 0; i < tris.size(); ++i) out << "5\n";   // VTK_TRIANGLE
    for (std::size_t i = 0; i < lines.size(); ++i) out << "3\n";  // VTK_LINE

    out << "POINT_DATA " << nodes.size() << "\n";
    out << "SCALARS temperature double 1\nLOOKUP_TABLE default\n";
    for (const auto& nd : nodes) {
        auto it = result.temperature.find(nd.id);
        if (it == result.temperature.end()) {
            error = "missing temperature for node " + std::to_string(nd.id);
            return false;
        }
        out << it->second << "\n";
    }

    out << "CELL_DATA " << ncells << "\n";
    out << "SCALARS region_tag int 1\nLOOKUP_TABLE default\n";
    for (const auto& t : tris) out << t.region_tag << "\n";
    for (const auto& l : lines) out << l.boundary_tag << "\n";
    out << "VECTORS heat_flux double\n";
    for (std::size_t e = 0; e < tris.size(); ++e) {
        if (e < result.heat_flux.size())
            out << result.heat_flux[e][0] << " " << result.heat_flux[e][1] << " 0\n";
        else
            out << "0 0 0\n";
    }
    for (std::size_t i = 0; i < lines.size(); ++i) out << "0 0 0\n";
    return true;
}

}  // namespace plate_heat233

