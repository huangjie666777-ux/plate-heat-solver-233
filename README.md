# plate_heat233

Steady-state 2-D heat conduction on a plate (unit thickness) with zoned
materials and thermal boundary conditions, solved with P1 triangular
finite elements. Pure C++17 library, no UI, no HTTP.

    -div(k grad T) = Q

Dependencies: Eigen 3.4.0 (vendored in third_party/eigen3), GCC 11.4,
GNU Make 4.3. Eigen include flag: -Ithird_party/eigen3

## Layout

- include/plate_heat233/mesh.h, src/mesh.cpp — Gmsh 2.2 ASCII mesh reader
  and topology validation
- include/plate_heat233/problem.h — materials, boundary conditions, Problem
- include/plate_heat233/result.h — solver output (temperatures, heat flux,
  residual, power balance)
- include/plate_heat233/solver.h, src/solver.cpp — assembly, Dirichlet
  elimination, sparse solve, post-processing
- include/plate_heat233/vtk.h, src/vtk.cpp — ASCII VTK export
- examples/two_material.cpp — two-material plate example
- tests/run_tests.cpp — self-tests

## Build, test, example

    make            # static library libplate_heat233.a + example
    make test       # run self-tests
    make example    # run the two-material example (writes examples/two_material.vtk)

## Mesh format (Gmsh 2.2 ASCII)

- Plane conforming mesh, coordinates in metres, z must be 0.
- Element types: 2-node lines (boundary) and 3-node triangles (domain).
- The first physical tag of a triangle selects the material region; the
  first physical tag of a line selects the boundary condition.
- Node/element ids may be non-contiguous; triangle winding is arbitrary.
- Rejected with a clear error: duplicate ids, references to unknown nodes,
  degenerate triangles, non-manifold edges, boundary lines that are not on
  the outer boundary, non-zero z, unsupported element types.
- Limits: 5000 nodes, 10000 triangles. Overlap-free input is assumed.

## Physics

- Per region tag: Material{k, Q} with k > 0 finite [W/(m K)] and finite
  volumetric source Q [W/m^3]. A tag missing its material fails.
- Per boundary tag, one of:
  - DirichletBC{T}      prescribed temperature [K]
  - FluxBC{q}           outward heat flux [W/m^2], positive out of the domain
  - ConvectionBC{h, Ta} outward flux = h (T - Ta), h > 0 finite
- Boundary tags without a condition are adiabatic.
- Conflicting prescribed temperatures at a shared node fail.
- Every connected component of the mesh must touch at least one Dirichlet
  or convection boundary, otherwise the temperature is not unique and the
  solve fails.

## Method

- P1 linear triangles; the sparse symmetric system is assembled with Eigen.
- Non-zero Dirichlet values are eliminated with right-hand-side adjustment
  (no penalty numbers, no finite-difference shortcuts); the reduced system
  stays symmetric and is solved with SimplicialLDLT.
- The solve reports success only if the linear solve succeeds and the
  maximum absolute residual on free nodes (original system) meets
  Problem::tolerance.
- Output: temperature per original node id, per-triangle heat flux
  -k grad T, max residual, total source power, outward power per boundary
  tag (Dirichlet power recovered from the reactions of the original,
  non-eliminated system), and the overall power imbalance.
- write_vtk exports an ASCII legacy VTK unstructured grid with the
  temperature point field and the heat-flux cell vectors.

## Example

examples/two_material.cpp generates a 0.1 m x 0.1 m plate split into a
copper half (k = 200) and a heated steel half (k = 50, Q = 1e5 W/m^3),
fixes the left edge at 300 K, convects on the right edge
(h = 500 W/(m^2 K), Ta = 290 K), solves and prints the power balance:

    source power [W]: 500
    boundary 10 outward power [W]: 76.9231
    boundary 11 outward power [W]: 423.077
    power imbalance [W]: ~1e-8

