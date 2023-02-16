[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../meshes/tet_cube.e
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  tally_type = cell
  solid_cell_level = 0
  tally_blocks = '1'
  solid_blocks = '1'
  power = 1000.0
  skinner = moab
[]

[UserObjects]
  [moab]
    type = NearestNodeNumberUO
    point = '0.0 0.0 0.0'
  []
[]

[Executioner]
  type = Steady
[]
