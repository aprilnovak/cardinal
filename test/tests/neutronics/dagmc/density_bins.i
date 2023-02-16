[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../meshes/tet_cube.e
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  tally_type = cell
  fluid_cell_level = 0
  fluid_blocks = '1'
  power = 1000.0
  skinner = moab
[]

[UserObjects]
  [moab]
    type = MoabSkinner
    temperature = temp
    temperature_min = 0.0
    temperature_max = 900.0
    n_temperature_bins = 1
    build_graveyard = true

    density = density
    density_min = 0.0
    density_max = 1000.0
    n_density_bins = 4
  []
[]

[Executioner]
  type = Steady
[]
