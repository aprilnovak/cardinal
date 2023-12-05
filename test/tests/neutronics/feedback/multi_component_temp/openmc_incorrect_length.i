[Mesh]
  type = FileMesh
  file = ../../meshes/pincell.e
[]

[Problem]
  type = OpenMCCellAverageProblem
  power = 500.0
  tally_blocks = '1'
  tally_type = cell
  cell_level = 1

  temperature_variables = 'solid_temp; fluid_temp'
  temperature_blocks = '1'
[]

[Executioner]
  type = Transient
[]
