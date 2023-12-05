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

  initial_properties = xml
  temperature_blocks = '1; 3'
  temperature_materials = 'mat'
[]

[Executioner]
  type = Transient
[]
