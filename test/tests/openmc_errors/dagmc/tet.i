[Mesh]
  [cube]
    type = FileMeshGenerator
    file = ../../neutronics/meshes/tet_cube.e
  []
  [add_block]
    type = SubdomainIDGenerator
    input = cube
    subdomain_id = 1
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  power = 100.0
  solid_blocks = '1'
  tally_blocks = '1'
  solid_cell_level = 0
  tally_type = cell
  tally_name = heat_source
  initial_properties = xml

  skinning_user_object = moab
[]

[Materials]
  [mat]
    type = OpenMCDensity
    density = 100.0
  []
[]

[UserObjects]
  [moab]
    type = MoabUserObject
    material_names = 'mat'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
