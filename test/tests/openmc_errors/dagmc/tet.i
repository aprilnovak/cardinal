[Mesh]
  [cube]
    type = FileMeshGenerator
    file = ../../neutronics/meshes/tet_cube.e
  []
  [add_block]
    type = SubdomainIDGenerator
    input = cube
    subdomain_id = 0
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  power = 100.0
  solid_blocks = '0'
  tally_blocks = '0'
  solid_cell_level = 0
  tally_type = cell
  tally_name = heat_source
  initial_properties = xml

  skinning_user_object = moab
[]

[UserObjects]
  [moab]
    type = MoabUserObject
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
