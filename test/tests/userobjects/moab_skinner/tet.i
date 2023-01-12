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
  type = FEProblem
  solve = false
[]

[UserObjects]
  [moab]
    type = MoabSkinner
    output_full = true
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
[]
