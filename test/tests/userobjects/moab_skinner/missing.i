[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../../neutronics/meshes/tet_cube.e
  []
[]

[Problem]
  type = FEProblem
  solve = false
[]

[Variables]
  [temp]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[UserObjects]
  [moab]
    type = MoabUserObject
    verbose = true
    material_names = "mat"

    temperature = "temp"
    n_temperature_bins = 5
    temperature_min = 400.0
    temperature_max = 650.0
  []
[]

[Executioner]
  type = Transient
[]

[Outputs]
  exodus = true
[]
