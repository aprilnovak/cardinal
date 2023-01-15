[Mesh]
  [cube]
    type = FileMeshGenerator
    file = ../../neutronics/meshes/tet_cube.e
  []
  [id1]
    type = ParsedSubdomainMeshGenerator
    input = cube
    combinatorial_geometry = 'y < 0.0'
    block_id = 1
  []
  [id2]
    type = ParsedSubdomainMeshGenerator
    input = id1
    combinatorial_geometry = 'y >= 0.0'
    block_id = 2
  []
[]

[Problem]
  type = FEProblem
  solve = false
[]

[AuxVariables]
  [temp]
    family = MONOMIAL
    order = CONSTANT
  []
  [mat_bins]
    family = MONOMIAL
    order = CONSTANT
  []
  [all_bins]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [temp]
    type = FunctionAux
    variable = temp
    function = temp
    execute_on = timestep_begin
  []
  [mat_bins]
    type = MoabSkinnedBins
    variable = mat_bins
    skinner = moab
    skin_by = material
  []
  [all_bins]
    type = MoabSkinnedBins
    variable = all_bins
    skinner = moab
  []
[]

[Functions]
  [temp]
    type = ParsedFunction
    value = 400+x*100+100*t
  []
[]

[Materials]
  [mat1]
    type = OpenMCDensity
    density = 5.0
    block = '1'
  []
  [mat2]
    type = OpenMCDensity
    density = 5.0
    block = '2'
  []
[]

[UserObjects]
  [moab]
    type = MoabUserObject
    verbose = true
    material_names = "mat2 mat1"

    temperature = "temp"
    n_temperature_bins = 5
    temperature_min = 400.0
    temperature_max = 650.0

    output_skins = true
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
  execute_on = final
  hide = 'temp'
[]
