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
  [rho]
    family = MONOMIAL
    order = CONSTANT
  []
  [temp]
    family = MONOMIAL
    order = CONSTANT
  []
  [temp_bins]
    family = MONOMIAL
    order = CONSTANT
  []
  [den_bins]
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
  [rho]
    type = FunctionAux
    variable = rho
    function = rho
    execute_on = timestep_begin
  []
  [temp_bins]
    type = MoabSkinnedBins
    variable = temp_bins
    skinner = moab
    skin_by = variable
  []
  [den_bins]
    type = MoabSkinnedBins
    variable = den_bins
    skinner = moab
    skin_by = density
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
    value = 400+x*100+t*100
  []
  [rho]
    type = ParsedFunction
    value = 400+z*100+t*100
  []
[]

[Materials]
  [mat1]
    type = OpenMCDensity
    density = 550.0
    block = '1'
  []
  [mat2]
    type = OpenMCDensity
    density = 550.0
    block = '2'
  []
[]

[UserObjects]
  [moab]
    type = MoabUserObject
    verbose = true
    material_names = "mat2 mat1"

    density = "rho"
    rel_den_min = -0.5
    rel_den_max =  0.5
    n_density_bins = 5

    temperature = "temp"
    temperature_min = 400.0
    temperature_max = 650.0
    n_temperature_bins = 5

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
