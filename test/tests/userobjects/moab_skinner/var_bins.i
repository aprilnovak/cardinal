[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../../neutronics/meshes/tet_cube.e
  []
  [id]
    type = SubdomainIDGenerator
    input = sphere
    subdomain_id = 1
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
  [temp_bins]
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
  [temp_bins]
    type = MoabSkinnedBins
    variable = temp_bins
    skinner = moab
    skin_by = variable
  []
[]

[Functions]
  [temp]
    type = ParsedFunction
    value = 400+x*100+100*t
  []
[]

[Materials]
  [mat]
    type = OpenMCDensity
    density = 5.0
    block = '1'
  []
[]

[UserObjects]
  [moab]
    type = MoabUserObject
    bin_varname = "temp"
    material_names = "mat"
    n_bins = 5
    var_min = 400.0
    var_max = 650.0
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
[]

[Outputs]
  exodus = true
  hide = 'temp'
[]
