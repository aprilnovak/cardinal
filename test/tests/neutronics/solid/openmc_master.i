[Mesh]
  [sphere]
    # Mesh of a single pebble with outer radius of 1.5 (cm)
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [multiple]
    type = CombinerGenerator
    inputs = sphere
    positions = '0 0 0
                 0 0 4
                 0 0 8'
  []
[]

[Variables]
  [temp]
    initial_condition = 300
  []
[]

[AuxVariables]
  [heat_source]
    family = MONOMIAL
    order = CONSTANT
  []
  [average_temp]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Kernels]
  [hc]
    type = HeatConduction
    variable = temp
  []
  [heat]
    type = CoupledForce
    variable = temp
    v = heat_source
  []
[]

[BCs]
  [outside]
    type = FunctionDirichletBC
    variable = temp
    boundary = '1'
    function = 'axial'
  []
[]

zmin = -1.5
zmax = 9.5

m=${fparse 1000/(zmax-zmin)}
b=${fparse 300-1000/(zmax-zmin)*zmin}

[Functions]
  [axial]
    type = ParsedFunction
    value = '${m}*z+${b}'
  []
[]

[AuxKernels]
  [average_temp]
    variable = average_temp
    type = SpatialUserObjectAux
    user_object = average_temp
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [hc]
    type = GenericConstantMaterial
    prop_values = '0.2' # 20 W/mK -> 0.2 W/cmK
    prop_names = 'thermal_conductivity'
  []
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type -pc_hypre_type'
  num_steps = 2
  petsc_options_value = 'hypre boomeramg'
  dt = 1.0
  nl_abs_tol = 1e-8
[]

[Outputs]
  exodus = true

  [console]
    type = Console
    hide = 'total_volume'
  []
[]

[UserObjects]
  [average_temp]
    type = CardinalNearestPointAverage
    variable = temp
    points_file = pebble_centers_rescaled.txt
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [total_volume]
    type = VolumePostprocessor
    execute_on = 'INITIAL'
  []
  [pebble_mesh_volume]
    type = LinearCombinationPostprocessor
    pp_names = 'total_volume'
    pp_coefs = '${fparse 1.0/3.0}'
    execute_on = 'INITIAL'
  []
  [power_from_openmc]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
  []

  # These are the max/min temperatures in the MOOSE temperature field; because
  # this field is averaged with a nearest point average before being sent to OpenMC,
  # these will be slightly different from the max/min temperatures sent to OpenMC
  [max_T]
    type = ElementExtremeValue
    variable = temp
    execute_on = 'initial timestep_end'
    value_type = max
  []
  [min_T]
    type = ElementExtremeValue
    variable = temp
    execute_on = 'initial timestep_end'
    value_type = min
  []
[]

[MultiApps]
  [openmc]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'openmc.i'
    execute_on = timestep_begin
  []
[]

[Transfers]
  [heat_source_from_openmc]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = openmc
    variable = heat_source
    source_variable = heat_source
  []
  [average_temp_to_openmc]
    type = NearestPointReceiverTransfer
    direction = to_multiapp
    multi_app = openmc
    from_uo = average_temp
    to_uo = average_temp
  []
[]
