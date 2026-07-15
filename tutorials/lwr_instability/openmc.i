!include common.i

[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  power = ${power}
  lowest_cell_level = 1
  scaling = 100
  temperature_blocks = 'fuel clad water'
  density_blocks = 'water'
  initial_properties = 'xml'

  max_batches = 1000
  batch_interval = 20

  [Tallies]
    [power]
      type = CellTally
      score = 'kappa_fission'
      trigger = rel_err
      trigger_threshold = 2e-2
      block = 'fuel'
    []
  []
[]

[AuxVariables]
  [q_bottom_half]
    family = MONOMIAL
    order = CONSTANT
    block = 'fuel'
  []
[]

[AuxKernels]
  [q_bottom_half]
    type = ParsedAux
    variable = q_bottom_half
    expression = 'if (z < ${fparse heated_length/2}, kappa_fission, 0.0)'
    coupled_variables = 'kappa_fission'
    use_xyzt = true
  []
[]

[Executioner]
  type = Transient
  num_steps = 10
[]

[Outputs]
  exodus = true
  csv = true
  hide = 'q_bottom_half total_power_in_bottom_half'
[]

[MultiApps]
  [conduction]
    type = TransientMultiApp
    input_files = 'conduction.i'
    sub_cycling = true
    execute_on = timestep_end
  []
  [subchannel]
    type = FullSolveMultiApp
    input_files = 'subchannel.i'
    max_procs_per_app = 1
    execute_on = timestep_begin
  []
[]

[Transfers]
  [power_to_conduction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = conduction
    source_variable = kappa_fission
    variable = heat_source
    from_postprocessors_to_be_preserved = openmc_power_integral
    to_postprocessors_to_be_preserved = conduction_power_integral
  []
  [solid_temperature_from_conduction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = conduction
    source_variable = T
    variable = temp
    to_blocks = 'fuel clad'
  []

  [linear_heat_rate_to_subchannel]
    type = MultiAppGeneralFieldNearestLocationTransfer
    variable = q_prime
    source_variable = q_prime
    from_multi_app = conduction
    to_multi_app = subchannel
    greedy_search = true
    use_bounding_boxes = false
    to_blocks = 'fuel_pins'
    #from_postprocessors_to_be_preserved = conduction_power_integral
    #to_postprocessors_to_be_preserved = power
  []
  [fluid_temperature_from_subchannel]
    type = MultiAppGeneralFieldNearestLocationTransfer
    source_variable = T
    variable = temp
    from_multi_app = subchannel
    greedy_search = true
    use_bounding_boxes = false
    to_blocks = 'water'
    from_blocks = 'subchannel'
  []
  [fluid_density_from_subchannel]
    type = MultiAppGeneralFieldNearestLocationTransfer
    source_variable = rho
    variable = density
    from_multi_app = subchannel
    greedy_search = true
    use_bounding_boxes = false
    to_blocks = 'water'
    from_blocks = 'subchannel'
  []
  [clad_surface_temperature_to_conduction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = subchannel
    to_multi_app = conduction
    source_variable = Tpin
    variable = T_wall
  []
[]

[Postprocessors]
  [openmc_power_integral]
    type = ElementIntegralVariablePostprocessor
    variable = kappa_fission
    execute_on = 'transfer timestep_end'
    block = 'fuel'
  []
  [total_power_in_bottom_half]
    type = ElementIntegralVariablePostprocessor
    variable = q_bottom_half
    block = 'fuel'
  []
  [fraction_power_in_bottom]
    type = ParsedPostprocessor
    expression = 'total_power_in_bottom_half / openmc_power_integral'
    pp_names = 'total_power_in_bottom_half openmc_power_integral'
  []
  [min_fuel_T]
    type = ElementExtremeValue
    variable = temp
    block = 'fuel'
    value_type = min
  []
  [max_fuel_T]
    type = ElementExtremeValue
    variable = temp
    block = 'fuel'
  []
  [average_rel_err]
    type = TallyRelativeError
    value_type = average
  []
  [max_rel_err]
    type = TallyRelativeError
    value_type = max
  []
[]

[UserObjects]
  [layered_average]
    type = NearestPointLayeredAverage
    variable = kappa_fission
    block = 'fuel'
    direction = z
    points = '0 0 0'
    num_layers = ${n_layers}
  []
[]

[VectorPostprocessors]
  [axial_power]
    type = SpatialUserObjectVectorPostprocessor
    userobject = layered_average
  []
[]
