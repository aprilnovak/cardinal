!include common.i

[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
  [delete_coolant]
    type = BlockDeletionGenerator
    input = file
    block = 'water'
  []
[]

[Variables]
  [T]
    initial_condition = ${T_inlet}
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [heat_source_fuel]
    type = CoupledForce
    variable = T
    v = heat_source
    block = 'fuel'
  []
[]

[Materials]
  [k_fuel]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = '(100/(7.5408+17.629*t/1000+3.6142*(t/1000)^2) + 6400/((t/1000)^(5/2))*exp(-16.35/(t/1000)))*1/(1-(2.6-0.5*t/1000)*0.05)'
    temperature = T
    block = 'fuel'
  []
  [dk_fuel_dT]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity_dT'
    prop_values = '0.0'
    block = 'clad'
  []
  [k_clad]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '20.0'
    block = 'clad'
  []
[]

[AuxVariables]
  [T_wall]
    initial_condition = ${T_inlet}
  []
  [heat_source]
    # OpenMC computes the variable as constant monomial, so we can receive the
    # field into exactly the same type here (not required, but this will be clearer to
    # visualize in paraview)
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1e8
    block = 'fuel'
  []
  [q]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0.0
    block = 'fuel'
  []
  [q_prime]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0.0
    block = 'fuel'
  []
[]

[BCs]
  [cladding_outer_bc]
    type = MatchedValueBC
    variable = T
    v = T_wall
    boundary = '3'
  []
[]

[Executioner]
  type = Transient
[]

[AuxKernels]
  [q]
    type = SpatialUserObjectAux
    variable = q
    user_object = q_prime_uo
    block = 'fuel'
  []
  [q_prime] # divide by height of each averaging layer to get W/m from W
    type = ParsedAux
    variable = q_prime
    coupled_variables = 'q'
    expression = 'q / ${fparse heated_length/n_layers}'
    block = 'fuel'
  []
[]

[UserObjects]
  [q_prime_uo]
    type = NearestPointLayeredIntegral
    variable = heat_source
    block = 'fuel'
    direction = z
    points = '${fparse -pin_pitch/2} ${fparse -pin_pitch/2} 0.0
              ${fparse  pin_pitch/2} ${fparse -pin_pitch/2} 0.0
              ${fparse  pin_pitch/2} ${fparse  pin_pitch/2} 0.0
              ${fparse -pin_pitch/2} ${fparse  pin_pitch/2} 0.0'
    num_layers = ${n_layers}
    execute_on = 'initial timestep_begin'
  []
[]

[Postprocessors]
  [conduction_power_integral]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    block = 'fuel'
    execute_on = 'transfer'
  []
[]

[Outputs]
  exodus = true
[]
