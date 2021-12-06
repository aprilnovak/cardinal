[Mesh]
  [circle]
    type = AnnularMeshGenerator
    nr = 10
    nt = 25
    rmin = 0
    rmax = 0.4e-2
    growth_r = -1.3
  []
  [cylinder]
    type = FancyExtruderGenerator
    input = circle
    heights = '0.05 0.7 0.05'
    num_layers = '5 10 5'
    direction = '0 0 1'
  []
  [transform]
    type = TransformGenerator
    input = cylinder
    transform = translate
    vector_value = '0 0 -0.4'
  []
[]

[Variables]
  [temperature]
    initial_condition = 500.0
  []
[]

[Kernels]
  [diffusion]
    type = HeatConduction
    variable = temperature
    diffusion_coefficient = thermal_conductivity
  []
  [source]
    type = CoupledForce
    variable = temperature
    v = source
  []
[]

[BCs]
  [interface]
    type = MatchedValueBC
    variable = temperature
    v = nek_temp
    boundary = '1'
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 0.01
  nl_abs_tol = 1e-8
[]

[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'nek.i'
    sub_cycling = true
    execute_on = timestep_end
  []
[]

[Transfers]
  [temperature]
    type = MultiAppNearestNodeTransfer
    source_variable = temp
    direction = from_multiapp
    multi_app = nek
    variable = nek_temp
  []
  [flux]
    type = MultiAppNearestNodeTransfer
    source_variable = flux
    direction = to_multiapp
    multi_app = nek
    variable = avg_flux
    source_boundary = '1'
  []
  [flux_integral]
    type = MultiAppPostprocessorTransfer
    to_postprocessor = flux_integral
    direction = to_multiapp
    from_postprocessor = flux_integral
    multi_app = nek
  []
[]

[AuxVariables]
  [flux]
    family = MONOMIAL
    order = CONSTANT
  []
  [nek_temp]
    initial_condition = 628.15
  []
  [source]
    initial_condition = 2e6
  []
[]

[AuxKernels]
  [flux]
    type = DiffusionFluxAux
    variable = flux
    diffusion_variable = temperature
    component = normal
    diffusivity = thermal_conductivity
    boundary = '1'
  []
[]

[Materials]
  [k]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '1.5'
  []
[]

[Postprocessors]
  [flux_integral]
    type = SideIntegralVariablePostprocessor
    variable = flux
    boundary = '1'
  []
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
  hide = 'flux_integral'

  [pg]
    type = PerfGraphOutput
    execute_on = 'timestep_end'
  []
[]

# Various postprocessors for determining runtime of different parts of the solve
[Postprocessors]
  [runtime] # total runtime of the simulation (both applications)
    type = PerfGraphData
    data_type = TOTAL
    section_name = 'Root'
    execute_on = 'timestep_end'
  []
  [moose_time]
    type = PerfGraphData
    data_type = TOTAL
    section_name = 'Transient::PicardSolve'
    execute_on = 'final timestep_end'
  []
  [nek_time]
    type = PerfGraphData
    data_type = TOTAL
    section_name = 'FEProblem::execMultiApps'
    must_exist = false
    execute_on = 'timestep_end'
 []
  [transfer_time]
    type = PerfGraphData
    data_type = TOTAL
    section_name = 'FEProblem::execMultiAppTransfers'
    must_exist = false
    execute_on = 'final timestep_end'
  []
[]
