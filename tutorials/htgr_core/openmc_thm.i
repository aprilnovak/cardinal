# This input file runs coupled OpenMC Monte Carlo transport, MOOSE heat
# conduction, and THM fluid flow and heat transfer.
# This input should be run with:
#
# cardinal-opt -i common_input.i openmc_thm.i

N = 2

num_layers_for_THM = 50
num_layers = 50

fluid_blocks = '102'
solid_blocks = '0 1 2 4 5'
tally_blocks = 'compacts'

[Mesh]
  [coolant_face]
    type = AnnularMeshGenerator
    nr = 1
    nt = 4
    rmin = 0.0
    rmax = ${fparse channel_diameter / 2.0}
    quad_subdomain_id = 101
    tri_subdomain_id = 102
  []
  [extrude]
    type = AdvancedExtruderGenerator
    input = coolant_face
    num_layers = ${num_layers_for_THM}
    direction = '0 0 1'
    heights = '${height}'
    top_boundary = '300' # inlet
    bottom_boundary = '400' # outlet
  []
  [repeat_fluid]
    type = CombinerGenerator
    inputs = extrude
    positions_file = coolant_channel_positions.txt
  []
  [solid]
    type = FileMeshGenerator
    file = solid_in.e
  []
  [add]
    type = CombinerGenerator
    inputs = 'solid repeat_fluid'
  []
[]

[AuxVariables]
  [cell_temperature]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [cell_temperature]
    type = CellTemperatureAux
    variable = cell_temperature
  []
  [density]
    type = FluidDensityAux
    variable = density
    p = ${outlet_P}
    T = temp
    fp = helium
    execute_on = 'timestep_begin linear'
  []
[]

[ICs]
  [temp]
    type = FunctionIC
    variable = thm_temp
    function = temp_ic
  []
  [solid_temp]
    type = FunctionIC
    variable = solid_temp
    function = temp_ic
  []
  [heat_source]
    type = ConstantIC
    variable = heat_source
    block = ${tally_blocks}
    value = ${fparse power / (pi * compact_diameter * compact_diameter / 4.0 * height * n_bundles * n_fuel_compacts_per_block)}
  []
[]

[FluidProperties]
  [helium]
    type = IdealGasFluidProperties
    molar_mass = 4e-3
    gamma = 1.668282 # should correspond to  Cp = 5189 J/kg/K
    k = 0.2556
    mu = 3.22639e-5
  []
[]

[Functions]
  [temp_ic]
    type = ParsedFunction
    value = '${inlet_T} + z / ${height} * ${power} / ${mdot} / ${fluid_Cp}'
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  output = 'unrelaxed_tally_std_dev'

  # optimizations to increase tracking rate
  check_tally_sum = false
  normalize_by_global_tally = false
  assume_separate_tallies = true

  power = ${power}
  scaling = 100.0
  solid_blocks = ${solid_blocks}
  fluid_blocks = ${fluid_blocks}
  tally_blocks = ${tally_blocks}
  tally_type = cell
  solid_cell_level = 2
  fluid_cell_level = 2

  relaxation = dufek_gudowski
  first_iteration_particles = 50000

  symmetry_plane_normal = '-1.0 0.0 0.0'
  symmetry_axis = '0.0 0.0 1.0'
  symmetry_angle = 30.0

  temperature_variables = 'solid_temp solid_temp solid_temp solid_temp solid_temp thm_temp'
  temperature_blocks = '${solid_blocks} ${fluid_blocks}'

  inactive_batches = 150
  batches = 500
[]

[Executioner]
  type = Transient
  dt = ${N}
  num_steps = 5

  steady_state_detection = true
  check_aux = true
  steady_state_tolerance = 1e-3
[]

[Postprocessors]
  [heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    execute_on = 'transfer initial timestep_end'
  []
  [max_tally_err]
    type = FissionTallyRelativeError
    value_type = max
  []
  [k]
    type = KEigenvalue
  []
  [k_std_dev]
    type = KStandardDeviation
  []
  [min_power]
    type = ElementExtremeValue
    variable = heat_source
    value_type = min
    block = ${tally_blocks}
  []
  [max_power]
    type = ElementExtremeValue
    variable = heat_source
    value_type = max
    block = ${tally_blocks}
  []
[]

[MultiApps]
  [bison]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'solid_thm.i'
    execute_on = timestep_begin
    sub_cycling = true
  []
[]

[Transfers]
  [solid_temp_to_openmc]
    type = MultiAppCopyTransfer
    source_variable = T
    variable = solid_temp
    from_direction = from_multiapp
  []
  [source_to_bison]
    type = MultiAppCopyTransfer
    source_variable = heat_source
    variable = power
    to_multi_app = bison
    from_postprocessors_to_be_preserved = heat_source
    to_postprocessors_to_be_preserved = power
  []
  [temp_from_thm]
    type = MultiAppCopyTransfer
    source_variable = thm_bulk_temp
    from_multi_app = bison
    variable = thm_temp
  []
[]

[UserObjects]
  [average_power_axial]
    type = LayeredAverage
    variable = heat_source
    direction = z
    num_layers = ${num_layers}
    block = ${tally_blocks}
  []
[]

[VectorPostprocessors]
  [power_avg]
    type = SpatialUserObjectVectorPostprocessor
    userobject = average_power_axial
  []
[]

[Outputs]
  [out]
    type = Exodus
    hide = 'density temp solid_temp thm_temp'
  []

  [csv]
    type = CSV
    file_base = 'csv_thm/openmc_thm'
  []

  [screen]
    type = Console
    hide = 'max_power min_power k_std_dev'
  []
[]
