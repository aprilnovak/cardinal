# Quantities copied from common_input.i, since MOOSE does not yet support
# passing multiple input files from the main app level to sub-apps

compact_diameter = 0.0127                # diameter of fuel compacts (m)
channel_diameter = 0.016                 # diameter of the coolant channels (m)
kernel_radius = 214.85e-6                # fissile kernel outer radius (m)
buffer_radius = 314.85e-6                # buffer outer radius (m)
iPyC_radius = 354.85e-6                  # inner PyC outer radius (m)
SiC_radius = 389.85e-6                   # SiC outer radius (m)
oPyC_radius = 429.85e-6                  # outer PyC outer radius (m)
buffer_k = 0.5                           # buffer thermal conductivity (W/m/K)
PyC_k = 4.0                              # PyC thermal conductivity (W/m/K)
SiC_k = 13.9                             # SiC thermal conductivity (W/m/K)
kernel_k = 3.5                           # fissil kernel thermal conductivity (W/m/K)
matrix_k = 15.0                          # graphite matrix thermal conductivity (W/m/K)
inlet_T = 598.0                          # inlet fluid temperature (K)
power = 200e6                            # full core power (W)
mdot = 117.3                             # fluid mass flowrate (kg/s)
height = 6.343                           # height of the full core (m)
fluid_Cp = 5189.0                        # fluid isobaric specific heat (J/kg/K)
triso_pf = 0.15                          # TRISO packing fraction (%)
n_bundles = 12                           # number of bundles in the full core
n_fuel_compacts_per_block = 210          # number of fuel compacts per assembly

# compute the volume fraction of each TRISO layer in a TRISO particle
# for use in computing average thermophysical properties
kernel_fraction = ${fparse kernel_radius^3 / oPyC_radius^3}
buffer_fraction = ${fparse (buffer_radius^3 - kernel_radius^3) / oPyC_radius^3}
ipyc_fraction = ${fparse (iPyC_radius^3 - buffer_radius^3) / oPyC_radius^3}
sic_fraction = ${fparse (SiC_radius^3 - iPyC_radius^3) / oPyC_radius^3}
opyc_fraction = ${fparse (oPyC_radius^3 - SiC_radius^3) / oPyC_radius^3}

num_layers_for_THM = 50                 # number of layers in THM mesh
num_layers_for_plots = 50               # number of layers to average fields over for plotting

fluid_blocks = '102'
solid_blocks = '0 1 2 4 5'

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

[Problem]
  type = FEProblem
  material_coverage_check = false
[]

[Variables]
  [T]
    initial_condition = ${inlet_T}
  []
[]

[Kernels]
  [diffusion]
    type = HeatConduction
    variable = T
    block = ${solid_blocks}
  []
  [source]
    type = CoupledForce
    variable = T
    v = power
    block = 'compacts'
  []
  [null]
    type = NullKernel
    variable = T
    block = ${fluid_blocks}
  []
[]

[BCs]
  [pin_outer]
    type = MatchedValueBC
    variable = T
    v = thm_temp
    boundary = 'fluid_solid_interface'
  []
  [vessel_outer]
    type = ConvectiveFluxFunction
    variable = T
    T_infinity = ${fparse 30 + 273}
    coefficient = 15.0
    boundary = 'vessel_outer'
  []
[]

[Functions]
  [k_graphite]
    type = ParsedFunction
    value = '${matrix_k}'
  []
  [k_TRISO]
    type = ParsedFunction
    value = '${kernel_fraction} * ${kernel_k} + ${buffer_fraction} * ${buffer_k} + ${fparse ipyc_fraction + opyc_fraction} * ${PyC_k} + ${sic_fraction} * ${SiC_k}'
  []
  [k_compacts]
    type = ParsedFunction
    value = '${triso_pf} * k_TRISO + ${fparse 1.0 - triso_pf} * k_graphite'
    vars = 'k_TRISO k_graphite'
    vals = 'k_TRISO k_graphite'
  []
  [k_b4c]
    type = ParsedFunction
    value = 13.62
  []
  [axial_fluid_temp]
    type = ParsedFunction
    value = '${inlet_T} + z / ${height} * ${power} / ${mdot} / ${fluid_Cp}'
  []
[]

[Materials]
  [graphite]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_graphite
    temp = T
    block = 'graphite reflector graphite_centers'
  []
  [compacts]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_TRISO
    temp = T
    block = 'compacts'
  []
  [poison]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_b4c
    temp = T
    block = 'poison'
  []
[]

[Postprocessors]
  [flux_integral] # evaluate the total heat flux for normalization
    type = SideDiffusiveFluxIntegral
    diffusivity = thermal_conductivity
    variable = T
    boundary = 'fluid_solid_interface'
  []
  [max_fuel_T]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'compacts'
  []
  [max_block_T]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'graphite'
  []
  [max_nonfuel_block_T]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'graphite_centers'
  []
  [max_reflector_T]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'reflector'
  []
  [power]
    type = ElementIntegralVariablePostprocessor
    variable = power
    block = 'compacts'
    execute_on = 'transfer initial'
  []
[]

[AuxVariables]
  [thm_temp]
  []
  [thm_bulk_temp]
    block = ${fluid_blocks}
  []
  [flux]
    family = MONOMIAL
    order = CONSTANT
  []
  [power]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${fparse power / (pi * compact_diameter * compact_diameter / 4.0 * height * n_bundles * n_fuel_compacts_per_block)}
  []
[]

[AuxKernels]
  [flux]
    type = DiffusionFluxAux
    diffusion_variable = T
    component = normal
    diffusivity = thermal_conductivity
    variable = flux
    boundary = 'fluid_solid_interface'
  []
[]

[ICs]
  [thm_temp]
    type = FunctionIC
    variable = thm_temp
    function = axial_fluid_temp
  []
[]

[MultiApps]
  [thm]
    type = FullSolveMultiApp
    app_type = CardinalApp
    input_files = 'thm.i'
    execute_on = timestep_end
    max_procs_per_app = 1
    bounding_box_padding = '${fparse 1.1 * channel_diameter / 2.0} ${fparse 1.1 * channel_diameter / 2.0} 0.0'
    positions_file = coolant_channel_positions.txt
    output_in_position = true
  []
[]

[Transfers]
  [q_wall_to_thm]
    type = MultiAppUserObjectTransfer
    variable = q_wall
    to_multi_app = thm
    user_object = q_wall_avg
  []
  [T_wall_from_thm]
    type = MultiAppUserObjectTransfer
    variable = thm_temp
    from_multi_app = thm
    user_object = T_wall
    boundary = 'fluid_solid_interface'
  []
  [T_bulk_from_thm]
    type = MultiAppUserObjectTransfer
    variable = thm_bulk_temp
    from_multi_app = thm
    user_object = T_bulk
    block = ${fluid_blocks}
  []
[]

[AuxVariables]
  [q_wall]
  []
[]

[AuxKernels]
  [q_wall]
    type = SpatialUserObjectAux
    variable = q_wall
    user_object = q_wall_avg
  []
[]

[UserObjects]
  [q_wall_avg]
    type = NearestPointLayeredSideAverage
    boundary = fluid_solid_interface
    variable = flux

    # Note: make this to match the num_elems in the channel
    direction = z
    num_layers = ${num_layers_for_THM}
    points_file = coolant_channel_positions.txt

    direction_min = 0.0
    direction_max = ${height}
  []
  [avg_fuel]
    type = LayeredAverage
    variable = T
    direction = z
    num_layers = ${num_layers_for_plots}
    block = 'compacts'
  []
  [avg_block]
    type = LayeredAverage
    variable = T
    direction = z
    num_layers = ${num_layers_for_plots}
    block = 'graphite graphite_centers'
  []
  [avg_reflector]
    type = LayeredAverage
    variable = T
    direction = z
    num_layers = ${num_layers_for_plots}
    block = 'reflector'
  []
  [avg_fluid_bulk]
    type = LayeredAverage
    variable = thm_bulk_temp
    direction = z
    num_layers = ${num_layers_for_plots}
    block = ${fluid_blocks}
  []
  [avg_fluid_wall]
    type = LayeredSideAverage
    variable = thm_temp
    direction = z
    num_layers = ${num_layers_for_plots}
    boundary = 'fluid_solid_interface'
  []
[]

[Executioner]
  type = Transient
  dt = 1.0
  num_steps = 5

  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-16
  petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -pc_hypre_type'
[]

[VectorPostprocessors]
  [fuel]
    type = SpatialUserObjectVectorPostprocessor
    userobject = avg_fuel
  []
  [block]
    type = SpatialUserObjectVectorPostprocessor
    userobject = avg_block
  []
  [reflector]
    type = SpatialUserObjectVectorPostprocessor
    userobject = avg_reflector
  []
  [fluid_bulk]
    type = SpatialUserObjectVectorPostprocessor
    userobject = avg_fluid_bulk
  []
  [fluid_wall]
    type = SpatialUserObjectVectorPostprocessor
    userobject = avg_fluid_wall
  []
[]

[Outputs]
  print_linear_residuals = false
  exodus = true

  [csv]
    type = CSV
    file_base = 'csv_thm/solid_out'
  []
[]
