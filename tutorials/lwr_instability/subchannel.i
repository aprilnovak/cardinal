!include common.i

# water specific heat, for estimating dT
Cp = 4108

flow_area = ${fparse 2 * pin_pitch * 2 * pin_pitch - 4 * pi * pin_diameter^2 / 4}
mass_flux_in = ${fparse mdot / flow_area}

[GlobalParams]

  # number of channels in the x-direction (number of pins in x plus 1)
  nx = 3

  # number of channels in the y-direction (number of pins in y plus 1)
  ny = 3

  # number of cells in the vertical direction
  n_cells = ${n_layers}

  n_blocks = 1
  pitch = ${pin_pitch}
  pin_diameter = ${pin_diameter}
  side_gap = 0.0
  heated_length = ${heated_length}

  # do not model any grid spacers
  spacer_z = '0.0'
  spacer_k = '0.0'
[]

[QuadSubChannelMesh]
  [assembly]
    type = SCMQuadAssemblyMeshGenerator
  []
[]

[FluidProperties]
  [water]
    type = Water97FluidProperties
  []
[]

[SubChannel]
  type = QuadSubChannel1PhaseProblem
  fp = water
  P_tol = 1e-6
  T_tol = 1e-5
  compute_density = true
  compute_viscosity = true
  compute_power = true
  P_out = ${P_outlet}
  implicit = true
  segregated = false
  friction_closure = 'MATRA'
  pin_HTC_closure = 'Dittus-Boelter'
  full_output = true
  mixing_closure ='constant_beta'
[]

[SCMClosures]
  [MATRA]
    type = SCMFrictionMATRA
  []
  [Dittus-Boelter]
    type = SCMHTCDittusBoelter
  []
  [constant_beta]
    type = SCMMixingConstantBeta
    beta = 0.006
    CT = 2.6
  []
[]

[ICs]
  [T_ic]
    type = ConstantIC
    variable = T
    value = ${T_inlet}
  []
  [Dpin_IC]
    type = ConstantIC
    variable = Dpin
    value = ${fparse pin_diameter}
  []
  [q_prime_IC]
    type = ConstantIC
    variable = q_prime
    value = ${fparse power/4/heated_length}
  []
  [Viscosity_ic]
    type = ViscosityIC
    variable = mu
    p = ${P_outlet}
    T = T
    fp = water
  []
  [rho_ic]
    type = RhoFromPressureTemperatureIC
    variable = rho
    p = ${P_outlet}
    T = T
    fp = water
  []
  [h_ic]
    type = SpecificEnthalpyFromPressureTemperatureIC
    variable = h
    p = ${P_outlet}
    T = T
    fp = water
  []
[]

[AuxKernels]
  [T_in_bc]
    type = ConstantAux
    variable = T
    boundary = inlet
    value = ${T_inlet}
    execute_on = 'timestep_begin'
  []
  [mdot_in_bc]
    type = SCMMassFlowRateAux
    variable = mdot
    boundary = inlet
    area = S
    mass_flux = ${mass_flux_in}
    execute_on = 'timestep_begin'
  []
[]

[Postprocessors]
  [power]
    type = SCMPinPowerPostprocessor
    execute_on = 'transfer'
  []
  [inlet_temp]
    type = SCMPlanarMean
    variable = T
    height = 0.0
  []
  [outlet_temp]
    type = SCMPlanarMean
    variable = T
    height = ${fparse heated_length}
  []
  [expected_dT]
    type = ParsedPostprocessor
    expression = 'power/${mdot}/${Cp}'
    pp_names = 'power'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

[Executioner]
  type = Steady
[]
