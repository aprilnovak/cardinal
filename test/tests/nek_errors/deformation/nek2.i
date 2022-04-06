[Mesh]
  type = NekRSMesh
  order = SECOND
  boundary = '2'
  parallel_type = replicated
  displacements = 'disp_x disp_y disp_z'
  moving_mesh = true
[]

[Problem]
  type = NekRSProblem
  casename = 'err2'
  minimize_transfers_in = true
[]

[AuxVariables]
  [dummy]
  []
[]

# This AuxVariable and AuxKernel is only here to get the postprocessors
# to evaluate correctly. This can be deleted after MOOSE issue #17534 is fixed.
[AuxKernels]
  [dummy]
    type = ConstantAux
    variable = dummy
    value = 0.0
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = NekTimeStepper
  []
  [Quadrature]
    type = GAUSS_LOBATTO
    order = SECOND
  []
[]

[Postprocessors]
  [nekbdry_icar]
    type = NekSideIntegral
    field = unity
    boundary = '2'
    use_displaced_mesh = true
    execute_on = INITIAL
  []
  [nekbdry_ar]
    type = NekSideIntegral
    field = unity
    boundary = '2'
    use_displaced_mesh = true
  []
[]

[Outputs]
  exodus = false
  execute_on = 'final'
  hide = 'flux_integral'
[]

