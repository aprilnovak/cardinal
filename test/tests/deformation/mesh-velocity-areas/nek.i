[Mesh]
  type = NekRSMesh
  order = SECOND
  boundary = '2'
  parallel_type = replicated
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh =true
[]

[Problem]
  type = NekRSProblem
  synchronization_interval = parent_app
  has_heat_source = false
  casename = 'pipe'
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Postprocessors]
  [nekbdry_icar]
    type = NekSideIntegral
    field = unity
    boundary = '2'
    execute_on = INITIAL
  []
  [nekbdry_ar]
    type = NekSideIntegral
    field = unity
    boundary = '2'
  []
[]

[Outputs]
  csv = true
  execute_on = 'final'
  show = 'nekbdry_ar nekbdry_icar'
[]
