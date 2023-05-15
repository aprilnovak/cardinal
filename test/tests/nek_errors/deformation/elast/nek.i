[Mesh]
  type = NekRSMesh
  volume = true
  parallel_type = replicated

  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  type = NekRSProblem
  casename = 'elast'
  has_heat_source = false
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = NekTimeStepper
  []
[]
