[Problem]
  type = NekRSStandaloneProblem
  casename = 'pyramid_low'
[]

[Mesh]
  type = NekRSMesh
  volume = true
  exact = true
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [num_elems]
    type = NekMeshInfoPostprocessor
    test_type = num_elems
  []
  [num_nodes]
    type = NekMeshInfoPostprocessor
    test_type = num_nodes
  []
  [num_elems_mirror]
    type = NumElems
  []
  [num_nodes_mirror]
    type = NumNodes
  []
[]
