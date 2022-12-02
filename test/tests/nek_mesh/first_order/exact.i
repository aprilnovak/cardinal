[Problem]
  type = NekRSStandaloneProblem
  casename = 'pyramid'
[]

[Mesh]
  type = NekRSMesh
  exact = true
  boundary = '1 3'
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

# The points provided to these postprocessors are the centroids of the elements that
# we wish to print the node coordinates for.
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
