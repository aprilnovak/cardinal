!include common.i

[Mesh]
  [pin]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    ring_radii = '${fparse pellet_diameter/2} ${fparse pin_diameter/2}'
    ring_intervals = '7 3'
    num_sectors_per_side = '6 6 6 6'
    polygon_size = ${fparse 0.5 * pin_pitch}

    ring_block_ids = '1 2'
    ring_block_names = 'fuel clad'
    background_block_ids = '3'
    background_block_names = 'water'
    quad_center_elements = true
    flat_side_up = true
  []
  [assembly]
    type = PatternedCartesianMeshGenerator
    inputs = 'pin'
    pattern = '0 0;
               0 0'
    square_size = ${fparse 2*pin_pitch}
    background_block_id = '3'
    background_block_name = 'water'
  []
  [extrude]
    type = AdvancedExtruderGenerator
    input = assembly
    direction = '0 0 1'
    num_layers = '${n_layers}'
    heights = '${heated_length}'
  []
[]
