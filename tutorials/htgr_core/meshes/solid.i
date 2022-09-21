# Mesh for the solid regions, run with:
#   cardinal-opt -i solid.i common_input.i --mesh-only

num_layers = 50 # number of axial layers in the mesh

[Mesh]
  [fuel_pin]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    polygon_size = ${fparse fuel_to_coolant_distance / 2.0}
    ring_radii =  '0.005499 ${fparse compact_diameter / 2.0}'
    ring_intervals = '1 1'
    num_sectors_per_side = '4 4 4 4 4 4'
    ring_block_ids = '2 2'
    ring_block_names = 'compacts compacts'
    background_block_ids = '1'
    background_block_names = 'graphite'
    background_intervals = 2
    quad_center_elements = true
  []
  [coolant_pin]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    polygon_size = ${fparse fuel_to_coolant_distance / 2.0}
    ring_radii = '${fparse channel_diameter / 2.0}'
    ring_intervals = '2'
    num_sectors_per_side = '4 4 4 4 4 4'
    ring_block_ids = '101 101'
    ring_block_names = 'coolant coolant'
    background_block_ids = '1'
    background_block_names = 'graphite'
    interface_boundary_id_shift = 100
    background_intervals = 2
    quad_center_elements = true
  []
  [poison_pin]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    polygon_size = ${fparse fuel_to_coolant_distance / 2.0}
    ring_radii = '${fparse compact_diameter / 2.0}'
    ring_intervals = '1'
    num_sectors_per_side = '4 4 4 4 4 4'
    ring_block_ids = '4'
    ring_block_names = 'poison'
    background_block_ids = '1'
    background_block_names = 'graphite'
    background_intervals = 2
    quad_center_elements = true
  []
  [graphite_pin]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    polygon_size = ${fparse fuel_to_coolant_distance / 2.0}
    ring_radii = '${fparse compact_diameter / 2.0}'
    ring_intervals = '1'
    num_sectors_per_side = '4 4 4 4 4 4'
    ring_block_ids = '1'
    ring_block_names = 'graphite'
    background_block_ids = '1'
    background_block_names = 'graphite'
    quad_center_elements = true
  []
  [bundle]
    type = PatternedHexMeshGenerator
    inputs = 'fuel_pin coolant_pin poison_pin graphite_pin'
    hexagon_size = ${fparse bundle_flat_to_flat / 2.0 + bundle_gap_width / 2.0}
    pattern = '2 0 1 0 0 1 0 0 1 0 2;
              0 1 0 0 1 0 0 1 0 0 1 0;
             1 0 0 1 0 0 1 0 0 1 0 0 1;
            0 0 1 0 0 1 0 0 1 0 0 1 0 0;
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0;
          1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1;
         0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0;
        0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0;
       1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1;
      0 0 1 0 0 1 0 0 1 3 3 1 0 0 1 0 0 1 0 0;
     2 1 0 0 1 0 0 1 0 3 3 3 0 1 0 0 1 0 0 1 2;
      0 0 1 0 0 1 0 0 1 3 3 1 0 0 1 0 0 1 0 0;
       1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1;
        0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0;
         0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0;
          1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1;
           0 1 0 0 1 0 0 1 0 0 1 0 0 1 0;
            0 0 1 0 0 1 0 0 1 0 0 1 0 0;
             1 0 0 1 0 0 1 0 0 1 0 0 1;
              0 1 0 0 1 0 0 1 0 0 1 0;
               2 0 1 0 0 1 0 0 1 0 2'

    background_intervals = 2
    background_block_id = '1'
    background_block_names = 'graphite'
  []
  [graphite_background]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    polygon_size = ${fparse fuel_to_coolant_distance / 2.0}
    num_sectors_per_side = '4 4 4 4 4 4'
    background_block_ids = '5'
    background_block_names = 'graphite_centers'
    quad_center_elements = false
  []
  [graphite_bundle]
    type = PatternedHexMeshGenerator
    inputs = 'graphite_background'
    hexagon_size = ${fparse bundle_flat_to_flat / 2.0 + bundle_gap_width / 2.0}
    pattern = '0 0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0 0 0 0'

    background_intervals = 1
    background_block_id = '1'
    background_block_names = 'graphite'
  []
  [core]
    type = PatternedHexMeshGenerator
    inputs = 'bundle graphite_bundle'
    pattern = '1 0 1;
              0 0 0 0;
             1 0 1 0 1;
              0 0 0 0;
               1 0 1'
    generate_core_metadata = true
    pattern_boundary = none
  []
  [delete_coolant]
    type = BlockDeletionGenerator
    input = core
    block = '101'
  []

  [excore]
    type = PeripheralRingMeshGenerator
    input = delete_coolant
    peripheral_layer_num = 10
    peripheral_ring_radius = ${fparse vessel_inner_diameter / 2.0}
    input_mesh_external_boundary = 10000
    peripheral_ring_block_id = 0
    peripheral_ring_block_name = 'reflector'
    #growth_factor = 1.4
    peripheral_radial_bias = 1.4
  []

  [extrude]
    type = FancyExtruderGenerator
    input = excore
    heights = ${height}
    num_layers = ${num_layers}
    direction = '0 0 1'
  []
  [rename_coolant_sideset]
    type = RenameBoundaryGenerator
    input = extrude
    old_boundary = 102
    new_boundary = 'fluid_solid_interface'
  []
  [rename_block_ids]
    type = RenameBlockGenerator
    input = rename_coolant_sideset
    old_block = '0 1 2 4 5'
    new_block = 'reflector graphite compacts poison graphite_centers'
  []
  [rename_outer_sideset]
    type = RenameBoundaryGenerator
    input = rename_block_ids
    old_boundary = 10000
    new_boundary = 'vessel_outer'
  []

  construct_side_list_from_node_list = true
[]
