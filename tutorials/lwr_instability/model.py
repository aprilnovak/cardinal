import openmc
import numpy as np
import openmc.universe
import common as specs

model = openmc.Model()

## The radius of a fuel pin
r_fuel = specs.pellet_diameter / 2 * 1e2

## The outer diameter of the cladding
r_clad = specs.pin_diameter / 2 * 1e2

## The pitch of a single lattice element.
pitch = specs.pin_pitch * 1e2

## The height of the fuel assemblies from the axial midplane.
core_height = specs.heated_length * 1e2

# Some discretization parameters.
core_axial_slices = specs.n_cell_layers

# Material definitions.

## Fuel region: UO2 at 3% enriched.
uo2 = openmc.Material(name = 'UO2 Fuel')
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cm3', 11.0)

## Borated water at 1000 ppm
h2o = openmc.model.borated_water(1000, specs.T_inlet, specs.P_outlet * 1e-6, temp_unit='K', press_unit='MPa')

## Zr clad.
zr = openmc.Material(name = 'Zr Cladding')
zr.add_element('Zr', 1.0)
zr.set_density('g/cm3', 6.6)

# Geometry definitions.
### The entire UO2 pincell.
fuel_pin_or = openmc.ZCylinder(r = r_fuel)
fuel_zr_or = openmc.ZCylinder(r = r_clad)
fuel_bb = openmc.model.RectangularPrism(width = pitch, height = pitch)

zr_clad_cell_4 = openmc.Cell(name = 'UO2 Pin Zr Clad')
zr_clad_cell_4.region = +fuel_pin_or & -fuel_zr_or
zr_clad_cell_4.fill = zr

h2o_bb_cell_4 = openmc.Cell(name = 'UO2 Pin Water Bounding Box')
h2o_bb_cell_4.region = +fuel_zr_or & -fuel_bb
h2o_bb_cell_4.fill = h2o

uo2_fuel_cell = openmc.Cell(name = 'UO2 Fuel Pin')
uo2_fuel_cell.region = -fuel_pin_or
uo2_fuel_cell.fill = uo2
uo2_u = openmc.Universe(cells=[uo2_fuel_cell, zr_clad_cell_4, h2o_bb_cell_4])

uo2_assembly_cells = [
    [uo2_u, uo2_u],
    [uo2_u, uo2_u]
]

## The assembly.
assembly_bb = openmc.model.RectangularPrism(width=2*pitch,height=2*pitch,boundary_type='reflective')

uo2_assembly = openmc.RectLattice(name = 'UO2 Assembly')
uo2_assembly.pitch = (pitch, pitch)
uo2_assembly.lower_left = (-2 * pitch / 2.0, -2 * pitch / 2.0)
uo2_assembly.universes = uo2_assembly_cells

core_z_planes = []
for z in np.linspace(0.0, core_height, core_axial_slices + 1):
  core_z_planes.append(openmc.ZPlane(z0 = z))
core_z_planes[0].boundary_type = 'vacuum'
core_z_planes[-1].boundary_type = 'vacuum'

uo2_assembly_cells = []
all_cells = []
for i in range(core_axial_slices):
  uo2_assembly_cells.append(openmc.Cell(name = 'UO2 Assembly Cell ' + str(i), region = -assembly_bb & +core_z_planes[i] & -core_z_planes[i + 1], fill = uo2_assembly))
  all_cells.append(uo2_assembly_cells[-1])

# Setup the model.
model.geometry = openmc.Geometry(all_cells)

## The simulation settings.
model.settings.source = [openmc.IndependentSource(space = openmc.stats.Box(lower_left = (-2 * pitch / 2.0, -2 * pitch / 2.0, 0.0),
                                                                                upper_right = (2 * pitch / 2.0, 2 * pitch / 2.0, core_height)))]

model.settings.batches = 200
model.settings.inactive = 100
model.settings.particles = 20000

model.settings.ptables   = True
model.settings.temperature['method']='interpolation'
model.settings.temperature['range'] = (300.0, 3000.0)
model.settings.temperature['default'] = specs.T_inlet

model.export_to_model_xml()
