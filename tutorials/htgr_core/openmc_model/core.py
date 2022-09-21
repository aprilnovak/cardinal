#!/bin/env python

from argparse import ArgumentParser
import math
import numpy as np
import matplotlib.pyplot as plt

import openmc
import sys
import os

# Get common input parameters shared by other physics
script_dir = os.path.dirname(__file__)
sys.path.append(script_dir)
import common_input as specs
import materials as mats

def coolant_temp(t_in, t_out, l, z):
    """
    Computes the coolant temperature based on an expected cosine power distribution
    for a specified temperature rise. The total core temperature rise is governed
    by energy conservation as dT = Q / m / Cp, where dT is the total core temperature
    rise, Q is the total core power, m is the mass flowrate, and Cp is the fluid
    isobaric specific heat. If you neglect axial heat conduction and assume steady
    state, then the temperature rise in a layer of fluid i can be related to the
    ratio of the power in that layer to the total power,
    dT_i / dT = Q_i / Q. We assume here a sinusoidal power distribution to get
    a reasonable estimate of an initial coolant temperature distribution.

    Parameters
    ----------

    t_in : float
        Inlet temperature of the channel
    t_out : float
        Outlet temperature of the channel
    l : float
        Length of the channel
    z : float or 1-D numpy.array
        Axial position where the temperature will be computed

    Returns
    -------
        float or 1-D numpy array of float depending on z
    """
    dT = t_out - t_in
    Q = 2 * l / math.pi
    Qi = (l - l * np.cos(math.pi * z / l)) / math.pi

    t = t_in + Qi / Q * dT

    return t

def coolant_density(t):
  """
  Computes the helium density from temperature assuming a fixed operating pressure.

  Parameters
  ----------

  t : float
    Fluid temperature

  Returns
  _______
    float or 1-D numpy array of float depending on t
  """

  p_in_bar = specs.outlet_P * 1.0e-5;
  return 48.14 * p_in_bar / (t + 0.4446 * p_in_bar / math.pow(t, 0.2));

# -------------- Unit Conversions: OpenMC requires cm -----------
m = 100.0

# estimate the outlet temperature using bulk energy conservation for steady state
coolant_outlet_temp = specs.power / specs.mdot / specs.fluid_Cp + specs.inlet_T

# geometry parameters
coolant_channel_diam = specs.channel_diameter * m
reactor_bottom = 0.0
reactor_height = specs.height * m
reactor_top = reactor_bottom + reactor_height
bundle_pitch = specs.bundle_flat_to_flat * m + specs.bundle_gap_width * m
cell_pitch = specs.fuel_to_coolant_distance * m
fuel_channel_diam = specs.compact_diameter * m
top_reflector_height = specs.top_reflector_thickness * m
bottom_reflector_height = specs.bottom_reflector_thickness * m
vessel_radius = specs.vessel_inner_diameter / 2.0 * m
total_height = reactor_height + top_reflector_height + bottom_reflector_height

def assembly(n_ax_zones, n_rings, n_azimuthal, n_inactive, n_active, add_entropy_mesh=False):
    axial_section_height = reactor_height / n_ax_zones

    # superimposed search lattice
    triso_lattice_shape = (3, 3, int(axial_section_height * 2))

    model = openmc.model.Model()

    m_b4c = openmc.Material(name='B4C')
    enrichment_10 = specs.B10_enrichment
    mass_10 = openmc.data.atomic_mass('B10')
    mass_11 = openmc.data.atomic_mass('B11')

    # number of atoms in one gram of boron mixture
    n_10 = enrichment_10 / mass_10
    n_11 = (1.0 - enrichment_10) / mass_11
    total_n = n_10 + n_11
    grams_10 = n_10 / total_n
    grams_11 = n_11 / total_n

    # now, figure out how much carbon needs to be in the poison to get
    # an overall specified B10 weight percent
    total_b10_weight_percent = specs.total_B10_wt_percent
    total_mass = grams_10 / total_b10_weight_percent
    carbon_mass = total_mass - grams_10 - grams_11

    m_b4c.add_nuclide('B10', grams_10 / total_mass, 'wo')
    m_b4c.add_nuclide('B11', grams_11 / total_mass, 'wo')
    m_b4c.add_element('C', carbon_mass / total_mass, 'wo')
    m_b4c.set_density('kg/m3', specs.B4C_density)

    m_b4c_high = openmc.Material(name='B4C')
    enrichment_10 = specs.B10_enrichment

    # number of atoms in one gram of boron mixture
    n_10 = enrichment_10 / mass_10
    n_11 = (1.0 - enrichment_10) / mass_11
    total_n = n_10 + n_11
    grams_10 = n_10 / total_n
    grams_11 = n_11 / total_n

    # now, figure out how much carbon needs to be in the poison to get
    # an overall specified B10 weight percent
    total_b10_weight_percent = specs.total_B10_wt_percent
    total_mass = grams_10 / total_b10_weight_percent
    carbon_mass = total_mass - grams_10 - grams_11

    m_b4c_high.add_nuclide('B10', grams_10 / total_mass, 'wo')
    m_b4c_high.add_nuclide('B11', grams_11 / total_mass, 'wo')
    m_b4c_high.add_element('C', carbon_mass / total_mass, 'wo')
    m_b4c_high.set_density('kg/m3', specs.B4C_density)

    # reflector is 40 percent helium by volume (arbitrary assumption), with helium
    # at the inlet conditions
    rho = coolant_density(specs.inlet_T)
    matrix_density = 1700
    reflector_porosity = 0.40
    n_helium = reflector_porosity * rho / 4.002602
    n_carbon = (1.0 - reflector_porosity) * matrix_density / 12.0107
    combined_density = rho * reflector_porosity + matrix_density * (1.0 - reflector_porosity)
    m_reflector = openmc.Material(name='reflector')
    m_reflector.add_element('He', n_helium / (n_helium + n_carbon))
    m_reflector.add_element('C', n_carbon / (n_helium + n_carbon))
    m_reflector.set_density('kg/m3', combined_density)

    # TRISO particle
    radius_pyc_outer   = specs.oPyC_radius * m

    s_fuel             = openmc.Sphere(r=specs.kernel_radius*m)
    s_c_buffer         = openmc.Sphere(r=specs.buffer_radius*m)
    s_pyc_inner        = openmc.Sphere(r=specs.iPyC_radius*m)
    s_sic              = openmc.Sphere(r=specs.SiC_radius*m)
    s_pyc_outer        = openmc.Sphere(r=radius_pyc_outer)
    c_triso_fuel       = openmc.Cell(name='c_triso_fuel'     , fill=mats.m_fuel,              region=-s_fuel)
    c_triso_c_buffer   = openmc.Cell(name='c_triso_c_buffer' , fill=mats.m_graphite_c_buffer, region=+s_fuel      & -s_c_buffer)
    c_triso_pyc_inner  = openmc.Cell(name='c_triso_pyc_inner', fill=mats.m_graphite_pyc,      region=+s_c_buffer  & -s_pyc_inner)
    c_triso_sic        = openmc.Cell(name='c_triso_sic'      , fill=mats.m_sic,               region=+s_pyc_inner & -s_sic)
    c_triso_pyc_outer  = openmc.Cell(name='c_triso_pyc_outer', fill=mats.m_graphite_pyc,      region=+s_sic       & -s_pyc_outer)
    c_triso_matrix     = openmc.Cell(name='c_triso_matrix'   , fill=mats.m_graphite_matrix,   region=+s_pyc_outer)
    u_triso            = openmc.Universe(cells=[c_triso_fuel, c_triso_c_buffer, c_triso_pyc_inner, c_triso_sic, c_triso_pyc_outer, c_triso_matrix])

    # Channel surfaces
    fuel_cyl = openmc.ZCylinder(r=0.5 * fuel_channel_diam)
    coolant_cyl = openmc.ZCylinder(r=0.5 * coolant_channel_diam)
    poison_cyl = openmc.ZCylinder(r=0.5 * fuel_channel_diam)
    graphite_cyl = openmc.ZCylinder(r=0.5 * fuel_channel_diam)

    # create a TRISO lattice for one axial section (to be used in the rest of the axial zones)
    # center the TRISO region on the origin so it fills lattice cells appropriately
    min_z = openmc.ZPlane(z0=-0.5 * axial_section_height)
    max_z = openmc.ZPlane(z0=0.5 * axial_section_height)

    # region in which TRISOs are generated
    r_triso = -fuel_cyl & +min_z & -max_z

    rand_spheres = openmc.model.pack_spheres(radius=radius_pyc_outer, region=r_triso, pf=specs.triso_pf)
    random_trisos = [openmc.model.TRISO(radius_pyc_outer, u_triso, i) for i in rand_spheres]

    # Make the surfaces that will be used to cut in the radial direction
    r = 0.5 * fuel_channel_diam
    area = r * r
    radii = []
    r_inner = 0.0

    # Whether to create the special case with non-equal sized elements
    special_case = True

    if (special_case and n_rings != 2):
      raise ValueError("Special case is only implemented for 2 rings!")

    if (special_case):
      nr = 4
      ring_area = area / nr
      for i in range(nr):
        r = math.sqrt(ring_area + r_inner**2)

        if (i > 1):
          radii.append(r)

        r_inner = r
    else:
      ring_area = area / n_rings
      for i in range(n_rings):
        r = math.sqrt(ring_area + r_inner**2)
        radii.append(r)
        r_inner = r

    if (n_rings == 1):
      fuel_rings = [fuel_cyl]
    else:
      fuel_rings = [openmc.ZCylinder(r=r0) for r0 in radii]

    llc, urc = r_triso.bounding_box
    pitch = (urc - llc) / triso_lattice_shape

    # insert TRISOs into a lattice to accelerate point location queries
    triso_lattice = openmc.model.create_triso_lattice(random_trisos, llc, pitch, triso_lattice_shape, mats.m_graphite_matrix)

    axial_coords = np.linspace(reactor_bottom, reactor_top, n_ax_zones + 1)
    lattice_univs = []

    planes = []
    wedges = []
    n_wedges = n_azimuthal + int(n_azimuthal / 2)

    if (n_azimuthal > 1):
      n_planes = n_azimuthal + 1 + int(n_azimuthal / 2)
      if (n_azimuthal > 0):
        for i in range(n_planes):
          theta = i * math.pi / (n_azimuthal + n_azimuthal / 2)
          slope = math.tan(theta)

          if (i == n_planes - 1):
            p = planes[0]
          elif (abs(theta - math.pi / 2.0) < 1e-3):
            p = openmc.XPlane()
          elif (slope < 0.0):
            p = openmc.Plane(a=-slope, b=1.0, c=0.0, d=0.0)
          else:
            p = openmc.Plane(a=slope, b=-1.0, c=0.0, d=0.0)
          planes.append(p)

      for i in range(n_wedges):
        if (i != n_wedges - 1):
          wedges.append(+planes[i] & -planes[i + 1])
          wedges.append(+planes[i + 1] & -planes[i])
        else:
          wedges.append(-planes[i] & -planes[0])
          wedges.append(+planes[0] & +planes[i])

    if (n_azimuthal == 1):
      planes.append(openmc.XPlane())
      planes.append(openmc.Plane(a=-1.0 / math.sqrt(3.0), b=1.0, c=0.0, d=0.0))
      planes.append(openmc.Plane(a=1.0 / math.sqrt(3.0), b=1.0, c=0.0, d=0.0))

      wedges.append(+planes[0] & -planes[1])
      wedges.append(+planes[1] & +planes[2])
      wedges.append(-planes[2] & -planes[0])

    triso_lattice_cell = openmc.Cell(fill=triso_lattice)
    triso_lattice_univ = openmc.Universe(cells=[triso_lattice_cell])

    # New approach:
    cells = []
    if (n_azimuthal > 0):
      for wedge in wedges:
        for i, annulus in enumerate(openmc.model.subdivide(fuel_rings)[:-1]):
          c = openmc.Cell(region=annulus & wedge, fill=triso_lattice_univ)
          c.temperature = 500.0
          cells.append(c)
    else:
      print('No wedges')
      for i, annulus in enumerate(openmc.model.subdivide(fuel_rings)[:-1]):
        c = openmc.Cell(region=annulus, fill=triso_lattice_univ)
        c.temperature = 500.0
        cells.append(c)

    fuel_ch_matrix_cell = openmc.Cell(region=+fuel_cyl, fill=mats.m_graphite_matrix)
    fuel_ch_matrix_cell.temperature = 500.0
    cells.append(fuel_ch_matrix_cell)
    fuel_u = openmc.Universe(cells=cells)

    half_lattice_univs = []
    rotated_half_lattice_univs = []
    homogeneous_lattice_univs = []
    reflector_lattice_univs = []
    empty_lattice_univs = []

    m_colors = {}

    for z_min, z_max in zip(axial_coords[0:-1], axial_coords[1:]):
        # use the middle of the axial section to compute the temperature and density
        ax_pos = 0.5 * (z_min + z_max)
        T = coolant_temp(specs.inlet_T, coolant_outlet_temp, reactor_height, ax_pos)

        # create solid cells, which don't require us to clone materials in order to set temperatures
        poison_cell = openmc.Cell(region=-poison_cyl, fill=m_b4c)
        poison_cell.temperature = T
        poison_cell_high = openmc.Cell(region=-poison_cyl, fill=m_b4c_high)
        poison_cell_high.temperature = T
        poison_cell_matrix = openmc.Cell(region=+poison_cyl, fill=mats.m_graphite_matrix)
        poison_cell_matrix.temperature = T

        poison_matrix_cell = openmc.Cell(region=+poison_cyl, fill=mats.m_graphite_matrix)
        poison_matrix_cell.temperature = T

        graphite_cell = openmc.Cell(fill=mats.m_graphite_matrix)
        graphite_cell_g = openmc.Cell(fill=mats.m_graphite_matrix)
        graphite_cell.temperature = T
        graphite_cell_g.temperature = T

        coolant_matrix_cell = openmc.Cell(region=+coolant_cyl, fill=mats.m_graphite_matrix)
        coolant_matrix_cell.temperature = T

        void_cell = openmc.Cell(fill=None)

        # create fluid cells, which require us to clone the material in order to be able to
        # set unique densities
        coolant_cell = openmc.Cell(region=-coolant_cyl, fill=mats.m_coolant)
        n_coolant_cells_in_half_block = int(specs.n_coolant_channels_per_block / 2 + 5)
        n_coolant_in_diagonal_half_block = int(specs.n_coolant_channels_per_block / 2 + 3)
        #n_cools = specs.n_coolant_channels_per_block + n_coolant_cells_in_half_block * 2
        n_cools = n_coolant_cells_in_half_block + n_coolant_in_diagonal_half_block
        coolant_cell.fill = [mats.m_coolant.clone() for i in range(n_cools)]

        for mat in range(len(coolant_cell.fill)):
          m_colors[coolant_cell.fill[mat]] = 'red'

        # Define a universe for each type of solid-only pin (fuel, poison, and graphite)
        f = fuel_u
        p = openmc.Universe(cells=[poison_cell, poison_matrix_cell])
        pp = openmc.Universe(cells=[poison_cell_high, poison_cell_matrix])
        c = openmc.Universe(cells=[coolant_cell, coolant_matrix_cell])
        g = openmc.Universe(cells=[graphite_cell])
        gg = openmc.Universe(cells=[graphite_cell_g])
        v = openmc.Universe(cells=[void_cell])

        d = [f] * 2

        #ring2 = ([f] + [c]) * 6
        #ring3 = ([c] + d) * 6
        #ring4 = (d + [c] + [f]) * 6
        #ring5 = ([f] + [c] + d + [c]) * 6
        #ring6 = ([c] + d + [c] + d) * 6
        #ring7 = (d + [c] + d + [c] + [f]) * 6
        #ring8 = ([f] + [c] + d + [c] + d + [c]) * 6
        #ring9 = ([c] + d + [c] + d + [c] + d) * 6
        #ring10 = ([p] + [f] + [c] + d + [c] + d + [c] + [f]) * 6
        #ring11 = [g] * 66

        # inner two rings where there aren't any fuel/compact/poison pins
        ring1 = [g] * 6
        ring0 = [g]

        ring2 = [v] * 4 + ([f] + [c]) * 3 + [f] + [v]
        ring3 = [v] * 6 + ([c] + d) * 3 + [c] + [v] * 2
        ring4 = [v] * 8 + (d + [c] + [f]) * 3 + [f] + [v] * 3
        ring5 = [v] * 10 + ([f] + [c] + d + [c]) * 3 + [f] + [v] * 4
        ring6 = [v] * 12 + ([c] + d + [c] + d) * 3 + [c] + [v] * 5
        ring7 = [v] * 14 + (d + [c] + d + [c] + [f]) * 3 + [f] + [v] * 6
        ring8 = [v] * 16 + ([f] + [c] + d + [c] + d + [c]) * 3 + [f] + [v] * 7
        ring9 = [v] * 18 + ([c] + d + [c] + d + [c] + d) * 3 + [c] + [v] * 8
        ring10 = [v] * 20 + ([p] + [f] + [c] + d + [c] + d + [c] + [f]) * 3 + [p] + [v] * 9
        ring11 = [v] * 22 + [g] * 11 * 3 + [g] + [v] * 10

        lattice_univs.append([ring11, ring10, ring9, ring8, ring7, ring6, ring5, ring4, ring3, ring2, ring1, ring0])

        homogeneous_lattice_univs.append([[g] * 66, [g] * 60, [g] * 54, [g] * 48, [g] * 42, [g] * 36, [g] * 30, [g] * 24, [g] * 18, [g] * 12, [g] * 6, [g]])
        reflector_lattice_univs.append([[gg]])
        empty_lattice_univs.append([[v]])

        # now make a "half" assembly
        side = [f] + [c]
        hring2 = side * 2 + [v] * 5 + [c] + side

        side = [c] + d
        hring3 = side + [c] + [f] + [v] * 9 + [f] + side

        side = d + [c] + [f]
        hring4 = side + [f] * 2 + [c] + [v] * 11 + [c] + [f] + side

        side = [f] + [c] + d + [c]
        hring5 = side + [f] + [c] + [f] + [v] * 15 + [f] + [c] + side

        side = [c] + d + [c] + d
        hring6 = side + [c] + [f] * 2 + [c] + [v] * 17 + [c] + [f] * 2 + side

        side = d + [c] + d + [c] + [f]
        hring7 = side + d + [c] + [f] + [v] * 21 + [f] + [c] + [f] + side

        side = [f] + [c] + d + [c] + d + [c]
        hring8 = side + [f] + [c] + d + [c] + [v] * 23 + [c] + d + [c] + side

        side = [c] + d + [c] + d + [c] + d
        hring9 = side + [c] + d + [c] + [f] + [v] * 27 + [f] + [c] + d + side

        side = [p] + [f] + [c] + d + [c] + d + [c] + [f]
        hring10 = side + [pp] + [f] + [c] + d + [c] + [v] * 29 + [c] + d + [c] + [f] + side
        ring11 = [g] * 66

        half_lattice_univs.append([ring11, hring10, hring9, hring8, hring7, hring6, hring5, hring4, hring3, hring2, ring1, ring0])

        # now make another "half" assembly, but rotated
        side = [f] + [c]
        fring2 = [v] * 5 + [c] + side * 3

        side = [c] + d
        fring3 = [v] * 8 + [f] + side * 2 + [c] + [f] + [v]

        side = d + [c] + [f]
        fring4 = [v] * 10 + [c] + [f] + side * 2 + [f] * 2 + [c] + [v]

        side = [f] + [c] + d + [c]
        fring5 = [v] * 13 + [f] + [c] + side * 2 + [f] + [c] + [f] + [v] * 2

        side = [c] + d + [c] + d
        fring6 = [v] * 15 + [c] + [f] * 2 + side * 2 + [c] + [f] * 2 + [c] + [v] * 2

        side = d + [c] + d + [c] + [f]
        fring7 = [v] * 18 + [f] + [c] + [f] + side * 2 + d + [c] + [f] + [v] * 3

        side = [f] + [c] + d + [c] + d + [c]
        fring8 = [v] * 20 + [c] + d + [c] + side * 2 + [f] + [c] + d + [c] + [v] * 3

        side = [c] + d + [c] + d + [c] + d
        fring9 = [v] * 23 + [f] + [c] + d + side * 2 + [c] + d + [c] + [f] + [v] * 4

        side = [f] + [c] + d + [c] + d + [c] + [f]
        fring10 = [v] * 25 + [c] + d + [c] + [f] + [pp] + side + [p] + side + [p] + [f] + [c] + d + [c] + [v] * 4

        ring11 = [g] * 66
        rotated_half_lattice_univs.append([ring11, fring10, fring9, fring8, fring7, fring6, fring5, fring4, fring3, fring2, ring1, ring0])

    # create a hexagonal lattice used in each axial zone to represent the fuel
    hex_lattice = openmc.HexLattice(name="Bundle cell lattice")
    hex_lattice.orientation = 'x'
    hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    hex_lattice.pitch = (cell_pitch, axial_section_height)
    hex_lattice.universes = lattice_univs

    # create a "half" hexagonal lattice used in each axial zone to represent the fuel
    half_hex_lattice = openmc.HexLattice(name="Half bundle cell lattice")
    half_hex_lattice.orientation = 'x'
    half_hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    half_hex_lattice.pitch = (cell_pitch, axial_section_height)
    half_hex_lattice.universes = half_lattice_univs

    # create a rotated "half" hexagonal lattice used in each axial zone to represent the fuel
    #rotated_half_hex_lattice = openmc.HexLattice(name="Rotated Half bundle cell lattice")
    #rotated_half_hex_lattice.orientation = 'x'
    #rotated_half_hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    #rotated_half_hex_lattice.pitch = (cell_pitch, axial_section_height)
    #rotated_half_hex_lattice.universes = rotated_half_lattice_univs

    # create a hexagonal lattice used in each axial zone to represent the homogeneous regions
    homogeneous_hex_lattice = openmc.HexLattice(name="Homogeneous cell lattice")
    homogeneous_hex_lattice.orientation = 'x'
    homogeneous_hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    homogeneous_hex_lattice.pitch = (cell_pitch, axial_section_height)
    homogeneous_hex_lattice.universes = homogeneous_lattice_univs

    # create a single-pin lattice in each axial zone to represent the homogeneous regions
    reflector_hex_lattice = openmc.HexLattice(name="Reflector cell lattice")
    reflector_hex_lattice.orientation = 'y'
    reflector_hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    reflector_hex_lattice.pitch = (bundle_pitch, axial_section_height)
    reflector_hex_lattice.universes = reflector_lattice_univs

    # create a void hexagonal lattice
    empty_hex_lattice = openmc.HexLattice(name="Empty cell lattice")
    empty_hex_lattice.orientation = 'y'
    empty_hex_lattice.center = (0.0, 0.0, 0.5 * (reactor_bottom + reactor_top))
    empty_hex_lattice.pitch = (bundle_pitch, axial_section_height)
    empty_hex_lattice.universes = empty_lattice_univs

    hexagon_volume = reactor_height * math.sqrt(3) / 2.0 * bundle_pitch**2
    coolant_channel_volume = math.pi * coolant_channel_diam**2 / 4.0 * reactor_height

    graphite_outer_cell = openmc.Cell(fill=mats.m_graphite_matrix)
    graphite_outer_cell.temperature = T
    inf_graphite_univ = openmc.Universe(cells=[graphite_outer_cell])
    hex_lattice.outer = inf_graphite_univ
    half_hex_lattice.outer = inf_graphite_univ
    #rotated_half_hex_lattice.outer = inf_graphite_univ
    homogeneous_hex_lattice.outer = inf_graphite_univ
    reflector_hex_lattice.outer = inf_graphite_univ
    empty_hex_lattice.outer = inf_graphite_univ

    fuel_bundle_wall = openmc.hexagonal_prism(bundle_pitch / math.sqrt(3.0), 'x')
    fuel_bundle = openmc.Cell(region=fuel_bundle_wall, fill=hex_lattice)
    fuel_bundle_univ = openmc.Universe(cells=[fuel_bundle])

    half_fuel_bundle = openmc.Cell(region=fuel_bundle_wall, fill=half_hex_lattice)
    half_fuel_bundle_univ = openmc.Universe(cells=[half_fuel_bundle])

    #rotated_half_fuel_bundle = openmc.Cell(region=fuel_bundle_wall, fill=rotated_half_hex_lattice)
    #rotated_half_fuel_bundle_univ = openmc.Universe(cells=[rotated_half_fuel_bundle])

    homogeneous_bundle = openmc.Cell(region=fuel_bundle_wall, fill=homogeneous_hex_lattice)
    homogeneous_bundle_univ = openmc.Universe(cells=[homogeneous_bundle])

    homogeneous_reflector = openmc.Cell(region=fuel_bundle_wall, fill=reflector_hex_lattice)
    homogeneous_reflector_univ = openmc.Universe(cells=[homogeneous_reflector])

    empty = openmc.Cell(region=fuel_bundle_wall, fill=empty_hex_lattice)
    empty_univ = openmc.Universe(cells=[empty])

    # create a lattice of fuel bundles
    core_lattice = openmc.HexLattice(name="Core lattice")
    core_lattice.orientation = 'y'
    core_lattice.center = (0.0, 0.0)
    core_lattice.pitch = [bundle_pitch]

    ring0 = [homogeneous_bundle_univ]
    ring1 = [half_fuel_bundle_univ] + [empty_univ] + [empty_univ] * 4
    ring2 = [homogeneous_reflector_univ] + [fuel_bundle_univ] + \
            [homogeneous_reflector_univ] + [empty_univ] * 9

    core_lattice.universes = [[homogeneous_reflector_univ] * 30, [homogeneous_reflector_univ] * 24, \
        [homogeneous_reflector_univ] * 18, ring2, ring1, ring0]
    core_lattice.outer = inf_graphite_univ

    # create additional axial regions
    axial_planes = [openmc.ZPlane(z0=coord) for coord in axial_coords]

    # axial planes
    min_z = axial_planes[0]
    max_z = axial_planes[-1]

    # add symmetry plane
    symmetry_plane = openmc.XPlane(boundary_type='reflective')
    second_symmetry_plane = openmc.Plane(a=-math.sqrt(3.0) / 2.0, b=0.5, c=0.0, d=0.0, boundary_type='reflective')

    # fill the unit cell with the hex lattice
    vessel_inner_surface = openmc.ZCylinder(r=vessel_radius, boundary_type='vacuum')
    outer_cell = openmc.Cell(region=-vessel_inner_surface & +min_z & -max_z & +symmetry_plane & +second_symmetry_plane, fill=core_lattice)

    # add the top and bottom reflector
    top_refl_z = reactor_height + top_reflector_height
    top_refl = openmc.ZPlane(z0=top_refl_z, boundary_type='vacuum')
    bottom_refl_z = -bottom_reflector_height
    bottom_refl = openmc.ZPlane(z0=bottom_refl_z, boundary_type='vacuum')

    top_refl_cell = openmc.Cell(region=-vessel_inner_surface & +max_z & -top_refl & +symmetry_plane & +second_symmetry_plane, fill=m_reflector)
    bottom_refl_cell = openmc.Cell(region=-vessel_inner_surface & -min_z & +bottom_refl & +symmetry_plane & +second_symmetry_plane, fill=m_reflector)
    top_refl_cell.temperature = specs.inlet_T
    bottom_refl_cell.temperature = coolant_outlet_temp

    model.geometry = openmc.Geometry([outer_cell, top_refl_cell, bottom_refl_cell])

    ### Settings ###
    settings = openmc.Settings()

    settings.particles = 150000
    settings.inactive = n_inactive
    settings.batches = settings.inactive + n_active
    settings.temperature['method'] = 'interpolation'
    settings.temperature['multipole'] = True
    settings.temperature['range'] = (294.0, 1500.0)

    l = 2 * bundle_pitch / math.sqrt(3.0)
    lower_left = (0.0, 0.0, reactor_bottom)
    upper_right = (l,  l, reactor_top)
    source_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    source = openmc.Source(space=source_dist)
    settings.source = source

    if (add_entropy_mesh):
        entropy_mesh = openmc.RegularMesh()
        entropy_mesh.lower_left = lower_left
        entropy_mesh.upper_right = upper_right
        entropy_mesh.dimension = (10, 10, 20)
        settings.entropy_mesh = entropy_mesh

    model.settings = settings

    m_colors[mats.m_fuel] = 'palegreen'
    m_colors[mats.m_graphite_c_buffer] = 'sandybrown'
    m_colors[mats.m_graphite_pyc] = 'orange'
    m_colors[mats.m_sic] = 'yellow'
    m_colors[mats.m_graphite_matrix] = 'darkblue'
    m_colors[m_b4c] = 'lightskyblue'
    m_colors[m_b4c_high] = 'black'
    m_colors[mats.m_graphite_reflector] = 'red'

    bundle_p_rounded = int(bundle_pitch)

    plot1          = openmc.Plot()
    plot1.filename = 'plot1'
    plot1.width    = (2.1 * vessel_radius, 2.1 * total_height)
    plot1.basis    = 'xz'
    plot1.origin   = (0.0, 0.0, reactor_height / 2.0)
    plot1.pixels   = (100 * 2 * bundle_p_rounded, int(100 * 3 * axial_section_height))
    plot1.color_by = 'cell'

    plot2          = openmc.Plot()
    plot2.filename = 'plot2'
    plot2.width    = (2.1 * vessel_radius, 2.1 * vessel_radius)
    plot2.basis    = 'xy'
    plot2.origin   = (0.0, 0.0, axial_section_height / 4.0)
    plot2.pixels   = (250 * bundle_p_rounded, 250 * bundle_p_rounded)
    plot2.color_by = 'material'
    plot2.colors   = m_colors

    plot3          = openmc.Plot()
    plot3.filename = 'plot3'
    plot3.width    = plot2.width
    plot3.basis    = plot2.basis
    plot3.origin   = plot2.origin
    plot3.pixels   = (500 * bundle_p_rounded, 500 * bundle_p_rounded)
    plot3.color_by = 'cell'
    plot3.level = 2

    model.plots = openmc.Plots([plot1, plot2, plot3])

    return model


def main():

    ap = ArgumentParser()
    ap.add_argument('-n', dest='n_axial', type=int, default=50,
                    help='Number of axial cell divisions')
    ap.add_argument('-s', '--entropy', action='store_true',
                    help='Whether to add a Shannon entropy mesh')
    ap.add_argument('-i', dest='n_inactive', type=int, default=100,
                    help='Number of inactive cycles')
    ap.add_argument('-a', dest='n_active', type=int, default=100,
                    help='Number of active cycles')
    ap.add_argument('-r', dest='n_rings', type=int, default=2,
                    help='Number of rings in fuel compacts')
    ap.add_argument('-t', dest='n_azimuthal', type=int, default=1,
                    help='Number of azimuthal divisions in fuel copmacts')

    args = ap.parse_args()

    model = assembly(args.n_axial, args.n_rings, args.n_azimuthal, args.n_inactive, args.n_active, args.entropy)
    model.export_to_xml()

if __name__ == "__main__":
    main()

