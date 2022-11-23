# Full-Core Coupling of OpenMC, MOOSE, and THM for an HTGR

This page contains all the input files used to create the full-core
[!ac](HTGR) model in:

A.J. Novak, D. Andrs, P. Shriwise, J. Fang, H. Yuan, D. Shaver, E. Merzari, P.K. Romano, and R.C. Martineau,
["Coupled Monte Carlo and Thermal-Fluid Modeling of High Temperature Gas Reactors Using Cardinal"](https://www.sciencedirect.com/science/article/pii/S0306454922003450) 177 (109310) *Annals of Nuclear Energy* (2022)

This model is a 12-bundle HTGR core surrounded by a graphite reflector. The files can be accessed by

```
cd tutorials/htgr_core
```

!media htgr_core.png
  id=htgr_core
  caption=Top-down view of the 12-bundle HTGR core geometry, colored by material
  style=width:80%;margin-left:auto;margin-right:auto

!alert! note
This model
is very large, and you will most likely run out of memory attempting to run it unless you
split the meshes with MOOSE's [distributed mesh](https://mooseframework.inl.gov/syntax/Mesh/splitting.html)
features. Even if you decrease the number of layers, based on how we build the
model in OpenMC (by repeating the same partial-height compact everywhere in the model),
this may not be enough to help you memory-wise (because you'd be storing *more* TRISO particle surfaces
in the repeated compact surface).
!alert-end!

## Meshes

The solid mesh (used for the OpenMC solid mesh mirror and the MOOSE heat conduction
solve) is generated with:

```
cd meshes
cardinal-opt -i common_input.i solid.i --mesh-only
```

This will generate the following solid mesh.

!media htgr_fullcore.png
  id=htgr_fullcore
  caption=HTGR full-core solid mesh
  style=width:80%;margin-left:auto;margin-right:auto

No special syntax is needed to generate the THM mesh, since THM automatically builds
its own (1-D) mesh based on the component syntax.

## THM Fluid Model

The THM model consists of a single HTGR flow channel in 1-D. In the solid input file
shown later in [#solid_htgr],
this model will be repeated for each of the 1296 coolant channels.

!listing tutorials/htgr_core/thm.i

## MOOSE Heat Conduction Model
  id=solid_htgr

The MOOSE heat conduction input file is as follows. Heat is produced in the fuel
compacts, while each coolant channel is represented with one THM sub-app.

!listing tutorials/htgr_core/solid_thm.i

## OpenMC Model

The OpenMC [!ac](CSG) model is generated with:

```
cd openmc_model
python core.py
```

This will generate the following [!ac](CSG) OpenMC geometry and the needed
XML files to run OpenMC.

!media htgr_openmc.png
  id=htgr_openmc
  caption=HTGR OpenMC [!ac](CSG) geometry
  style=width:80%;margin-left:auto;margin-right:auto

Then, we structure OpenMC as the main application as follows.

!listing tutorials/htgr_core/openmc_thm.i
