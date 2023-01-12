#pragma once

#include "GeneralUserObject.h"

#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTagConventions.hpp"

class MoabSkinner : public GeneralUserObject
{
public:

  MoabSkinner(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void initialize() override;

  /// Helper method to create MOAB tags
  virtual void createTags();

  /**
   * Create a map from the [Mesh] DOFs to MOAB entity handles
   * @return mapping from DOFs to MOAB entity handles
   */
  virtual std::map<dof_id_type,moab::EntityHandle> createNodes();

  /**
   * Set length multiplier to convert from [Mesh] to OpenMC's centimeters
   * @param[in] s scaling factor
   */
  virtual void setScaling(const Real & s) { _scaling = s; }

protected:
  /// Faceting tolerence needed by DAGMC
  const Real & _faceting_tol;

  /// Geometry tolerence needed by DAGMC
  const Real & _geom_tol;

  /// MOAB interface
  std::shared_ptr<moab::Interface> _moab;

  /// Pointer to a moab skinner for finding temperature surfaces
  std::unique_ptr<moab::Skinner> _skinner;

  /// Pointer for topology for setting surface sense
  std::unique_ptr<moab::GeomTopoTool> _gtt;

  /// An entitiy handle to represent the set of all tets
  moab::EntityHandle _meshset;

  /// Tag for dimension for geometry
  moab::Tag _geometry_dimension_tag;

  /// Tag for entity set ID
  moab::Tag _id_tag;

  /// Tag for faceting tolerance
  moab::Tag _faceting_tol_tag;

  /// Tag needed by DAGMC
  moab::Tag _geometry_resabs_tag;

  /// Tag for type of entity set
  moab::Tag _category_tag;

  /// Tag for name of entity set
  moab::Tag _name_tag;

  /// Length multiplier to convert from [Mesh] to OpenMC's unit of centimeters
  Real _scaling;
};
