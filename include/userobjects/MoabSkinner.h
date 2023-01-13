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
  virtual void execute() override;
  virtual void finalize() override {};

  /// Helper method to create MOAB tags
  virtual void createTags();

  /**
   * Create a map from the [Mesh] DOFs to MOAB entity handles
   * @return mapping from DOFs to MOAB entity handles
   */
  virtual std::map<dof_id_type,moab::EntityHandle> createNodes();

  /// Helper method to create MOAB elements
  void createElems(std::map<dof_id_type,moab::EntityHandle>& node_id_to_handle);

  /**
   * Set length multiplier to convert from [Mesh] to OpenMC's centimeters
   * @param[in] s scaling factor
   */
  virtual void setScaling(const Real & s) { _scaling = s; }

  /// Clear the maps between entity handles and DOF IDs
  virtual void clearElemMaps();

  /** Return all sets of node indices; for second-order (TET10) elements, this
   * splits each tet into 8 sub-elements to be compatible with MOAB's first-order TET4.
   * @param[in] type element type
   * @return node indices for each tet
   */
  virtual std::vector<std::vector<unsigned int>> getTetSets(ElemType type);

  /**
   * Add an element to the libMesh-MOAB mapping
   * @param[in] id libMesh ID
   * @param[in] ent entity handle
   */
  virtual void addElem(dof_id_type id, moab::EntityHandle ent);

  /// Write the full MOAB mesh to _output_base_full_<int>.h5m
  virtual void writeFullMesh() const;

  /// Write the MOAB mesh
  virtual void write();

  /** Helper function to wrap moab::tag_set_data for a string
   * TODO: meaning of params?
   */
  void setTagData(moab::Tag tag, moab::EntityHandle ent, std::string data, unsigned int SIZE);

  /// Helper function to wrap moab::tag_set_data for a generic pointer
  void setTagData(moab::Tag tag, moab::EntityHandle ent, void* data);

  /// Generic method to set the tags that DAGMC requires
  void setTags(moab::EntityHandle ent, std::string name, std::string category, unsigned int id, int dim);

  /// Helper method to create MOAB group entity set
  void createGroup(unsigned int id, std::string name, moab::EntityHandle & group_set);

  /// Helper method to create MOAB volume entity set
  void createVol(unsigned int id, moab::EntityHandle & volume_set, moab::EntityHandle group_set);


protected:

  /// Encode the whether the surface normal faces into or out of the volume
  enum Sense { BACKWARDS=-1, FORWARDS=1};

  /// \brief Encode MOAB information about volumes needed when creating surfaces
  struct VolData
  {
    moab::EntityHandle vol;
    Sense sense;
  };

  /// Helper method to create MOAB surface entity set
  void createSurf(unsigned int id,moab::EntityHandle& surface_set, moab::Range& faces,  std::vector<VolData> & voldata);

  /// Add parent-child metadata relating a surface to its volume
  void updateSurfData(moab::EntityHandle surface_set, VolData data);

  /// Faceting tolerence needed by DAGMC
  const Real & _faceting_tol;

  /// Geometry tolerence needed by DAGMC
  const Real & _geom_tol;

  /// Whether to output the full MOAB mesh
  const bool & _output_full;

  /// Base filename to write full (entire) MOAB meshes
  const std::string & _output_base_full;

  /// File index for writing MOAB meshes
  unsigned int _n_write;

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

  /// Map from libmesh ID to MOAB element entity handles
  std::map<dof_id_type, std::vector<moab::EntityHandle>> _id_to_elem_handles;

  /// First tet entity handle
  moab::EntityHandle _offset;

  /// Save some topological data: map from surface handle to vol handle and sense
  std::map<moab::EntityHandle, std::vector<VolData>> _surfs_to_vols;

  /// Number of nodes per MOAB Tet
  const unsigned int NODES_PER_MOAB_TET = 4;
};
