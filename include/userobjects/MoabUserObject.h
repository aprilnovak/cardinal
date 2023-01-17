#pragma once

#include "GeneralUserObject.h"
#include "MaterialBase.h"

#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTagConventions.hpp"
#include "libmesh/mesh_function.h"

// TODO: see what I can make const

/// Convenience struct
struct MOABMaterialProperties{
  double rel_density;
  double temp;
};

/**
 * \brief Skins the [Mesh] according to individual bins for temperature, density, and subdomain ID
 *
 * Skins a [Mesh] according to temperature, density, and subdomain. The surfaces bounding
 * those grouped elements are then generated, providing geometry information needed for DAGMC
 * to then track particles on this new geometry.
 */
class MoabUserObject : public GeneralUserObject
{
public:

  MoabUserObject(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void execute() override;

  virtual void initialize() override;

  virtual void finalize() override {};

  virtual void threadJoin(const UserObject & /*uo*/) override {};

  /// Perform the skinning operation
  virtual void update();

  /**
   * Whether the element is owned by this rank
   * @return whether element is owned by this rank
   */
  bool isLocalElem(const Elem * elem) const;

  /**
   * Get the bin index for the temperature
   * @param[in] elem element
   * @return temperature bin index
   */
  virtual unsigned int getTemperatureBin(const Elem * const elem) const;

  /**
   * Get the bin index for the density
   * @param[in] elem element
   * @return density bin index
   */
  virtual unsigned int getDensityBin(const Elem * const elem) const;

  /**
   * Get the bin index for the subdomain
   * @param[in] elem element
   * @return subdomain bin index
   */
  virtual unsigned int getSubdomainBin(const Elem * const elem) const { return _blocks.at(elem->subdomain_id()); }

  /**
   * Set the length multiplier to get from [Mesh] units into centimeters
   * @param[in] scale multiplier
   */
  virtual void setScaling(const Real & scale) { _scaling = scale; }

  /**
   * Set the verbosity level
   * @param[in] verbose whether to print diagnostic information
   */
  virtual void setVerbosity(const bool & verbose) { _verbose = verbose; }

  /**
   * Set whether the libMesh mesh is fixed
   * @param[in] fixed whether mesh is fixed
   */
  virtual void setFixedMesh(const bool & fixed) { _fixed_mesh = fixed; }

  /**
   * Indicate whether this userobject is run by itself (for testing purposes)
   * or controlled by some other class.
   */
  virtual void makeDependentOnExternalAction() { _standalone = false; }

  /**
   * Get variable number in the auxiliary system
   * @param[in] name variable name
   * @param[in] param_name parameter name, for printing a helpful error message
   * @return variable number
   */
  unsigned int getAuxiliaryVariableNumber(const std::string & name, const std::string & param_name) const;

  /// TODO Clear mesh data
  void reset();

  /// TODO Retrieve a list of original material names and properties
  void getMaterialProperties(std::vector<std::string>& mat_names_out,
                             std::vector<double>& initial_densities,
                             std::vector<std::string>& tails,
                             std::vector<MOABMaterialProperties>& properties);

  /// MOAB interface
  std::shared_ptr<moab::Interface> _moab;

  /**
   * Get total bin index given individual indices for the temperature, density, and subdomain bins
   * @param[in] temp_bin temperature bin
   * @param[in] density_bin density bin
   * @param[in] subdomain_bin subdomain ID bin
   * @return total bin index
   */
  virtual unsigned int getBin(const unsigned int & temp_bin, const unsigned int & density_bin, const unsigned int & subdomain_bin) const;

protected:
  std::unique_ptr<NumericVector<Number>> _serialized_solution;

  /// Whether to print diagnostic information
  bool _verbose;

  /// Whether the libMesh mesh is fixed in time (e.g. no spatial deformations or refinement)
  bool _fixed_mesh;

  /// Whether to build a graveyard as two additional cube surfaces surrounding the mesh.
  const bool & _build_graveyard;

  /// Name of the temperature variable
  const std::string & _temperature_name;

  /// Lower bound of temperature bins
  const Real & _temperature_min;

  /// Upper bound of temperature bins
  const Real & _temperature_max;

  /// Number of temperature bins
  const unsigned int & _n_temperature_bins;

  /// Temperature bin width
  const Real _temperature_bin_width;

  /// Whether elements are binned by density (in addition to temperature and block)
  const bool _bin_by_density;

  /// material names
  std::vector<std::string> mat_names;
  /// OpenMC material names
  std::vector<std::string> openmc_mat_names;

  /// Faceting tolerence needed by DAGMC
  const Real & _faceting_tol;

  /// Geometry tolerence needed by DAGMC
  const Real & _geom_tol;

  /// Multiplier on bounding box for inner surface of graveyard
  const Real & _graveyard_scale_inner;

  /// Multiplier on bounding box for outer surface of graveyard
  const Real & _graveyard_scale_outer;

  /// Whether to output the MOAB mesh skins to a .h5m file
  const bool & _output_skins;

  /// Whether to output the MOAB mesh to a .h5m file
  const bool & _output_full;

  /// Length multiplier to get from [Mesh] units into OpenMC's centimeters
  Real _scaling;

  /// Count number of times file has been written to
  unsigned int _n_write;

  /// Whether this class runs by itself, or is controlled by an external class
  bool _standalone;

  /// Encode the whether the surface normal faces into or out of the volume
  enum Sense { BACKWARDS=-1, FORWARDS=1};

  /// \brief Encode MOAB information about volumes needed when creating surfaces
  struct VolData{
    moab::EntityHandle vol;
    Sense sense;
  };

  /// Get a modifyable reference to the underlying libmesh mesh.
  MeshBase& mesh();

  /// Get a modifyable reference to the underlying libmesh equation systems
  EquationSystems & systems();

  /// Get a modifyable reference to the underlying libmesh system
  System& system(std::string var_now);

  /// Helper method to create MOAB elements
  void createMOABElems();

  /// Helper method to create MOAB tags
  virtual void createTags();

  /// Helper method to create MOAB group entity set
  moab::ErrorCode createGroup(unsigned int id, std::string name,moab::EntityHandle& group_set);

  /// Helper method to create MOAB volume entity set
  moab::ErrorCode createVol(unsigned int id,moab::EntityHandle& volume_set,moab::EntityHandle group_set);

  /// Helper method to create MOAB surface entity set
  void createSurf(unsigned int id,moab::EntityHandle& surface_set, moab::Range& faces,  std::vector<VolData> & voldata);

  /// Helper method to create MOAB surfaces with no overlaps
  moab::ErrorCode createSurfaces(moab::Range& reversed, VolData& voldata, unsigned int& surf_id);

  /**
   * Create a MOAB surface from a bounding box
   * @param[in] box bounding box
   * @param[in]
   * @param[in]
   * @param[in]
   * @param[in] factor
   */
  void createSurfaceFromBox(const BoundingBox& box, VolData& voldata, unsigned int& surf_id, bool normalout, const Real & factor);

  /**
   * Create MOAB nodes from a bounding box
   * @param[in] box bounding box
   * @param[in] factor multiplicative factor to resize the bounding box sides
   * @return nodes
   */
  std::vector<moab::EntityHandle> createNodesFromBox(const BoundingBox & box, const Real & factor) const;

  /// Create 3 tri faces stemming from one corner of a cude (an open tetrahedron)
  void createCornerTris(const std::vector<moab::EntityHandle> & verts,
                                   unsigned int corner, unsigned int v1,
                                   unsigned int v2 ,unsigned int v3,
                                   bool normalout, moab::Range &surface_tris);

  /// Create MOAB tri surface element
  moab::EntityHandle createTri(const std::vector<moab::EntityHandle> & vertices,unsigned int v1, unsigned int v2 ,unsigned int v3);

  /// Add parent-child metadata relating a surface to its volume
  moab::ErrorCode updateSurfData(moab::EntityHandle surface_set,VolData data);

  /// Generic method to set the tags that DAGMC requires
  moab::ErrorCode setTags(moab::EntityHandle ent,std::string name, std::string category, unsigned int id, int dim);

  /// Helper function to wrap moab::tag_set_data for a string
  moab::ErrorCode setTagData(moab::Tag tag, moab::EntityHandle ent, std::string data, unsigned int SIZE);

  /// Helper function to wrap moab::tag_set_data for a generic pointer
  moab::ErrorCode setTagData(moab::Tag tag, moab::EntityHandle ent, void* data);

  /**
   * Get the node numberings for the MOAB TET4 elements to build for each [Mesh] element
   * @param[in] type element type
   */
  const std::vector<std::vector<unsigned int>> & getTetSets(ElemType type) const;

  /// Build the graveyard (needed by OpenMC)
  void buildGraveyard(unsigned int & vol_id, unsigned int & surf_id);

  /// Store a mapping from [Mesh] subdomain IDs to an index, to be used for binning by block ID
  virtual void findBlocks();

  /// Helper method to convert between elem / solution indices
  dof_id_type elem_to_soln_index(const Elem& elem,unsigned int iSysNow, unsigned int iVarNow);

  /// Sort all the elements in the [Mesh] into bins for temperature, density, and subdomain.
  virtual void sortElemsByResults();

  /// Group the binned elems into local temperature regions and find their surfaces
  bool findSurfaces();

  /// Group a given bin into local regions
  /// NB elems in param is a copy, localElems is a reference
  void groupLocalElems(std::set<dof_id_type> elems, std::vector<moab::Range>& localElems);

  /// Map density and temp bin indices onto a linearised index
  int getMatBin(int iVarBin, int iDenBin, int n_temperature_binsIn, int _n_density_binsIn);
  /// Map density and temp bin indices onto a linearised index
  /// with default parameters for number of bins
  int getMatBin(int iVarBin, int iDenBin){
    return getMatBin(iVarBin,iDenBin,_n_temperature_bins,_n_density_bins);
  }

  /// Clear the containers of elements grouped into bins of constant temp
  void resetContainers();

  /// Clear MOAB entity sets
  bool resetMOAB();

  /// Find the surfaces for the provided range and add to group
  bool findSurface(const moab::Range& region,moab::EntityHandle group, unsigned int & vol_id, unsigned int & surf_id,moab::EntityHandle& volume_set);

  /// Write MOAB volume and/or skin meshes to file
  virtual void write();

  /// Moab skinner for finding temperature surfaces
  std::unique_ptr<moab::Skinner> skinner;

  /// Topology tool for setting surface sense
  std::unique_ptr<moab::GeomTopoTool> gtt;

  /// Map from libmesh id to MOAB element entity handles
  std::map<dof_id_type,std::vector<moab::EntityHandle> > _id_to_elem_handles;

  /// Save the first tet entity handle
  moab::EntityHandle offset;

  /// Name of the MOOSE variable containing the density
  std::string den_var_name;

  /// Lower bound of density bins
  Real _density_min;

  /// Upper bound of density bins
  Real _density_max;

  /// Density bin width
  Real _density_bin_width;

  /// Number of density bins
  unsigned int _n_density_bins;

  /// Number of block bins
  unsigned int _n_block_bins;

  /// Mapping from total bin ID to a set of elements sorted into that bin
  std::vector<std::set<dof_id_type>> _elem_bins;

  /// A place to store the entire solution
  // N.B. For big problems this is going to be a memory bottleneck
  // TODO: We will need to come up with a better solution
  // Map is from systems index
  std::map<unsigned int, std::unique_ptr<NumericVector<Number> > > serial_solutions;

  // Materials data

  /// Blocks in the [Mesh]
  std::map<SubdomainID, unsigned int> _blocks;

  /// Entity handle to represent the set of all tets
  moab::EntityHandle _meshset;

  /// Save some topological data: map from surface handle to vol handle and sense
  std::map<moab::EntityHandle, std::vector<VolData> > surfsToVols;

  // Some moab tags

  /// Tag for dimension for geometry
  moab::Tag geometry_dimension_tag;
  /// Tag for entitiy set ID
  moab::Tag id_tag;
  /// Tag for faceting tolerance
  moab::Tag faceting_tol_tag;
  /// Tag needed by DAGMC
  moab::Tag geometry_resabs_tag;
  /// Tag for type of entity set
  moab::Tag category_tag;
  /// Tag for name of entity set
  moab::Tag name_tag;

  /// Number of nodes per MOAB tet (which are first order, so TET4)
  const unsigned int NODES_PER_MOAB_TET = 4;

  /// Bounds of the temperature bins
  std::vector<Real> _temperature_bin_bounds;

  /// Bounds of the density bins
  std::vector<Real> _density_bin_bounds;

  /// Node ordering for a TET4 MOAB element, based on libMesh node numberings
  std::vector<std::vector<unsigned int>> _tet4_nodes;

  /**
   * Node ordering for eight TET4 MOAB elements, based on libMesh node numberings
   * for a TET10 element. We re-build the libMesh element into first-order MOAB elements.
   */
  std::vector<std::vector<unsigned int>> _tet10_nodes;

  /// Auxiliary variable number for temperature
  unsigned int _temperature_var_num;

  /// Auxiliary variable number for density
  unsigned int _density_var_num;

  const int INVALID_POINT_LOCATOR = -1;
};
