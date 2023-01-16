#ifdef ENABLE_DAGMC

#include "MoabUserObject.h"
#include "DisplacedProblem.h"
#include "VariadicTable.h"
#include "AuxiliarySystem.h"
#include "BinUtility.h"
#include "GeometryUtility.h"
#include "UserErrorChecking.h"

#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("CardinalApp", MoabUserObject);

InputParameters
MoabUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addParam<bool>("verbose", false, "Whether to print diagnostic information");

  // temperature binning
  params.addRequiredParam<std::string>("temperature", "Temperature variable by which to bin elements");
  params.addRangeCheckedParam<Real>("temperature_min", 0.0, "temperature_min >= 0.0",
    "Lower bound of temperature bins");
  params.addRequiredParam<Real>("temperature_max", "Upper bound of temperature bins");
  params.addRequiredRangeCheckedParam<unsigned int>("n_temperature_bins", "n_temperature_bins > 0",
    "Number of temperature bins");

  // density binning
  params.addParam<std::string>("density", "Density variable by which to bin elements");
  params.addRangeCheckedParam<Real>("density_min", 0.0, "density_min >= 0.0", "Lower bound of density bins");
  params.addParam<Real>("density_max", "Upper bound of density bins");
  params.addRangeCheckedParam<unsigned int>("n_density_bins", "n_density_bins > 0", "Number of density bins");

  // Mesh metadata
  params.addParam<std::vector<std::string> >("material_names", "List of MOOSE material names");
  params.addParam<std::vector<std::string> >("material_openmc_names", std::vector<std::string>(), "List of OpenMC material names");

  params.addRangeCheckedParam<Real>("faceting_tol", 1e-4, "faceting_tol > 0", "Faceting tolerance for DagMC");
  params.addRangeCheckedParam<Real>("geom_tol", 1e-6, "geom_tol > 0", "Geometry tolerance for DagMC");

  params.addParam<bool>("build_graveyard", false, "Whether to build a graveyard around the geometry");
  params.addRangeCheckedParam<Real>("graveyard_scale_inner", 1.01, "graveyard_scale_inner > 1",
    "Multiplier on mesh bounding box to form inner graveyard surface");
  params.addParam<Real>("graveyard_scale_outer", 1.10,
    "Multiplier on mesh bounding box to form outer graveyard surface");

  params.addParam<bool>("output_skins", false, "Whether the skinned MOAB mesh (skins generated from the "
    "libMesh [Mesh]) should be written to a file. The files will be named moab_skins_<n>.h5m, where <n> "
    "is the time step index. You can then visualize these files by running 'mbconvert'.");
  params.addParam<bool>("output_full",  false, "Whether the MOAB mesh (copied from the libMesh [Mesh]) should "
    "be written to a file. The files will be named moab_full_<n>.h5m, where <n> is the time step index. "
    "You can then visualize these files by running 'mbconvert'.");

  return params;
}

MoabUserObject::MoabUserObject(const InputParameters & parameters) :
  GeneralUserObject(parameters),
  _serialized_solution(NumericVector<Number>::build(_communicator).release()),
  _verbose(getParam<bool>("verbose")),
  _build_graveyard(getParam<bool>("build_graveyard")),
  _temperature_name(getParam<std::string>("temperature")),
  _temperature_min(getParam<Real>("temperature_min")),
  _temperature_max(getParam<Real>("temperature_max")),
  _n_temperature_bins(getParam<unsigned int>("n_temperature_bins")),
  _temperature_bin_width((_temperature_max - _temperature_min) / _n_temperature_bins),
  _bin_by_density(isParamValid("density")),
  mat_names(getParam<std::vector<std::string> >("material_names")),
  openmc_mat_names(getParam<std::vector<std::string> >("material_openmc_names")),
  _faceting_tol(getParam<Real>("faceting_tol")),
  _geom_tol(getParam<Real>("geom_tol")),
  _graveyard_scale_inner(getParam<double>("graveyard_scale_inner")),
  _graveyard_scale_outer(getParam<double>("graveyard_scale_outer")),
  _output_skins(getParam<bool>("output_skins")),
  _output_full(getParam<bool>("output_full")),
  _scaling(1.0),
  _n_write(0)
{
  // Create MOAB interface
  _moab = std::make_shared<moab::Core>();

  // Create a skinner
  skinner = std::make_unique<moab::Skinner>(_moab.get());

  // Create a geom topo tool
  gtt = std::make_unique<moab::GeomTopoTool>(_moab.get());

  if (_bin_by_density)
  {
    checkRequiredParam(parameters, "density_min", "binning by density");
    checkRequiredParam(parameters, "density_max", "binning by density");
    checkRequiredParam(parameters, "n_density_bins", "binning by density");

    _density_min = getParam<Real>("density_min");
    _density_max = getParam<Real>("density_max");
    _n_density_bins = getParam<unsigned int>("n_density_bins");
    den_var_name = getParam<std::string>("density");
    _density_bin_width = (_density_max - _density_min) / _n_density_bins;

    if (_density_max < _density_min)
      paramError("density_max", "'density_max' must be greater than 'density_min'");
  }
  else
  {
    checkUnusedParam(parameters, "density_min", "not binning by density");
    checkUnusedParam(parameters, "density_max", "not binning by density");
    checkUnusedParam(parameters, "n_density_bins", "not binning by density");

    _n_density_bins = 1;
  }

  if (_build_graveyard)
  {
    if (_graveyard_scale_outer < _graveyard_scale_inner)
      paramError("graveyard_scale_outer", "'graveyard_scale_outer' must be greater than 'graveyard_scale_inner'!");
  }
  else
  {
    checkUnusedParam(parameters, "graveyard_scale_inner", "'build_graveyard' is false");
    checkUnusedParam(parameters, "graveyard_scale_outer", "'build_graveyard' is false");
  }

  // get variable numbers
  _temperature_var_num = getAuxiliaryVariableNumber(_temperature_name, "temperature");
  if (_bin_by_density)
  {
    if (_temperature_name == den_var_name)
      mooseError("The 'temperature' and 'density' variables cannot be the same!");

    _density_var_num = getAuxiliaryVariableNumber(den_var_name, "density");
  }

  // If no alternative names were provided for openmc materials
  // assume they are the same as in MOOSE
  if(openmc_mat_names.empty()){
    openmc_mat_names = mat_names;
  }
  if(openmc_mat_names.size() != mat_names.size() ){
    mooseError("If both are provided, the vectors material_names and material_openmc_names should have identical lengths.");
  }

  if (_temperature_max <= _temperature_min)
    paramError("temperature_max", "'temperature_max' must be greater than 'temperature_min'");

  for (unsigned int i = 0; i < _n_temperature_bins + 1; ++i)
    _temperature_bin_bounds.push_back(_temperature_min + i * _temperature_bin_width);

  for (unsigned int i = 0; i < _n_density_bins + 1; ++i)
    _density_bin_bounds.push_back(_density_min + i * _density_bin_width);

  _tet4_nodes.push_back({0,1,2,3});

  _tet10_nodes.push_back({0,4,6,7});
  _tet10_nodes.push_back({1,5,4,8});
  _tet10_nodes.push_back({2,6,5,9});
  _tet10_nodes.push_back({7,8,9,3});
  _tet10_nodes.push_back({4,9,7,8});
  _tet10_nodes.push_back({4,5,9,8});
  _tet10_nodes.push_back({4,7,9,6});
  _tet10_nodes.push_back({4,9,5,6});
}

unsigned int
MoabUserObject::getAuxiliaryVariableNumber(const std::string & name, const std::string & param_name) const
{
  if (!_fe_problem.getAuxiliarySystem().hasVariable(name))
    paramError(param_name, "Cannot find auxiliary variable '", name , "'!");

  return _fe_problem.getAuxiliarySystem().getFieldVariable<Real>(0, name).number();
}

MeshBase&
MoabUserObject::mesh()
{
  if(_fe_problem.haveDisplaced()){
    return _fe_problem.getDisplacedProblem()->mesh().getMesh();
  }
  return _fe_problem.mesh().getMesh();
}

EquationSystems&
MoabUserObject::systems()
{
  return _fe_problem.es();
}

System&
MoabUserObject::system(std::string var_now)
{
  return _fe_problem.getSystem(var_now);
}

void
MoabUserObject::initialize()
{
  // Fetch spatial dimension from libMesh
  int dim = mesh().spatial_dimension();

  // Set spatial dimension in MOAB
  moab::ErrorCode  rval = _moab->set_dimension(dim);
  if(rval!=moab::MB_SUCCESS)
    mooseError("Failed to set MOAB dimension");

  //Create a meshset
  rval = _moab->create_meshset(moab::MESHSET_SET,meshset);
  if(rval!=moab::MB_SUCCESS)
    mooseError("Failed to create mesh set");

  rval = createTags();
  if(rval!=moab::MB_SUCCESS)
    mooseError("Could not set up tags");

  std::map<dof_id_type,moab::EntityHandle> node_id_to_handle;
  rval = createNodes(node_id_to_handle);
  if(rval!=moab::MB_SUCCESS)
    mooseError("Could not create nodes");

  createElems(node_id_to_handle);

  findBlocks();
}

void
MoabUserObject::execute()
{
  std::cout << "execute" << std::endl;
  // Clear MOAB mesh data from last timestep
  reset();

  _serialized_solution->init(_fe_problem.getAuxiliarySystem().sys().n_dofs(), false, SERIAL);
  _fe_problem.getAuxiliarySystem().solution().localize(*_serialized_solution);

  // Re-initialise the mesh data
  initialize();

  // Sort libMesh elements into bins
  sortElemsByResults();

  // Find the surfaces of local temperature regions
  findSurfaces();
  std::cout << "done executing" << std::endl;
}

void
MoabUserObject::findBlocks()
{
  _blocks.clear();

  int i = 0;
  for (const auto & b : _fe_problem.mesh().meshSubdomains())
    _blocks[b] = i++;

  _n_block_bins = _blocks.size();
}

moab::ErrorCode
MoabUserObject::createNodes(std::map<dof_id_type,moab::EntityHandle>& node_id_to_handle)
{
  moab::ErrorCode rval(moab::MB_SUCCESS);

  // Clear prior results.
  node_id_to_handle.clear();

  // Init array for MOAB node coords
  double 	coords[3];

  // TODO think about how the mesh is distributed...
  // Iterate over nodes in libmesh
  auto itnode = mesh().nodes_begin();
  auto endnode = mesh().nodes_end();
  for( ; itnode!=endnode; ++itnode){
    // Fetch a const ref to node
    const Node& node = **itnode;

    // Fetch coords (and scale to correct units)
    coords[0]=_scaling*double(node(0));
    coords[1]=_scaling*double(node(1));
    coords[2]=_scaling*double(node(2));

    // Fetch ID
    dof_id_type id = node.id();

    // Add node to MOAB database and get handle
    moab::EntityHandle ent(0);
    rval = _moab->create_vertex(coords,ent);
    if(rval!=moab::MB_SUCCESS){
      node_id_to_handle.clear();
      return rval;
    }

    // Save mapping of ids.
    node_id_to_handle[id] = ent;

  }

  return rval;
}

void
MoabUserObject::createElems(std::map<dof_id_type,moab::EntityHandle>& node_id_to_handle)
{
  moab::ErrorCode rval(moab::MB_SUCCESS);

  // Clear prior results.
  clearElemMaps();

  moab::Range all_elems;

  // Iterate over elements in the mesh
  for (const auto & elem : mesh().active_element_ptr_range())
  {
    auto nodeSets = getTetSets(elem->type());

    // Get the connectivity
    std::vector< dof_id_type > conn_libmesh;
    elem->connectivity(0,libMesh::IOPackage::VTK,conn_libmesh);
    if(conn_libmesh.size()!=elem->n_nodes())
      mooseError("Element connectivity is inconsistent");

    // Loop over sub tets
    for(const auto& nodeSet: nodeSets){

      // Set MOAB connectivity
      std::vector<moab::EntityHandle> conn(NODES_PER_MOAB_TET);
      for(unsigned int iNode=0; iNode<NODES_PER_MOAB_TET;++iNode){

        // Get the elem node index of the ith node of the sub-tet
        unsigned int nodeIndex = nodeSet.at(iNode);

        if(nodeIndex >= conn_libmesh.size())
          mooseError("Element index is out of range");

        // Get node's entity handle
        if(node_id_to_handle.find(conn_libmesh.at(nodeIndex)) ==
           node_id_to_handle.end())
          mooseError("Could not find node entity handle");

        conn[iNode]=node_id_to_handle[conn_libmesh.at(nodeIndex)];
      }

      // Create an element in MOAB database
      moab::EntityHandle ent(0);
      rval = _moab->create_element(moab::MBTET,conn.data(),NODES_PER_MOAB_TET,ent);
      if(rval!=moab::MB_SUCCESS){
        std::string err="Could not create MOAB element: rval = "
          +std::to_string(rval);
        mooseError(err);
      }

      // Save mapping between libMesh ids and moab handles
      addElem(elem->id(),ent);

      // Save the handle for adding to entity sets
      all_elems.insert(ent);
    }
  }

  // Add the elems to the full meshset
  rval = _moab->add_entities(meshset,all_elems);
  if(rval!=moab::MB_SUCCESS){
    std::string err="Could not create meshset: rval = "
      +std::to_string(rval);
    mooseError(err);
  }

  // Save the first elem
  offset = all_elems.front();
}

const std::vector<std::vector<unsigned int>> &
MoabUserObject::getTetSets(ElemType type) const
{
  if (type == TET4)
    return _tet4_nodes;
  else if (type == TET10)
    return _tet10_nodes;
  else
    mooseError("The MoabUserObject can only be used with a tetrahedral [Mesh]!");
}


moab::ErrorCode
MoabUserObject::createTags()
{
  // Create some tags for later use
  moab::ErrorCode rval = moab::MB_SUCCESS;

  // First some built-in MOAB tag types
  rval = _moab->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS) return rval;

  rval = _moab->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS) return rval;

  rval = _moab->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS)  return rval;

  rval = _moab->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE, name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS)  return rval;

  // Some tags needed for DagMC
  rval = _moab->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS)  return rval;

  rval = _moab->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if(rval!=moab::MB_SUCCESS)  return rval;

  // Set the values for DagMC faceting / geometry tolerance tags on the mesh entity set
  rval = _moab->tag_set_data(faceting_tol_tag, &meshset, 1, &_faceting_tol);
  if(rval!=moab::MB_SUCCESS)  return rval;

  rval = _moab->tag_set_data(geometry_resabs_tag, &meshset, 1, &_geom_tol);
  return rval;
}

moab::ErrorCode
MoabUserObject::createGroup(unsigned int id, std::string name,moab::EntityHandle& group_set)
{
  // Create a new mesh set
  moab::ErrorCode rval = _moab->create_meshset(moab::MESHSET_SET,group_set);
  if(rval!=moab::MB_SUCCESS) return rval;

  // Set the tags for this material
  return setTags(group_set,name,"Group",id,4);
}


moab::ErrorCode
MoabUserObject::createVol(unsigned int id,moab::EntityHandle& volume_set,moab::EntityHandle group_set)
{
  moab::ErrorCode rval = _moab->create_meshset(moab::MESHSET_SET,volume_set);
  if(rval!=moab::MB_SUCCESS) return rval;

  rval =  setTags(volume_set,"","Volume",id,3);
  if(rval != moab::MB_SUCCESS) return rval;

  // Add the volume to group
  rval = _moab->add_entities(group_set, &volume_set,1);
  if(rval != moab::MB_SUCCESS) return rval;

  return rval;
}

void
MoabUserObject::createSurf(unsigned int id,moab::EntityHandle& surface_set, moab::Range& faces,  std::vector<VolData> & voldata)
{
  // Create meshset
  _moab->create_meshset(moab::MESHSET_SET,surface_set);

  // Set tags
  setTags(surface_set,"","Surface",id,2);

  // Add tris to the surface
  _moab->add_entities(surface_set,faces);

  // Create entry in map
  surfsToVols[surface_set] = std::vector<VolData>();

  // Add volume to list associated with this surface
  for(const auto & data : voldata)
    updateSurfData(surface_set,data);
}

moab::ErrorCode
MoabUserObject::updateSurfData(moab::EntityHandle surface_set,VolData data)
{
  // Add the surface to the volume set
  moab::ErrorCode rval = _moab->add_parent_child(data.vol,surface_set);
  if(rval != moab::MB_SUCCESS) return rval;

  // Set the surfaces sense
  rval = gtt->set_sense(surface_set,data.vol,int(data.sense));
  if(rval != moab::MB_SUCCESS) return rval;

  // Save
  surfsToVols[surface_set].push_back(data);

  return moab::MB_SUCCESS;
}


moab::ErrorCode
MoabUserObject::setTags(moab::EntityHandle ent, std::string name, std::string category, unsigned int id, int dim)
{

  moab::ErrorCode rval;

  // Set the name tag
  if(name!=""){
    rval = setTagData(name_tag,ent,name,NAME_TAG_SIZE);
    if(rval!=moab::MB_SUCCESS) return rval;
  }

  // Set the category tag
  if(category!=""){
    rval = setTagData(category_tag,ent,category,CATEGORY_TAG_SIZE);
    if(rval!=moab::MB_SUCCESS) return rval;
  }

  // Set the dimension tag
  rval = setTagData(geometry_dimension_tag,ent,&dim);
  if(rval!=moab::MB_SUCCESS) return rval;

  // Set the id tag
  rval = setTagData(id_tag,ent,&id);
  return rval;

}

moab::ErrorCode
MoabUserObject::setTagData(moab::Tag tag, moab::EntityHandle ent, std::string data, unsigned int SIZE)
{
  auto namebuf= new char[SIZE];
  memset(namebuf,'\0', SIZE); // fill C char array with null
  strncpy(namebuf,data.c_str(),SIZE-1);
  moab::ErrorCode rval = _moab->tag_set_data(tag,&ent,1,namebuf);
  // deallocate memory
  delete[] namebuf;
  return rval;
}

moab::ErrorCode
MoabUserObject::setTagData(moab::Tag tag, moab::EntityHandle ent, void* data)
{
  return _moab->tag_set_data(tag,&ent,1,data);
}

void
MoabUserObject::clearElemMaps()
{
  _id_to_elem_handles.clear();
  offset=0;
}

void
MoabUserObject::addElem(dof_id_type id,moab::EntityHandle ent)
{
  if(_id_to_elem_handles.find(id)==_id_to_elem_handles.end())
    _id_to_elem_handles[id]=std::vector<moab::EntityHandle>();

  (_id_to_elem_handles[id]).push_back(ent);
}

void MoabUserObject::getMaterialProperties(std::vector<std::string>& mat_names_out,
                                           std::vector<double>& initial_densities,
                                           std::vector<std::string>& tails,
                                           std::vector<MOABMaterialProperties>& properties)
{
  // We shouldn't be calling this if we didn't provide any materials
  if(openmc_mat_names.empty())
    mooseError("No material names were provided.");

  // Set the list of materials names we expect to find in openmc
  mat_names_out=openmc_mat_names;

  tails.clear();
  properties.clear();

  // Loop over density bins
  for(unsigned int iDen=0; iDen<_n_density_bins; iDen++){

    // Retrieve the relative density
    double rel_den = bin_utility::midpoint(iDen, _density_bin_bounds);

    // Loop over temperature bins
    for(unsigned int iVar=0; iVar<_n_temperature_bins; iVar++){
      // Retrieve the average bin temperature
      double temp = bin_utility::midpoint(iVar, _temperature_bin_bounds);

      // Get the name modifier
      int iNewMatBin = getMatBin(iVar,iDen);
      std::string tail ="_"+std::to_string(iNewMatBin);
      tails.push_back(tail);

      // Save material properties
      MOABMaterialProperties mat_props;
      mat_props.temp = temp;
      mat_props.rel_density = rel_den;
      properties.push_back(mat_props);
    }
  }
}


dof_id_type
MoabUserObject::elem_to_soln_index(const Elem& elem,unsigned int iSysNow,  unsigned int iVarNow)
{
    // Expect only one component, but check anyay
  unsigned int n_components = elem.n_comp(iSysNow,iVarNow);
  if(n_components != 1){
    throw std::runtime_error("Unexpected number of expected solution components");
  }

  // Get the degree of freedom number
  dof_id_type soln_index = elem.dof_number(iSysNow,iVarNow,0);

  return soln_index;
}

void
MoabUserObject::sortElemsByResults()
{
   // Clear any prior data
  resetContainers();

  // accumulate information for printing diagnostics
  std::vector<unsigned int> n_block_hits(_n_block_bins, 0);
  std::vector<unsigned int> n_temp_hits(_n_temperature_bins, 0);
  std::vector<unsigned int> n_density_hits(_n_density_bins, 0);

  for (unsigned int e = 0; e < _fe_problem.mesh().nElem(); ++e)
  {
    const Elem * const elem = _fe_problem.mesh().queryElemPtr(e);
    if (!elem)
      continue;

    Point p = elem->vertex_average();

    // bin by subdomain ID
    auto iMat = getSubdomainBin(elem);
    n_block_hits[iMat] += 1;

    // bin by density
    auto iDenBin = getDensityBin(elem);
    n_density_hits[iDenBin] += 1;

    // bin by temperature
    auto iBin = getTemperatureBin(elem);
    n_temp_hits[iBin] += 1;

    // Sort elem into a bin
    auto iSortBin = getBin(iBin, iDenBin, iMat);
    _elem_bins.at(iSortBin).insert(elem->id());
  }

  VariadicTable<unsigned int, std::string, unsigned int> vtt({"Bin", "Range (K)", "# Elems"});
  VariadicTable<unsigned int, std::string, unsigned int> vtd({"Bin", "Range (kg/m3)", "# Elems"});

  for (unsigned int i = 0; i < _n_temperature_bins; ++i)
    vtt.addRow(i, std::to_string(_temperature_bin_bounds[i]) + " to " +
      std::to_string(_temperature_bin_bounds[i + 1]), n_temp_hits[i]);

  for (unsigned int i = 0; i < _n_density_bins; ++i)
    vtd.addRow(i, std::to_string(_density_bin_bounds[i]) + " to " +
      std::to_string(_density_bin_bounds[i + 1]), n_density_hits[i]);

  if (_verbose)
  {
    _console << "Mapping of Elements to Temperature Bins:" << std::endl;
    vtt.print(_console);

    if (_bin_by_density)
    {
      _console << "\nMapping of Elements to Density Bins:" << std::endl;
      vtd.print(_console);
    }
  }
}

unsigned int
MoabUserObject::getTemperatureBin(const Elem * const elem) const
{
  auto dof = elem->dof_number(_fe_problem.getAuxiliarySystem().number(), _temperature_var_num, 0);
  auto value = (*_serialized_solution)(dof);

  // TODO: add option to truncate instead
  if (value < _temperature_min)
    mooseError("Variable '", _temperature_name, "' has value below minimum range of bins. "
      "Please decrease 'temperature_min'.\n\n"
      "  value: ", value, "\n  temperature_min: ", _temperature_min);

  if (value > _temperature_max)
    mooseError("Variable '", _temperature_name, "' has value above maximum range of bins. "
      "Please increase 'temperature_max'.\n\n"
      "  value: ", value, "\n  temperature_max: ", _temperature_max);

  return bin_utility::linearBin(value, _temperature_bin_bounds);
}

unsigned int
MoabUserObject::getDensityBin(const Elem * const elem) const
{
  if (!_bin_by_density)
    return 0;

  auto dof = elem->dof_number(_fe_problem.getAuxiliarySystem().number(), _density_var_num, 0);
  auto value = (*_serialized_solution)(dof);

  // TODO: add option to truncate instead
  if (value < _density_min)
    mooseError("Variable '", den_var_name, "' has value below minimum range of bins. "
      "Please decrease 'density_min'.\n\n"
      "  value: ", value, "\n  density_min: ", _density_min);

  if (value > _density_max)
    mooseError("Variable '", den_var_name, "' has value above maximum range of bins. "
      "Please increase 'density_max'.\n\n"
      "  value: ", value, "\n  density_max: ", _density_max);

  return bin_utility::linearBin(value, _density_bin_bounds);
}

bool
MoabUserObject::findSurfaces()
{

  moab::ErrorCode rval = moab::MB_SUCCESS;
  try{
    // Find all neighbours in mesh
    mesh().find_neighbors();

    // Counter for volumes
    unsigned int vol_id=0;

    // Counter for surfaces
    unsigned int surf_id=0;

    // Loop over material bins
    for(unsigned int iMat=0; iMat<_n_block_bins; iMat++){

      // Get the base material name:
      std::string mat_name = "mat:"+openmc_mat_names.at(iMat);

      // Loop over density bins
      for(unsigned int iDen=0; iDen<_n_density_bins; iDen++){

        // Loop over temperature bins
        for(unsigned int iVar=0; iVar<_n_temperature_bins; iVar++){

          // Update material name
          std::string updated_mat_name=mat_name;
          int iNewMatBin = getMatBin(iVar,iDen);
          updated_mat_name+="_"+std::to_string(iNewMatBin);

          // Create a material group
          // Todo set temp in metadata?
          int iSortBin = getBin(iVar,iDen,iMat);
          moab::EntityHandle group_set;
          unsigned int group_id = iSortBin+1;
          rval = createGroup(group_id,updated_mat_name,group_set);
          if(rval != moab::MB_SUCCESS) return false;

          // Sort elems in this mat-density-temp bin into local regions
          std::vector<moab::Range> regions;
          groupLocalElems(_elem_bins.at(iSortBin),regions);

          // Loop over all regions and find surfaces
          for(const auto & region : regions){
            moab::EntityHandle volume_set;
            if(!findSurface(region,group_set,vol_id,surf_id,volume_set)){
              return false;
            }

          } // End loop over local regions

        } // End loop over temperature bins
      } // End loop over density bins
    } // End loop over materials

    // Finally, build a graveyard
    if (_build_graveyard)
      buildGraveyard(vol_id,surf_id);
  }
  catch(std::exception &e){
    std::cerr<<e.what()<<std::endl;
    return false;
  }

  // Write MOAB volume and/or skin meshes to file
  write();

  return true;
}

void
MoabUserObject::write()
{
  // Only write to file on root process
  if (processor_id() != 0)
    return;

  if (_output_skins)
  {
    // Generate list of surfaces to write
    std::vector<moab::EntityHandle> surfs;
    for(const auto & itsurf : surfsToVols)
      surfs.push_back(itsurf.first);

    std::string filename = "moab_skins_" + std::to_string(_n_write) +".h5m";

    if (_verbose)
      _console << "Writing MOAB skins to "<< filename << "...";

    _moab->write_mesh(filename.c_str(), surfs.data(), surfs.size());

    if (_verbose)
      _console << "done" << std::endl;
  }

  if (_output_full)
  {
    std::string filename = "moab_full_" + std::to_string(_n_write) +".h5m";

    if (_verbose)
      _console << "Writing MOAB mesh to "<< filename << std::endl;

    _moab->write_mesh(filename.c_str());

    if (_verbose)
      _console << "done" << std::endl;
  }

  _n_write++;
}

void
MoabUserObject::groupLocalElems(std::set<dof_id_type> elems, std::vector<moab::Range>& localElems)
{
  while(!elems.empty()){

    // Create a new local range of moab handles
    moab::Range local;

    // Retrieve and remove the fisrt elem
    auto it = elems.begin();
    dof_id_type next = *it;
    elems.erase(it);

    std::set<dof_id_type> neighbors;
    neighbors.insert(next);

    while(!neighbors.empty()){

      std::set<dof_id_type> new_neighbors;

      // Loop over all the new neighbors
      for(auto& next : neighbors){

        // Get the MOAB handles, and add to local set
        // (May be more than one if this libMesh elem has sub-tetrahedra)
        if(_id_to_elem_handles.find(next)==_id_to_elem_handles.end()){
          mooseError("No entity handles found for libmesh id.");
        }
        std::vector<moab::EntityHandle> ents = _id_to_elem_handles[next];
        for(const auto ent : ents){
          local.insert(ent);
        }

        // Get the libMesh element
        Elem& elem = mesh().elem_ref(next);

        // How many nearest neighbors (general element)?
        unsigned int NN = elem.n_neighbors();

        // Loop over neighbors
        for(unsigned int i=0; i<NN; i++){

          const Elem * nnptr = elem.neighbor_ptr(i);
          // If on boundary, some may be null ptrs
          if(nnptr == nullptr) continue;

          dof_id_type idnn = nnptr->id();

          // Select only those that are in the current bin
          if(elems.find(idnn)!= elems.end()){
            new_neighbors.insert(idnn);
            // Remove from those still available
            elems.erase(idnn);
          }

        }// End loop over new neighbors

      }// End loop over previous neighbors

      // Found all the new neighbors, done with current set.
      neighbors = new_neighbors;

    }
    // Done, no more local neighbors in the current bin.

    // Save this moab range of local neighbors
    localElems.push_back(local);
  }
  // Done, assigned all elems in bin to a local range.
 }

void
MoabUserObject::resetContainers()
{
  unsigned int nSortBins = _n_block_bins*_n_density_bins*_n_temperature_bins;
  _elem_bins.clear();
  _elem_bins.resize(nSortBins);

  // Update the serial solutions
  for(const auto& sol :  serial_solutions){
    System & sys = 	systems().get_system(sol.first);

    // Check if solution vector size has changed, e.g. due to mesh refinement
    if(sys.n_dofs() != sol.second->size()){
      // clear
      sol.second->init(0,false,SERIAL);
      // resize
      sol.second->init(sys.n_dofs(),false,SERIAL);
    }

    sys.solution->localize(*sol.second);
  }
}

void
MoabUserObject::reset()
{
  // Clear data
  _moab.reset(new moab::Core());

  // Create a skinner and geometry topo tool
  skinner.reset(new moab::Skinner(_moab.get()));
  gtt.reset(new moab::GeomTopoTool(_moab.get()));

  // Clear entity set maps
  surfsToVols.clear();
}

unsigned int
MoabUserObject::getBin(const unsigned int & iVarBin, const unsigned int & iDenBin, const unsigned int & iMat) const
{
  return _n_temperature_bins * (_n_density_bins * iMat + iDenBin) + iVarBin;
}

int
MoabUserObject::getMatBin(int iVarBin, int iDenBin, int n_temperature_binsIn, int _n_density_binsIn)
{

  if(iDenBin<0 || iDenBin >= _n_density_binsIn ){
    std::string err = "Relative density of material fell outside of binning range";
    mooseError(err);
  }
  if(iVarBin<0 || iVarBin >= n_temperature_binsIn ){
    std::string err = "Relative temperature of material fell outside of binning range";
    mooseError(err);
  }

  int _n_block_bins = _n_density_binsIn*_n_temperature_bins;
  int iMatBin= n_temperature_binsIn*iDenBin + iVarBin;
  if(iMatBin<0 || iMatBin >= _n_block_bins){
    mooseError("Cannot find material bin index.");
  }
  return iMatBin;
}

bool
MoabUserObject::findSurface(const moab::Range& region,moab::EntityHandle group, unsigned int & vol_id, unsigned int & surf_id,moab::EntityHandle& volume_set)
{

  moab::ErrorCode rval;

  // Create a volume set
  vol_id++;
  rval = createVol(vol_id,volume_set,group);
  if(rval != moab::MB_SUCCESS) return false;

  // Find surfaces from these regions
  moab::Range tris; // The tris of the surfaces
  moab::Range rtris;  // The tris which are reversed with respect to their surfaces
  rval = skinner->find_skin(0,region,false,tris,&rtris);
  if(rval != moab::MB_SUCCESS) return false;
  if(tris.size()==0 && rtris.size()==0) return false;


  // Create surface sets for the forwards tris
  VolData vdata = {volume_set,Sense::FORWARDS};
  rval = createSurfaces(tris,vdata,surf_id);
  if(rval != moab::MB_SUCCESS) return false;

  // Create surface sets for the reversed tris
  vdata.sense =Sense::BACKWARDS;
  rval = createSurfaces(rtris,vdata,surf_id);
  if(rval != moab::MB_SUCCESS) return false;

  return true;
}


moab::ErrorCode
MoabUserObject::createSurfaces(moab::Range& faces, VolData& voldata, unsigned int& surf_id){

  moab::ErrorCode rval = moab::MB_SUCCESS;

  if(faces.empty()) return rval;

  // Loop over the surfaces we have already created
  for ( const auto & surfpair : surfsToVols ) {

    // Local copies of surf/vols
    moab::EntityHandle surf = surfpair.first;
    std::vector<VolData> vols = surfpair.second;

    // First get the entities in this surface
    moab::Range tris;
    rval = _moab->get_entities_by_handle(surf,tris);
    if(rval!=moab::MB_SUCCESS) return rval;

    // Find any tris that live in both surfs
    moab::Range overlap = moab::intersect(tris,faces);
    if(!overlap.empty()) {

      // Check if the tris are a subset or the entire surf
      if(tris.size()==overlap.size()){
        // Whole surface -> Just need to update the volume relationships
        rval = updateSurfData(surf,voldata);
      }
      else{
        // If overlap is subset, subtract shared tris from this surface and create a new shared surface
        rval = _moab->remove_entities(surf,overlap);
        if(rval!=moab::MB_SUCCESS) return rval;

        // Append our new volume to the list that share this surf
        vols.push_back(voldata);

        // Create a new shared surface
        moab::EntityHandle shared_surf;
        surf_id++;
        createSurf(surf_id,shared_surf,overlap,vols);
      }

      // Subtract from the input list
      for( auto& shared : overlap ){
        faces.erase(shared);
      }
      if(faces.empty()) break;
    }
  }

  if(!faces.empty()){
    moab::EntityHandle surface_set;
    std::vector<VolData> voldatavec(1,voldata);
    surf_id++;
    createSurf(surf_id,surface_set,faces,voldatavec);
  }

  return rval;
}

void
MoabUserObject::buildGraveyard( unsigned int & vol_id, unsigned int & surf_id)
{
  moab::ErrorCode rval(moab::MB_SUCCESS);

  // Create the graveyard set
  moab::EntityHandle graveyard;
  unsigned int id = _n_block_bins * _n_temperature_bins * _n_density_bins + 1;
  std::string mat_name = "mat:Graveyard";
  createGroup(id,mat_name,graveyard);

  // Create a volume set
  moab::EntityHandle volume_set;
  vol_id++;
  createVol(vol_id,volume_set,graveyard);

  // Set up for the volume data to pass to surfs
  VolData vdata = {volume_set,Sense::FORWARDS};

  // Find a bounding box
  BoundingBox bbox =  MeshTools::create_bounding_box(mesh());

  // Create inner surface with normals pointing into box
  createSurfaceFromBox(bbox,vdata,surf_id,false,_graveyard_scale_inner);

  // Create outer surface with face normals pointing out of the box
  createSurfaceFromBox(bbox,vdata,surf_id,true,_graveyard_scale_outer);
}

void
MoabUserObject::createSurfaceFromBox(const BoundingBox& box, VolData& voldata, unsigned int& surf_id, bool normalout, const Real & factor)
{
  std::vector<moab::EntityHandle> vert_handles = createNodesFromBox(box, factor);

  // Create the tris in 4 groups of 3 (4 open tetrahedra)
  moab::Range tris;
  createCornerTris(vert_handles,0,1,2,4,normalout,tris);
  createCornerTris(vert_handles,3,2,1,7,normalout,tris);
  createCornerTris(vert_handles,6,4,2,7,normalout,tris);
  createCornerTris(vert_handles,5,1,4,7,normalout,tris);

  moab::EntityHandle surface_set;
  std::vector<VolData> voldatavec(1,voldata);
  surf_id++;
  createSurf(surf_id,surface_set,tris,voldatavec);
}

std::vector<moab::EntityHandle>
MoabUserObject::createNodesFromBox(const BoundingBox & box, const Real & factor) const
{
  std::vector<moab::EntityHandle> vert_handles;

  // Fetch the vertices of the box
  auto verts = geom_utility::boxCorners(box, factor);

  // Array to represent a coord in MOAB
  double coord[3];

  // Create the vertices in MOAB and get the handles
  for(const auto & vert : verts)
  {
    coord[0] = vert(0) * _scaling;
    coord[1] = vert(1) * _scaling;
    coord[2] = vert(2) * _scaling;

    moab::EntityHandle ent;
    _moab->create_vertex(coord, ent);
    vert_handles.push_back(ent);
  }

  return vert_handles;
}

void
MoabUserObject::createCornerTris(const std::vector<moab::EntityHandle> & verts,
                                 unsigned int corner,
                                 unsigned int v1, unsigned int v2 ,unsigned int v3,
                                 bool normalout, moab::Range &surface_tris)
{
  // Create 3 tris stemming from one corner (i.e. an open tetrahedron)
  // Assume first is the central corner, and the others are labelled clockwise looking down on the corner
  unsigned int indices[3] = {v1,v2,v3};

  //Create each tri by a cyclic permutation of indices
  for(unsigned int i=0; i<3; i++)
  {
    // v1,v2 = 0,1; 1,2; 2;0
    int i1 = indices[i%3];
    int i2 = indices[(i+1)%3];
    if(normalout) // anti-clockwise: normal points outwards
      surface_tris.insert(createTri(verts,corner,i2,i1));
    else // clockwise: normal points inwards
      surface_tris.insert(createTri(verts,corner,i1,i2));
  }
}

moab::EntityHandle
MoabUserObject::createTri(const std::vector<moab::EntityHandle> & vertices,unsigned int v1, unsigned int v2 ,unsigned int v3)
{
  moab::EntityHandle triangle;
  moab::EntityHandle connectivity[3] = { vertices[v1],vertices[v2],vertices[v3] };
  _moab->create_element(moab::MBTRI,connectivity,3,triangle);
  return triangle;
}

bool
MoabUserObject::isLocalElem(const Elem * elem) const
{
  if (!elem)
  {
    // we should only not be able to find an element if the mesh is distributed
    libmesh_assert(!_mesh.is_serial());
    return false;
  }

  if (elem->processor_id() == _communicator.rank())
    return true;

  return false;
}

#endif
