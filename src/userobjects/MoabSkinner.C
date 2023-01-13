#include "MoabSkinner.h"
#include "UserErrorChecking.h"

#include "libmesh/enum_io_package.h"

// TODO: add ifdefs

registerMooseObject("CardinalApp", MoabSkinner);

InputParameters
MoabSkinner::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRangeCheckedParam<Real>("faceting_tol", 1e-4, "faceting_tol > 0", "Faceting tolerance for DagMC");
  params.addRangeCheckedParam<Real>("geom_tol", 1e-6, "geom_tol > 0", "Geometry tolerance for DagMC");

  params.addParam<bool>("output_full",  false, "Whether MOAB should write full mesh data to file");
  params.addParam<std::string>("output_base_full", "moab_full",
    "Base filename for full mesh file writes (will be appended by an integer for each unique write)");

  return params;
}

MoabSkinner::MoabSkinner(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _faceting_tol(getParam<Real>("faceting_tol")),
    _geom_tol(getParam<Real>("geom_tol")),
    _output_full(getParam<bool>("output_full")),
    _output_base_full(getParam<std::string>("output_base_full")),
    _n_write(0),
    _scaling(1.0)
{
  // Create MOAB interface
  _moab = std::make_shared<moab::Core>();

  // Create a skinner
  _skinner = std::make_unique<moab::Skinner>(_moab.get());

  // Create a geom topo tool
  _gtt = std::make_unique<moab::GeomTopoTool>(_moab.get());

  if (_output_full)
  {
    if (_output_base_full.empty())
      mooseError("'output_base_full' cannot be an empty string!");
  }
  else
    checkUnusedParam(parameters, "output_base_full", "'output_full' is false");
}

void
MoabSkinner::initialize()
{
  // Set spatial dimension in MOAB
  // TODO: What does "spatial_dimension" mean to MOAB? In libMesh, it conveys how many dimensions
  // of space the mesh is in (so a 2-D plane mesh could be in 3 spatial dimensions if it's tilted
  // with respect to the Cartesian axes)
  moab::ErrorCode rval = _moab->set_dimension(_fe_problem.mesh().getMesh().spatial_dimension());
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set MOAB mesh dimension"); // TODO: ever hit?

  rval = _moab->create_meshset(moab::MESHSET_SET, _meshset);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to create MOAB mesh set"); // TODO: ever hit?

  createTags();

  // Loop over the [Mesh] and build all the nodes as MOAB vertices and all the elements
  // as MOAB elements
  auto node_id_to_handle = createNodes();
  createElems(node_id_to_handle);
}

void
MoabSkinner::createTags()
{
  // TODO: will we ever hit any of these errors? If not, remove them.

  // First some built-in MOAB tag types
  auto rval = _moab->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
    _geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up geometry dimension tag"); // TODO

  rval = _moab->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
    _id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up id tag"); // TODO

  rval = _moab->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE,
    _category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up category tag"); // TODO

  rval = _moab->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
    _name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up name tag"); // TODO

  // Some tags needed for DagMC
  rval = _moab->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
    _faceting_tol_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up faceting tolerance tag"); // TODO

  rval = _moab->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE,
    _geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up geometry resabs tag"); // TODO

  // Set the values for DagMC faceting / geometry tolerance tags on the mesh entity set
  rval = _moab->tag_set_data(_faceting_tol_tag, &_meshset, 1, &_faceting_tol);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set faceting tolerance tag on mesh entity"); // TODO: ever hit?

  rval = _moab->tag_set_data(_geometry_resabs_tag, &_meshset, 1, &_geom_tol);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set geometry tolerance tag on mesh entity"); // TODO: ever hit?
}

// TODO: can delete this function, and just place it into createElems
std::map<dof_id_type, moab::EntityHandle>
MoabSkinner::createNodes()
{
  std::map<dof_id_type, moab::EntityHandle> node_id_to_handle;

  double coords[3];

  for (const auto & node : _fe_problem.mesh().getMesh().node_ptr_range())
  {
    coords[0] = _scaling * (*node)(0);
    coords[1] = _scaling * (*node)(1);
    coords[2] = _scaling * (*node)(2);

    // Add node to MOAB database and get handle
    moab::EntityHandle ent(0);
    _moab->create_vertex(coords, ent);

    node_id_to_handle[node->id()] = ent;
  }

  return node_id_to_handle;
}

void
MoabSkinner::clearElemMaps()
{
  _id_to_elem_handles.clear();
  _offset = 0;
}

std::vector<std::vector<unsigned int>>
MoabSkinner::getTetSets(ElemType type)
{
  std::vector<std::vector<unsigned int>> perms;

  // TODO: save these as member variables
  if (type == TET4)
    perms.push_back({0,1,2,3});
  else if (type == TET10)
  {
    perms.push_back({0,4,6,7});
    perms.push_back({1,5,4,8});
    perms.push_back({2,6,5,9});
    perms.push_back({7,8,9,3});
    perms.push_back({4,9,7,8});
    perms.push_back({4,5,9,8});
    perms.push_back({4,7,9,6});
    perms.push_back({4,9,5,6});
  }
  else
    mooseError("MoabSkinner requires the [Mesh] to be tetrahedral elements!"); // TODO: add test

  return perms;
}

void
MoabSkinner::addElem(dof_id_type id, moab::EntityHandle ent)
{
  if (_id_to_elem_handles.find(id) == _id_to_elem_handles.end())
    _id_to_elem_handles[id] = std::vector<moab::EntityHandle>();

  _id_to_elem_handles[id].push_back(ent);
}

void
MoabSkinner::createElems(std::map<dof_id_type, moab::EntityHandle> & node_id_to_handle)
{
  // TODO: can just paste function here, to delete clearElemMaps()
  clearElemMaps();

  moab::Range all_elems;

  for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range())
  {
    // Get all sub-tetrahedra node sets for this element type
    auto nodeSets = getTetSets(elem->type());

    // Get the connectivity
    std::vector< dof_id_type > conn_libmesh;
    elem->connectivity(0, libMesh::IOPackage::VTK, conn_libmesh);
    if (conn_libmesh.size() != elem->n_nodes())
      mooseError("Element connectivity is inconsistent"); // TODO: I don't think this could ever be hit

    // Loop over sub tets
    for (const auto& nodeSet: nodeSets)
    {
      // Set MOAB connectivity
      std::vector<moab::EntityHandle> conn(NODES_PER_MOAB_TET);
      for (unsigned int i = 0; i < NODES_PER_MOAB_TET; ++i)
      {
        // Get the elem node index of the ith node of the sub-tet
        unsigned int nodeIndex = nodeSet.at(i);

        if (nodeIndex >= conn_libmesh.size()) // TODO: I don't think this could be hit
          mooseError("Element index is out of range");

        // Get node's entity handle
        if (node_id_to_handle.find(conn_libmesh.at(nodeIndex)) == node_id_to_handle.end())
          mooseError("Could not find node entity handle"); // TODO: ever hit?

        conn[i] = node_id_to_handle[conn_libmesh.at(nodeIndex)];
      }

      // Create an element in MOAB database
      moab::EntityHandle ent(0);
      auto rval = _moab->create_element(moab::MBTET, conn.data(), NODES_PER_MOAB_TET, ent);
      if (rval != moab::MB_SUCCESS)
        mooseError("Failed to create MOAB element"); // TODO: ever hit?

      // Save mapping between libMesh ids and moab handles
      addElem(elem->id(), ent);

      // Save the handle for adding to entity sets
      all_elems.insert(ent);
    }
  }

  // Add the elems to the full meshset
  auto rval = _moab->add_entities(_meshset, all_elems);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to create MOAB mesh set");  // TODO: test?

  // Save the first elem
  _offset = all_elems.front();
}

void
MoabSkinner::writeFullMesh() const
{
  // only the root can write the file
  if (processor_id() != 0)
    return;

  std::string filename = _output_base_full + "_" + std::to_string(_n_write) +".h5m";
  _console << "Writing MOAB mesh to "<< filename << "... ";

  // TODO: would be nice to just write VTK or Exodus directly, instead of needing the mbconvert.
  // We currently test that we write the files, but this could be more rigorous if we can check
  // the actual mesh contents too.
  auto rval = _moab->write_mesh(filename.c_str());
  if (rval != moab::MB_SUCCESS) // TODO: ever hit?
    mooseError("Failed to write full MOAB mesh");

  _console << "done" << std::endl;
}

void
MoabSkinner::execute()
{
  write();
}

void
MoabSkinner::write()
{
  if (_output_full)
    writeFullMesh();

  _n_write++;
}

void
MoabSkinner::setTags(moab::EntityHandle ent, std::string name, std::string category, unsigned int id, int dim)
{
  if (name != "")
    setTagData(_name_tag, ent, name, NAME_TAG_SIZE);

  if (category != "")
    setTagData(_category_tag, ent, category, CATEGORY_TAG_SIZE);

  setTagData(_geometry_dimension_tag, ent, &dim);

  setTagData(_id_tag, ent, &id);
}

void
MoabSkinner::setTagData(moab::Tag tag, moab::EntityHandle ent, std::string data, unsigned int size)
{
  auto namebuf = new char[size];
  memset(namebuf,'\0', size);
  strncpy(namebuf, data.c_str(), size - 1);
  auto rval = _moab->tag_set_data(tag, &ent, 1, namebuf);

  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set tag data for tag '" + data + "'"); // TODO: hit?

  delete[] namebuf;
}

moab::ErrorCode
MoabSkinner::setTagData(moab::Tag tag, moab::EntityHandle ent, void * data)
{
  // TODO: can delete this function and just insert where called
  return _moab->tag_set_data(tag, &ent, 1, data);
}
