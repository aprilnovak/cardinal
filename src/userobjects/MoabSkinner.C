#include "MoabSkinner.h"

InputParameters
MoabSkinner::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRangeCheckedParam<Real>("faceting_tol", 1e-4, "faceting_tol > 0", "Faceting tolerance for DagMC");
  params.addRangeCheckedParam<Real>("geom_tol", 1e-6, "geom-tol > 0", "Geometry tolerance for DagMC");

  return params;
}

MoabSkinner::MoabSkinner(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _faceting_tol(getParam<Real>("faceting_tol")),
    _geom_tol(getParam<Real>("geom_tol"))
{
  // Create MOAB interface
  _moab = std::make_shared<moab::Core>();

  // Create a skinner
  _skinner = std::make_unique<moab::Skinner>(_moab.get());

  // Create a geom topo tool
  _gtt = std::make_unique<moab::GeomTopoTool>(_moab.get());
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

  /// Loop over the [Mesh] and build all the nodes as MOAB vertices
  auto node_id_to_handle = createNodes();
}

void
MoabSkinner::createTags()
{
  // TODO: will we ever hit any of these errors? If not, remove them.

  // First some built-in MOAB tag types
  auto rval = _moab->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
    _geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up geometry dimension tag");

  rval = _moab->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
    _id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up id tag");

  rval = _moab->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE,
    _category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up category tag");

  rval = _moab->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
    _name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up name tag");

  // Some tags needed for DagMC
  rval = _moab->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
    _faceting_tol_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up faceting tolerance tag");

  rval = _moab->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE,
    _geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set up geometry resabs tag");

  // Set the values for DagMC faceting / geometry tolerance tags on the mesh entity set
  rval = _moab->tag_set_data(_faceting_tol_tag, &_meshset, 1, &_faceting_tol);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set faceting tolerance tag on mesh entity"); // TODO: ever hit?

  rval = _moab->tag_set_data(_geometry_resabs_tag, &_meshset, 1, &_geom_tol);
  if (rval != moab::MB_SUCCESS)
    mooseError("Failed to set geometry tolerance tag on mesh entity"); // TODO: ever hit?
}

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
