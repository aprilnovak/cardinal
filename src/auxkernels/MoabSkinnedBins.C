/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#ifdef ENABLE_DAGMC

#include "MoabSkinnedBins.h"

registerMooseObject("CardinalApp", MoabSkinnedBins);

InputParameters
MoabSkinnedBins::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("skinner", "MOAB mesh skinner");

  MooseEnum skin_type("variable material density all", "all");
  params.addParam<MooseEnum>("skin_by", skin_type, "Which skin distribution to display");
  params.addClassDescription("Display the mapping of mesh elements to the skinned bins created by the MOAB skinner");
  return params;
}

MoabSkinnedBins::MoabSkinnedBins(const InputParameters & parameters) :
  AuxKernel(parameters),
  _skin_by(getParam<MooseEnum>("skin_by"))
{
  const UserObject & base = getUserObjectBase("skinner");
  _skinner = dynamic_cast<const MoabUserObject *>(&base);
  if (!_skinner)
    paramError("skinner", "This userobject must be of type MoabUserObject!");

  // TODO: test for correct variable type, since each element is only going to fall into one bin,
  // the applied variable should be constant monomial

  // TODO: check that the skinner has enabled the requested option
}

Real
MoabSkinnedBins::computeValue()
{
  if (_skin_by == "variable")
    return temperatureBin();
  else if (_skin_by == "material")
    return materialBin();
  else if (_skin_by == "density")
    return densityBin();
  else if (_skin_by == "all")
    return _skinner->getBin(temperatureBin(), densityBin(), materialBin());
  else
    mooseError("Unhandled skin_type enum in MoabSkinnedBins!");
}

int
MoabSkinnedBins::temperatureBin() const
{
  Point pt = _current_elem->vertex_average();
  return _skinner->getTemperatureBin(pt);
}

int
MoabSkinnedBins::materialBin() const
{
  // TODO: make this more efficient
  auto mat_blocks = _skinner->getMaterialBlocks();

  auto block = _current_elem->subdomain_id();
  for (unsigned int i = 0; i < mat_blocks.size(); ++i)
  {
    if (mat_blocks[i].find(block) == mat_blocks[i].end())
      return i;
  }

  mooseError("Failed to map block to a material!");
}

int
MoabSkinnedBins::densityBin() const
{
  auto iMat = materialBin();

  Point pt = _current_elem->vertex_average();
  return _skinner->getDensityBin(pt, iMat);
}

#endif
