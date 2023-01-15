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

#pragma once

#include "AuxKernel.h"
#include "MooseEnum.h"
#include "MoabUserObject.h"

/**
 * Auxkernel to display the mapping of [Mesh] elements to the spatial
 * bins created by a MOAB mesh skinner.
 */
class MoabSkinnedBins : public AuxKernel
{
public:
  MoabSkinnedBins(const InputParameters & parameters);

  static InputParameters validParams();

  int variableBin() const;
  int materialBin() const;
  int densityBin() const;

protected:
  virtual Real computeValue();

  const MoabUserObject * _skinner;

  const MooseEnum _skin_by;
};
