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

#include "CardinalEnums.h"

MooseEnum getNekOrderEnum()
{
  return MooseEnum("first second", "first");
}

MooseEnum getBinnedVelocityComponentEnum()
{
  return MooseEnum("normal user");
}

MooseEnum getNekFieldEnum()
{
  return MooseEnum("velocity_component velocity_x velocity_y velocity_z velocity temperature pressure unity");
}

MooseEnum getOperationEnum()
{
  return MooseEnum("max min", "max");
}

MooseEnum getTallyTypeEnum()
{
  return MooseEnum("cell mesh");
}

MooseEnum getEigenvalueEnum()
{
  return MooseEnum("collision absorption tracklength combined", "combined");
}

MooseEnum getChannelTypeEnum()
{
  return MooseEnum("interior edge corner");
}

MooseEnum getRelaxationEnum()
{
  return MooseEnum("constant robbins_monro dufek_gudowski none", "none");
}
