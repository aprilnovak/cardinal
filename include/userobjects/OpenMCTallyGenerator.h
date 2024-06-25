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

#include "GeneralUserObject.h"

/**
 * User object to modify the nuclides in an OpenMC tally.
 */
class OpenMCTallyGenerator : public GeneralUserObject
{
public:
  static InputParameters validParams();

  OpenMCTallyGenerator(const InputParameters & parameters);

  /// We don't want this user object to execute in MOOSE's control
  virtual void execute() override;

  virtual void initialize() override {}
  virtual void finalize() override {}

protected:
  bool _nuclides_uo {false};

  std::vector<int32_t> _created_tally_ids;

};
