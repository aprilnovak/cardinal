#pragma once

#include "GeneralUserObject.h"

class MoabSkinner : public GeneralUserObject
{
 public:

  MoabSkinner(const InputParameters & parameters);

  static InputParameters validParams();
};
