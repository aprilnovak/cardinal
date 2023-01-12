#include "MoabSkinner.h"

InputParameters
MoabSkinner::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  return params;
}

MoabSkinner::MoabSkinner(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}
