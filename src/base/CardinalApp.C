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

#include "CardinalApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "OpenMCSyntax.h"
#include "NekSyntax.h"

#ifdef ENABLE_SAM_COUPLING
#include "SamApp.h"
#endif

#ifdef ENABLE_SOCKEYE_COUPLING
#include "SockeyeApp.h"
#endif

#ifdef ENABLE_THM_COUPLING
#include "THMApp.h"
#endif

#ifdef ENABLE_SODIUM
#include "SodiumApp.h"
#endif

#ifdef ENABLE_POTASSIUM
#include "PotassiumApp.h"
#endif

#ifdef ENABLE_IAPWS95
#include "IAPWS95App.h"
#endif

registerKnownLabel("CardinalApp");

InputParameters
CardinalApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // only used for Nek wrappings - if used with another application (OpenMC wrapping
  // or just plain MOOSE-type apps), these are unused
  params.addCommandLineParam<std::string>(
    "nekrs_setup",  "--nekrs-setup [nekrs_setup]",
    "Specify NekRS setup file (basename for .par, .re2, .udf, and .oudf files)");
  params.addCommandLineParam<int>(
    "nekrs_buildonly",  "--nekrs-buildonly [#procs]",
    "#procs to build NekRS if pre-compiling");
  params.addCommandLineParam<int>(
    "nekrs_cimode",  "--nekrs-cimode [id]",
    "CI test ID for NekRS");

  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

CardinalApp::CardinalApp(InputParameters parameters) : MooseApp(parameters)
{
  ModulesApp::registerAll(_factory, _action_factory, _syntax);
  CardinalApp::registerAll(_factory, _action_factory, _syntax);

#ifdef ENABLE_SAM_COUPLING
  SamApp::registerApps();
  SamApp::registerAll(_factory, _action_factory, _syntax);
#endif

#ifdef ENABLE_SOCKEYE_COUPLING
  SockeyeApp::registerApps();
  SockeyeApp::registerAll(_factory, _action_factory, _syntax);
#endif

#ifdef ENABLE_THM_COUPLING
  THMApp::registerApps();
  THMApp::registerAll(_factory, _action_factory, _syntax);
#endif

#ifdef ENABLE_SODIUM
  SodiumApp::registerApps();
  SodiumApp::registerAll(_factory, _action_factory, _syntax);
#endif

#ifdef ENABLE_POTASSIUM
  PotassiumApp::registerApps();
  PotassiumApp::registerAll(_factory, _action_factory, _syntax);
#endif

#ifdef ENABLE_IAPWS95
  IAPWS95App::registerApps();
  IAPWS95App::registerAll(_factory, _action_factory, _syntax);
#endif
}

void
CardinalApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"CardinalApp"});
  Registry::registerActionsTo(af, {"CardinalApp"});

  /* register custom execute flags, action syntax, etc. here */
  OpenMC::associateSyntax(s, af);
  Nek::associateSyntax(s, af);
}

void
CardinalApp::registerApps()
{
  registerApp(CardinalApp);

#ifdef ENABLE_SAM_COUPLING
  SamApp::registerApps();
#endif

#ifdef ENABLE_SOCKEYE_COUPLING
  SockeyeApp::registerApps();
#endif

#ifdef ENABLE_THM_COUPLING
  THMApp::registerApps();
#endif

#ifdef ENABLE_SODIUM
  SodiumApp::registerApps();
#endif

#ifdef ENABLE_POTASSIUM
  PotassiumApp::registerApps();
#endif

#ifdef ENABLE_IAPWS95
  IAPWS95App::registerApps();
#endif
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
CardinalApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  CardinalApp::registerAll(f, af, s);
}

extern "C" void
CardinalApp__registerApps()
{
  CardinalApp::registerApps();
}
