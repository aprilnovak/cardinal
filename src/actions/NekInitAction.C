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

#include "NekInitAction.h"
#include "NekInterface.h"
#include "nekrs.hpp"

#include "MooseApp.h"
#include "CommandLine.h"

registerMooseAction("CardinalApp", NekInitAction, "nek_init");

int NekInitAction::_n_cases = 0;

InputParameters
NekInitAction::validParams()
{
  InputParameters params = MooseObjectAction::validParams();
  params.addParam<std::string>("type", "Problem type");
  params.addParam<std::string>("casename", "Case name for the NekRS input files; "
    "this is <case> in <case>.par, <case>.udf, <case>.oudf, and <case>.re2. "
    "Can also be provided on the command line with --nekrs-setup, which will override this setting");
  return params;
}

NekInitAction::NekInitAction(const InputParameters & parameters)
  : MooseObjectAction(parameters),
    _casename_in_input_file(isParamValid("casename"))
{
}

void
NekInitAction::act()
{
  if (_type == "NekRSProblem" || _type == "NekProblem" || _type == "NekRSStandaloneProblem")
  {
    std::shared_ptr<CommandLine> cl = _app.commandLine();
    bool casename_on_command_line = cl->search("nekrs_setup");

    if (!casename_on_command_line && !_casename_in_input_file)
      mooseError("All inputs using 'NekRSProblem' or 'NekProblem' must pass '--nekrs-setup <case>' on "
        "the command line\nor set casename = '<case>' in the [Problem] block in the Nek-wrapped input file!");

    std::string setup_file;
    if (casename_on_command_line)
      cl->search("nekrs_setup", setup_file);
    else
      setup_file = getParam<std::string>("casename");

    std::string cache_dir;

    // we need to set the default values here because it seems that the default values
    // that can be set via addCommandLineParam in CardinalApp aren't propagated through the 'search' function
    int size_target = 0;
    int ci_mode = 0;

    cl->search("nekrs_buildonly", size_target);
    cl->search("nekrs_cimode", ci_mode);

    int build_only = size_target > 0 ? 1: 0;

    nekrs::buildOnly(build_only);

    MPI_Comm comm = *static_cast<const MPI_Comm *>(&_communicator.get());

    // If this MPI communicator has already created one case, then we cannot also create a
    // second NekRS case. For instance, if you have 4 Nek sub-apps, but only 3 processes, then
    // NekRS doesn't like trying to set up a case with a communicator, then immediately try
    // to set up another case with the same communicator. If you un-comment this error message,
    // you'll get an error when reading the mesh file that is basically missing a "/".
    // This error is probably due to NekRS relying a lot on static variables.
    if (_n_cases > 0)
    {
      int size;
      MPI_Comm_size(comm, &size);
      mooseError("NekRS does not currently support setting up multiple cases with the same "
                 "MPI communicator.\nThat is, you need at least one MPI process in a master "
                 "application per Nek sub-application.\n\n"
                 "The MPI communicator has " + std::to_string(size) + " ranks and is trying to "
                 "construct " + std::to_string(_n_cases + 1) + "+ cases.");
    }

    nekrs::setup(comm, build_only, size_target, ci_mode, cache_dir, setup_file,
      "" /* backend */, "" /* device ID */);

    _n_cases++;
  }

  // setup actions only needed if coupling with MOOSE
  if (_type == "NekRSProblem" || _type == "NekProblem")
  {
    // First check we should do is that a temperature variable exists, or else many
    // of our indexes into `nrs->cds` would give seg faults
    if (!nekrs::hasTemperatureVariable())
      mooseError("To properly transfer temperature and heat flux between nekRS and MOOSE, "
        "your nekRS model must include a solution for temperature.\n\nDid you forget the "
        "TEMPERATURE block in the .par file?");

    // The legacy approach is based on Nek5000 backend, which initialized
    // the flux scratch space in the `.udf` file. Here, rather than force everyone
    // to update to using `NekRSProblem` and `NekRSMesh`, we will use the detection
    // of a used scratch space to figure out whether the legacy approach is being
    // used, and then warn people that it's going to be deprecated.
    bool scratch_available = nekrs::scratchAvailable();

    if (!scratch_available)
      mooseDeprecated("The 'NekProblem' and 'NekMesh' have been replaced by the "
        "'NekRSProblem' and 'NekRSMesh' classes! Please update your inputs files to "
        "use these new classes, as they will become deprecated in the future.\n\n"
        "How to update:\n"
        " - swap 'NekProblem' for 'NekRSProblem' and 'NekMesh' for 'NekRSMesh'\n"
        " - remove the manual allocation of scratch space (nrs->usrwrk and nrs->o_usrwrk) in udf\n"
        " - ensure that the sideset IDs are included in the .re2 file\n"
        " - specify what boundary ID you want to couple nekRS through, and set it with 'boundary' in 'NekRSMesh'\n"
        " - rename the 'total_flux' postprocessor to 'flux_integral'\n");

    // TODO: once all legacy inputs are updated, this error can be included without breaking anyone's inputs
    //if (!scratch_available)
    //  mooseError("The nrs_t.usrwrk and nrs_t.o_usrwrk arrays are reserved in Cardinal "
    //    "for transferring the heat flux fom MOOSE to nekRS!\n\nThis memory allocation is done "
    //    "automatically, but you have already initialized these arrays, perhaps for a different "
    //    "purpose than transferring flux. Please remove "
    //    "any UDF usage of these arrays to continue.");

    // For the new approach, initialize the space for the user scratch space on the host and device.
    // nekRS does not have a dedicated array to store a "flux" boundary condition, so we
    // will make use of these arrays to write incoming flux values from MOOSE.
    if (scratch_available)
     nekrs::initializeScratch();

    // Initialize default dimensional scales assuming a dimensional run is performed.
    // We need to do this construction here so that any tests that *dont* use NekRSProblem
    // (such as many of the mesh tests) have these scales initialized. If nekRS does indeed
    // run in nondimensional form, then NekRSProblem is needed to specify those scales,
    // and there we will again call initializeDimensionalScales.
    nekrs::solution::initializeDimensionalScales(1.0 /* U_ref */, 0.0 /* T_ref */,
      1.0 /* dT_ref */, 1.0 /* L_ref */, 1.0 /* rho_ref */, 1.0 /* Cp_ref */);
  }
}
