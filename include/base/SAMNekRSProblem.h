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

#include "NekRSProblemBase.h"
#include "NekTimeStepper.h"
#include "NekRSMesh.h"
#include "Transient.h"

#include <memory>

/**
 * \brief TODO description.
 *
 * ...
 * 
 */
class SAMNekRSProblem : public NekRSProblemBase
{
public:
  SAMNekRSProblem(const InputParameters & params);

  static InputParameters validParams();

  /**
   * \brief Write nekRS's solution at the last output step
   *
   * If Nek is not the master app, the number of time steps it takes is
   * controlled by the master app. Depending on the settings in the `.par` file,
   * it becomes possible that nekRS may not write an output file on the simulation's
   * actual last time step, because Nek may not know when that last time step is.
   * Therefore, here we can force nekRS to write its output.
   **/
  ~SAMNekRSProblem();

  /**
   * \brief Perform some sanity checks on the problem setup
   *
   * This function performs checks like making sure that a transient executioner is
   * used to wrap nekRS, that no time shift has been requested to the start of nekRS,
   * that the correct NekTimeStepper is used, etc.
   */
  virtual void initialSetup() override;

  /// Send boundary velocity to nekRS
  void sendBoundaryVelocityToNek();

  /// Send boundary temperature to nekRS
  void sendBoundaryTemperatureToNek();

//  /// Get boundary temperature from nekRS
//  void getBoundaryTemperatureFromNek();

  virtual void syncSolutions(ExternalProblem::Direction direction) override;

  virtual void addExternalVariables() override;

  /**
 * Whether data should be synchronized in to nekRS
 * \return whether inward data synchronization should occur
 */
  virtual bool synchronizeIn();

  /**
 * Whether data should be synchronized out of nekRS
 * \return whether outward data synchronization should occur
 */
  virtual bool synchronizeOut();



  virtual bool movingMesh() const override { return _moving_mesh; }

protected:
//  virtual void addTemperatureVariable() override { return; }

  std::unique_ptr<NumericVector<Number>> _serialized_solution;

  /// Whether the problem is a moving mesh problem i.e. with on-the-fly mesh deformation enabled
  const bool & _moving_mesh;

  const bool & _minimize_transfers_in;

  const bool & _minimize_transfers_out;


  /// Specify type of boundary/boundaries present for SAM-NekRS coupling 
  const bool & _SAMtoNekRS;
  const bool & _NekRStoSAM;
  const bool & _SAMtoNekRS_temperature;

  /// Velocity boundary condition coming from SAM to NekRS
  const PostprocessorValue * _SAMtoNekRS_velocity = nullptr;

  /// Temperature boundary condition coming from SAM to NekRS
  const PostprocessorValue * _SAMtoNekRS_temp = nullptr;


//  /**
//   * \brief Total surface-integrated flux coming from the coupled MOOSE app.
//   *
//   * The mesh used for the MOOSE app may be very different from the mesh used by nekRS.
//   * Elements may be much finer/coarser, and one element on the MOOSE app may not be a
//   * clear subset/superset of the elements on the nekRS mesh. Therefore, to ensure
//   * conservation of energy, we send the total flux integral to nekRS for internal
//   * normalization of the heat flux applied on the nekRS mesh.
//   */
//  const PostprocessorValue * _flux_integral = nullptr;

//  /// MOOSE flux interpolated onto the (boundary) data transfer mesh
//  double * _flux_face = nullptr;

  /// Postprocessor containing the signal of when a synchronization has occurred
  const PostprocessorValue * _transfer_in = nullptr;

  /// flag to indicate whether this is the first pass to serialize the solution
  static bool _first;
};
