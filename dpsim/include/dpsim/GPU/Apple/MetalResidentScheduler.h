/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#pragma once

#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator4OrderVBR.h>
#include <dpsim/Scheduler.h>

namespace DPsim {

/// Experimental device-resident EMT scheduler for the IEEE118 paper-reduced
/// SSN configuration. It intentionally rejects components for which no Metal
/// timestep kernel exists.
class MetalResidentScheduler final : public Scheduler {
public:
  MetalResidentScheduler(
      CPS::EMT::LinearSSNGroup::Ptr group,
      std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>>
          generators,
      CPS::EMT::SimNode::List nodes, Real timeStep, UInt expectedSteps,
      CPS::String tracePath, Bool wholeRunFastPath = true,
      Bool useSparseInputGpuLu = false, Bool useGpuGeneratedIcb = false);
  ~MetalResidentScheduler() override;

  void createSchedule(const CPS::Task::List &tasks, const Edges &inEdges,
                      const Edges &outEdges) override;
  void step(Real time, Int timeStepCount) override;
  Bool runSimulation(Real startTime, Int firstStep, UInt stepCount) override;
  void stop() override;

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
