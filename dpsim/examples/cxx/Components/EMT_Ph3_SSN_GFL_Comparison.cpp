// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <cmath>
#include <stdexcept>
#include <utility>

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

/// \brief EMT comparison example for the refactored SSN_GFL inverter.
///
/// This example intentionally uses the same topology, parameters, power-flow
/// initialization, load-switch disturbance, and logger signal names as the
/// AvVoltSourceInverterStateSpace reference example. Consequently, the CSV
/// results of both examples can be compared directly.
class Example_SSN_GFL {
public:
  Example_SSN_GFL()
      : mTimeStepEMT(100e-6), mFinalTimeEMT(10.0), mLoadSwitchTime(3.0),
        mSystemFrequency(50.0), mSystemOmega(2.0 * PI * mSystemFrequency),
        mGridVoltageRMSLineToLine(400.0), mLineResistance(0.3),
        mLineInductance(0.1e-3), mLineCapacitance(1e-6), mLineConductance(1e-6),
        mLoadActivePower(10000.0), mSwitchOpenResistance(1e12),
        mSwitchClosedResistance(1e-6), mLf(2e-3), mCf(10e-6), mRf(0.2),
        mRc(0.2), mKpPLL(0.25), mKiPLL(0.2), mOmegaCutoff(mSystemOmega),
        mPref(10000.0), mQref(5000.0), mKpPowerCtrl(0.05), mKiPowerCtrl(0.2),
        mKpCurrCtrl(0.25), mKiCurrCtrl(1.0) {}

  void run() const {
    const String simNameBase = "EMT_Ph3_SSN_GFL";

    const auto systemPF = runPowerFlow(simNameBase + "_PF");
    runEMTSimulation(simNameBase + "_EMT", systemPF);
  }

private:
  SystemTopology runPowerFlow(const String &simName) const {
    Logger::setLogDir("logs/" + simName);

    auto nGrid = SimNode<Complex>::make("nGrid", PhaseType::Single);
    auto nPcc = SimNode<Complex>::make("nPcc", PhaseType::Single);

    auto slack = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::info);
    slack->setParameters(mGridVoltageRMSLineToLine);
    slack->setBaseVoltage(mGridVoltageRMSLineToLine);
    slack->modifyPowerFlowBusType(PowerflowBusType::VD);

    auto line = SP::Ph1::PiLine::make("Line", Logger::Level::info);
    line->setParameters(mLineResistance, mLineInductance, mLineCapacitance,
                        mLineConductance);
    line->setBaseVoltage(mGridVoltageRMSLineToLine);

    // SSN_GFL uses the same controller equations as the reference inverter.
    // It therefore regulates filter-side power at Vc through Rc. Convert the
    // requested filter-side reference to the corresponding PCC-side injection
    // used for power-flow initialization.
    //
    // The switchable EMT load is intentionally excluded from the initial
    // power-flow operating point.
    const auto [pPccRef, qPccRef] =
        pccPowerFromFilterPowerReference(mPref, mQref);

    // Negative load represents injected active and reactive power at the PCC.
    //
    // The component name is identical to the EMT inverter name so
    // initWithPowerflow() can transfer the operating point.
    auto inverterInjection =
        SP::Ph1::Load::make("INV_SSN_GFL", Logger::Level::info);
    inverterInjection->setParameters(-pPccRef, -qPccRef,
                                     mGridVoltageRMSLineToLine);
    inverterInjection->modifyPowerFlowBusType(PowerflowBusType::PQ);

    slack->connect({nGrid});
    line->connect({nGrid, nPcc});
    inverterInjection->connect({nPcc});

    auto system =
        SystemTopology(mSystemFrequency, SystemNodeList{nGrid, nPcc},
                       SystemComponentList{slack, line, inverterInjection});

    auto logger = DataLogger::make(simName);
    logger->logAttribute("v_grid_pf", nGrid->attribute("v"));
    logger->logAttribute("v_pcc_pf", nPcc->attribute("v"));

    const Real timeStepPF = mFinalTimeEMT;
    const Real finalTimePF = mFinalTimeEMT + timeStepPF;

    Simulation sim(simName, Logger::Level::info);
    sim.setSystem(system);
    sim.setTimeStep(timeStepPF);
    sim.setFinalTime(finalTimePF);
    sim.setDomain(Domain::SP);
    sim.setSolverType(Solver::Type::NRP);
    sim.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
    sim.doInitFromNodesAndTerminals(false);
    sim.addLogger(logger);
    sim.run();

    return system;
  }

  void runEMTSimulation(const String &simName,
                        const SystemTopology &systemPF) const {
    Logger::setLogDir("logs/" + simName);

    auto nGrid = SimNode<Real>::make("nGrid", PhaseType::ABC);
    auto nPcc = SimNode<Real>::make("nPcc", PhaseType::ABC);
    auto nLoad = SimNode<Real>::make("nLoad", PhaseType::ABC);

    auto slack = EMT::Ph3::NetworkInjection::make("Slack", Logger::Level::info);

    auto line = EMT::Ph3::PiLine::make("Line", Logger::Level::info);
    line->setParameters(
        Math::singlePhaseParameterToThreePhase(mLineResistance),
        Math::singlePhaseParameterToThreePhase(mLineInductance),
        Math::singlePhaseParameterToThreePhase(mLineCapacitance),
        Math::singlePhaseParameterToThreePhase(mLineConductance));

    // Refactored GFL component. The setParameters() signature and all parameter
    // meanings intentionally match AvVoltSourceInverterStateSpace.
    auto inverter = EMT::Ph3::SSN_GFL::make("INV_SSN_GFL", Logger::Level::info);
    inverter->setParameters(mLf, mCf, mRf, mRc, mSystemOmega, mKpPLL, mKiPLL,
                            mOmegaCutoff, mPref, mQref, mKpPowerCtrl,
                            mKiPowerCtrl, mKpCurrCtrl, mKiCurrCtrl);

    auto loadSwitch = EMT::Ph3::Switch::make("LoadSwitch", Logger::Level::info);
    loadSwitch->setParameters(
        Math::singlePhaseParameterToThreePhase(mSwitchOpenResistance),
        Math::singlePhaseParameterToThreePhase(mSwitchClosedResistance), false);

    auto loadResistor =
        EMT::Ph3::Resistor::make("LoadResistor", Logger::Level::info);
    loadResistor->setParameters(
        Math::singlePhaseParameterToThreePhase(loadResistance()));

    slack->connect({nGrid});
    line->connect({nGrid, nPcc});

    // terminal 0 = GND, terminal 1 = nPcc
    inverter->connect({EMT::SimNode::GND, nPcc});

    // The resistor load is connected to ground and initially isolated from
    // the PCC by an open switch. Closing the switch at t = 3 s creates the
    // same disturbance as in the reference example.
    loadSwitch->connect({nPcc, nLoad});
    loadResistor->connect({nLoad, EMT::SimNode::GND});

    auto system = SystemTopology(
        mSystemFrequency, SystemNodeList{nGrid, nPcc, nLoad},
        SystemComponentList{slack, line, inverter, loadSwitch, loadResistor});

    system.initWithPowerflow(systemPF, Domain::EMT);

    auto logger = DataLogger::make(simName);

    // Keep the same signal names as the reference example. This allows direct
    // comparison or overlay of the two CSV result sets.
    logger->logAttribute("v_pcc", nPcc->attribute("v"));
    logger->logAttribute("i_inv", inverter->attribute("i_intf"));
    logger->logAttribute("i_load", loadResistor->attribute("i_intf"));
    logger->logAttribute("vc_d", inverter->attribute("vc_d"));
    logger->logAttribute("vc_q", inverter->attribute("vc_q"));
    logger->logAttribute("p_inst", inverter->attribute("p_inst"));
    logger->logAttribute("q_inst", inverter->attribute("q_inst"));
    logger->logAttribute("omega_pll", inverter->attribute("omega_pll"));

    Simulation sim(simName, Logger::Level::info);
    sim.setSystem(system);
    sim.addLogger(logger);
    sim.setDomain(Domain::EMT);
    sim.setSolverType(Solver::Type::MNA);

    // SSN_GFL is locally linearized at every time step, exactly as required by
    // the time-varying dq/abc transformation and controller operating point.
    sim.doSystemMatrixRecomputation(true);
    sim.doInitFromNodesAndTerminals(true);

    sim.addEvent(SwitchEvent3Ph::make(mLoadSwitchTime, loadSwitch, true));

    sim.setTimeStep(mTimeStepEMT);
    sim.setFinalTime(mFinalTimeEMT);
    sim.run();
  }

  Real loadResistance() const {
    return mGridVoltageRMSLineToLine * mGridVoltageRMSLineToLine /
           mLoadActivePower;
  }

  std::pair<Real, Real>
  pccPowerFromFilterPowerReference(Real pFilterRef, Real qFilterRef) const {
    const Real vPccPeakPhase = RMS3PH_TO_PEAK1PH * mGridVoltageRMSLineToLine;

    if (std::abs(mRc) < 1e-12 || vPccPeakPhase < 1e-9)
      return {pFilterRef, qFilterRef};

    // With peak-valued phase phasors:
    //
    // P_filter =
    //     P_pcc + Rc / (1.5 * |U|^2) * (P_pcc^2 + Q_pcc^2)
    //
    // Q_filter = Q_pcc
    const Real qPccRef = qFilterRef;
    const Real a = mRc / (1.5 * vPccPeakPhase * vPccPeakPhase);

    const Real discriminant =
        1.0 + 4.0 * a * (pFilterRef - a * qPccRef * qPccRef);

    if (discriminant < 0.0) {
      throw std::runtime_error(
          "No feasible PCC power found for the requested filter-side power "
          "reference, Rc, and PCC voltage estimate.");
    }

    const Real sqrtDisc = std::sqrt(discriminant);
    const Real p1 = (-1.0 + sqrtDisc) / (2.0 * a);
    const Real p2 = (-1.0 - sqrtDisc) / (2.0 * a);

    // Select the root close to pFilterRef. This is the physically relevant
    // branch for a small coupling resistance.
    const Real pPccRef =
        std::abs(p1 - pFilterRef) < std::abs(p2 - pFilterRef) ? p1 : p2;

    return {pPccRef, qPccRef};
  }

private:
  Real mTimeStepEMT;
  Real mFinalTimeEMT;
  Real mLoadSwitchTime;

  Real mSystemFrequency;
  Real mSystemOmega;
  Real mGridVoltageRMSLineToLine;

  Real mLineResistance;
  Real mLineInductance;
  Real mLineCapacitance;
  Real mLineConductance;

  Real mLoadActivePower;
  Real mSwitchOpenResistance;
  Real mSwitchClosedResistance;

  Real mLf;
  Real mCf;
  Real mRf;
  Real mRc;

  Real mKpPLL;
  Real mKiPLL;

  Real mOmegaCutoff;
  Real mPref;
  Real mQref;
  Real mKpPowerCtrl;
  Real mKiPowerCtrl;
  Real mKpCurrCtrl;
  Real mKiCurrCtrl;
};

int main(int argc, char *argv[]) {
  Example_SSN_GFL example;
  example.run();
  return 0;
}
