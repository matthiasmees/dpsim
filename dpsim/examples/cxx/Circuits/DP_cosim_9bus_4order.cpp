#include "../GeneratorFactory.h"
#include "../GridParameters.h"

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

CPS::CIM::Examples::NineBus::ScenarioConfig ninebus;
CPS::CIM::Examples::Components::GovernorKundur::Parameters govKundur;
CPS::CIM::Examples::Components::ExcitationSystemEremia::Parameters excEremia;

void IEEE9Bus_Simulation(String simName, Real timeStep,
                                      Real finalTime) {
  // ----- POWER FLOW FOR INITIALIZATION -----
  String simNamePF = simName + "_PF";
  Logger::setLogDir("logs/" + simNamePF);

  // ----- INITIALIZE COMPONENTS -----

  // Nodes
  auto n1PF = SimNode<Complex>::make("BUS1", PhaseType::Single);
  auto n2PF = SimNode<Complex>::make("BUS2", PhaseType::Single);
  auto n3PF = SimNode<Complex>::make("BUS3", PhaseType::Single);
  auto n4PF = SimNode<Complex>::make("BUS4", PhaseType::Single);
  auto n5PF = SimNode<Complex>::make("BUS5", PhaseType::Single);
  auto n6PF = SimNode<Complex>::make("BUS6", PhaseType::Single);
  auto n7PF = SimNode<Complex>::make("BUS7", PhaseType::Single);
  auto n8PF = SimNode<Complex>::make("BUS8", PhaseType::Single);
  auto n9PF = SimNode<Complex>::make("BUS9", PhaseType::Single);

  auto gen1PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen1.Name, Logger::Level::debug);
  gen1PF->setParameters(ninebus.gen1.RatedPower, ninebus.gen1.RatedVoltage,
                        ninebus.gen1.InitialPower, ninebus.gen1.InitialVoltage,
                        ninebus.gen1.BusType);
  gen1PF->setBaseVoltage(ninebus.gen1.RatedVoltage);

  auto gen2PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen2.Name, Logger::Level::debug);
  gen2PF->setParameters(ninebus.gen2.RatedPower, ninebus.gen2.RatedVoltage,
                        ninebus.gen2.InitialPower, ninebus.gen2.InitialVoltage,
                        ninebus.gen2.BusType);
  gen2PF->setBaseVoltage(ninebus.gen2.RatedVoltage);

  auto gen3PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen3.Name, Logger::Level::debug);
  gen3PF->setParameters(ninebus.gen3.RatedPower, ninebus.gen3.RatedVoltage,
                        ninebus.gen3.InitialPower, ninebus.gen3.InitialVoltage,
                        ninebus.gen3.BusType);
  gen3PF->setBaseVoltage(ninebus.gen3.RatedVoltage);

  // Loads
  auto load5PF = SP::Ph1::Load::make(ninebus.load5.Name, Logger::Level::debug);
  load5PF->setParameters(ninebus.load5.RealPower, ninebus.load5.ReactivePower,
                         ninebus.load5.BaseVoltage);
  load5PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load6PF = SP::Ph1::Load::make(ninebus.load6.Name, Logger::Level::debug);
  load6PF->setParameters(ninebus.load6.RealPower, ninebus.load6.ReactivePower,
                         ninebus.load6.BaseVoltage);
  load6PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load8PF = SP::Ph1::Load::make(ninebus.load8.Name, Logger::Level::debug);
  load8PF->setParameters(ninebus.load8.RealPower, ninebus.load8.ReactivePower,
                         ninebus.load8.BaseVoltage);
  load8PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // Transmission Lines

  auto line54PF =
      SP::Ph1::PiLine::make(ninebus.line54.Name, Logger::Level::debug);
  line54PF->setParameters(ninebus.line54.Resistance, ninebus.line54.Inductance,
                          ninebus.line54.Capacitance,
                          ninebus.line54.Conductance);
  line54PF->setBaseVoltage(ninebus.line54.BaseVoltage);

  auto line64PF =
      SP::Ph1::PiLine::make(ninebus.line64.Name, Logger::Level::debug);
  line64PF->setParameters(ninebus.line64.Resistance, ninebus.line64.Inductance,
                          ninebus.line64.Capacitance,
                          ninebus.line64.Conductance);
  line64PF->setBaseVoltage(ninebus.line64.BaseVoltage);

  auto line75PF =
      SP::Ph1::PiLine::make(ninebus.line75.Name, Logger::Level::debug);
  line75PF->setParameters(ninebus.line75.Resistance, ninebus.line75.Inductance,
                          ninebus.line75.Capacitance,
                          ninebus.line75.Conductance);
  line75PF->setBaseVoltage(ninebus.line75.BaseVoltage);

  auto line96PF =
      SP::Ph1::PiLine::make(ninebus.line96.Name, Logger::Level::debug);
  line96PF->setParameters(ninebus.line96.Resistance, ninebus.line96.Inductance,
                          ninebus.line96.Capacitance,
                          ninebus.line96.Conductance);
  line96PF->setBaseVoltage(ninebus.line96.BaseVoltage);

  auto line78PF =
      SP::Ph1::PiLine::make(ninebus.line78.Name, Logger::Level::debug);
  line78PF->setParameters(ninebus.line78.Resistance, ninebus.line78.Inductance,
                          ninebus.line78.Capacitance,
                          ninebus.line78.Conductance);
  line78PF->setBaseVoltage(ninebus.line78.BaseVoltage);

  auto line89PF =
      SP::Ph1::PiLine::make(ninebus.line89.Name, Logger::Level::debug);
  line89PF->setParameters(ninebus.line89.Resistance, ninebus.line89.Inductance,
                          ninebus.line89.Capacitance,
                          ninebus.line89.Conductance);
  line89PF->setBaseVoltage(ninebus.line89.BaseVoltage);

  // Transformers

  // Transformer between bus 1 and bus 4
  auto transf14PF =
      SP::Ph1::Transformer::make(ninebus.transf14.Name, Logger::Level::debug);
  transf14PF->setParameters(
      ninebus.transf14.VoltageLVSide, ninebus.transf14.VoltageHVSide,
      ninebus.transf14.Ratio, 0.0, // No phase shift (ratioPhase = 0.0)
      ninebus.transf14.Resistance, ninebus.transf14.Inductance);
  transf14PF->setBaseVoltage(ninebus.transf14.VoltageHVSide);

  // Transformer between bus 2 and bus 7
  auto transf27PF =
      SP::Ph1::Transformer::make(ninebus.transf27.Name, Logger::Level::debug);
  transf27PF->setParameters(
      ninebus.transf27.VoltageLVSide, ninebus.transf27.VoltageHVSide,
      ninebus.transf27.Ratio, 0.0, ninebus.transf27.Resistance,
      ninebus.transf27.Inductance);
  transf27PF->setBaseVoltage(ninebus.transf27.VoltageHVSide);

  // Transformer between bus 3 and bus 9
  auto transf39PF =
      SP::Ph1::Transformer::make(ninebus.transf39.Name, Logger::Level::debug);
  transf39PF->setParameters(
      ninebus.transf39.VoltageLVSide, ninebus.transf39.VoltageHVSide,
      ninebus.transf39.Ratio, 0.0, ninebus.transf39.Resistance,
      ninebus.transf39.Inductance);
  transf39PF->setBaseVoltage(ninebus.transf39.VoltageHVSide);

  // ----- CONNECT COMPONENTS TO SYSTEM TOPOLOGY -----

  gen1PF->connect({n1PF});
  gen2PF->connect({n2PF});
  gen3PF->connect({n3PF});

  load5PF->connect({n5PF});
  load6PF->connect({n6PF});
  load8PF->connect({n8PF});

  line54PF->connect({n5PF, n4PF});
  line64PF->connect({n6PF, n4PF});
  line75PF->connect({n7PF, n5PF});
  line96PF->connect({n9PF, n6PF});
  line78PF->connect({n7PF, n8PF});
  line89PF->connect({n8PF, n9PF});

  transf14PF->connect({n1PF, n4PF});
  transf27PF->connect({n2PF, n7PF});
  transf39PF->connect({n3PF, n9PF});

  // ----- CREATE SYSTEM TOPOLOGY -----
  auto systemPF = SystemTopology(
      ninebus.nomFreq, SystemNodeList{n1PF, n2PF, n3PF, n4PF, n5PF, n6PF, n7PF, n8PF, n9PF},
      SystemComponentList{gen1PF, gen2PF, gen3PF, load5PF, load6PF, load8PF,
                          line54PF, line64PF, line75PF, line96PF, line78PF,
                          line89PF, transf14PF, transf27PF, transf39PF});

  systemPF.renderToFile("logs/"+simNamePF+".svg");

  // ----- LOGGING -----
  auto loggerPF = DataLogger::make(simNamePF, Logger::Level::debug);
  // Log node voltages
  loggerPF->logAttribute("v_bus1", n1PF->attribute("v"));
  loggerPF->logAttribute("v_bus2", n2PF->attribute("v"));
  loggerPF->logAttribute("v_bus3", n3PF->attribute("v"));
  loggerPF->logAttribute("v_bus4", n4PF->attribute("v"));
  loggerPF->logAttribute("v_bus5", n5PF->attribute("v"));
  loggerPF->logAttribute("v_bus6", n6PF->attribute("v"));
  loggerPF->logAttribute("v_bus7", n7PF->attribute("v"));
  loggerPF->logAttribute("v_bus8", n8PF->attribute("v"));
  loggerPF->logAttribute("v_bus9", n9PF->attribute("v"));
  // Log node powers
  loggerPF->logAttribute("s_bus1", n1PF->attribute("s"));
  loggerPF->logAttribute("s_bus2", n2PF->attribute("s"));
  loggerPF->logAttribute("s_bus3", n3PF->attribute("s"));
  loggerPF->logAttribute("s_bus4", n4PF->attribute("s"));
  loggerPF->logAttribute("s_bus5", n5PF->attribute("s"));
  loggerPF->logAttribute("s_bus6", n6PF->attribute("s"));
  loggerPF->logAttribute("s_bus7", n7PF->attribute("s"));
  loggerPF->logAttribute("s_bus8", n8PF->attribute("s"));
  loggerPF->logAttribute("s_bus9", n9PF->attribute("s"));

  // ----- SIMULATION SETUP AND RUN -----
  Simulation simPF(simNamePF, Logger::Level::debug);
  simPF.setSystem(systemPF);
  simPF.setTimeStep(timeStep);
  simPF.setFinalTime(1*timeStep);
  simPF.setDomain(Domain::SP);
  simPF.setSolverType(Solver::Type::NRP);
  simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Simulation);
  simPF.addLogger(loggerPF);
  simPF.run();


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


// ----- DYNAMIC SIMULATION -----

String simNameDP = simName + "_DP";
Logger::setLogDir("logs/" + simNameDP);

// Nodes
auto n1DP = SimNode<Real>::make("BUS1", PhaseType::ABC);
auto n2DP = SimNode<Real>::make("BUS2", PhaseType::ABC);
auto n3DP = SimNode<Real>::make("BUS3", PhaseType::ABC);
auto n4DP = SimNode<Real>::make("BUS4", PhaseType::ABC);
auto n5DP = SimNode<Real>::make("BUS5", PhaseType::ABC);
auto n6DP = SimNode<Real>::make("BUS6", PhaseType::ABC);
auto n7DP = SimNode<Real>::make("BUS7", PhaseType::ABC);
auto n8DP = SimNode<Real>::make("BUS8", PhaseType::ABC);
auto n9DP = SimNode<Real>::make("BUS9", PhaseType::ABC);

// Generator 1 Initialization
auto gen1DP = DP::Ph3::SynchronGenerator4OrderVBR::make(ninebus.gen1.Name, Logger::Level::debug);

gen1DP->setOperationalParametersPerUnit(
    ninebus.gen1.RatedPower,        // nomPower [VA]
    ninebus.gen1.RatedVoltage,      // nomVolt [V]
    ninebus.nomFreq,                // nomFreq [Hz]
    ninebus.gen1.H,
    ninebus.gen1.Xd,
    ninebus.gen1.Xq,
    ninebus.gen1.Xa,
    ninebus.gen1.XdPrime,
    ninebus.gen1.XqPrime,
    ninebus.gen1.TdoPrime,
    ninebus.gen1.TqoPrime
);

// Add Exciter for Generator 1
gen1DP->addExciter(
    ninebus.exc1.TA,      // Ta: Voltage regulator time constant
    ninebus.exc1.KA,      // Ka: Voltage regulator gain
    ninebus.exc1.TE,      // Te: Exciter time constant
    ninebus.exc1.KE,      // Ke: Exciter field proportional constant
    ninebus.exc1.TF,      // Tf: Stabilizer time constant
    ninebus.exc1.KF,      // Kf: Stabilizer gain
    0.01                  // Tr: Measurement time constant (placeholder)
);

 // Extract relevant powerflow results
//CPS::Real initElecPower_G1 = Math::abs(gen1PF->getApparentPower());
CPS::Real T4 = 1.0;
CPS::Real T5 = 1.0;

std::shared_ptr<Signal::TurbineGovernorType1> turbineGovernor1 = Signal::TurbineGovernorType1::make("Gen1_TurbineGovernor", CPS::Logger::Level::debug);

turbineGovernor1->setParameters(
    ninebus.gov1.T2,                            // T3: Transient gain time constant
    T4,                                        // T4: Power fraction time constant
    T5,                                        // T5: Reheat time constant
    ninebus.gov1.T3,                            // Tc: Servo time constant
    ninebus.gov1.T1,                            // Ts: Governor time constant
    ninebus.gov1.R,                             // R: Droop
    ninebus.gov1.Vmin,                          // Pmin: Maximum torque
    ninebus.gov1.Vmax,                          // Pmax: Minimum torque
    1.0// initElecPower_G1/ninebus.gen1.RatedPower   // OmRef: Speed reference
);

gen1DP->addGovernor(turbineGovernor1);

// Generator 2 Initialization
auto gen2DP = DP::Ph3::SynchronGenerator4OrderVBR::make(ninebus.gen2.Name, Logger::Level::debug);

gen2DP->setOperationalParametersPerUnit(
    ninebus.gen2.RatedPower,        // nomPower [VA]
    ninebus.gen2.RatedVoltage,      // nomVolt [V]
    ninebus.nomFreq,                // nomFreq [Hz]
    ninebus.gen2.H,
    ninebus.gen2.Xd,
    ninebus.gen2.Xq,
    ninebus.gen2.Xa,
    ninebus.gen2.XdPrime,
    ninebus.gen2.XqPrime,
    ninebus.gen2.TdoPrime,
    ninebus.gen2.TqoPrime
);
// Add Exciter for Generator 2
gen2DP->addExciter(
    ninebus.exc2.TA,      // Ta: Voltage regulator time constant
    ninebus.exc2.KA,      // Ka: Voltage regulator gain
    ninebus.exc2.TE,      // Te: Exciter time constant
    ninebus.exc2.KE,      // Ke: Exciter field proportional constant
    ninebus.exc2.TF,      // Tf: Stabilizer time constant
    ninebus.exc2.KF,      // Kf: Stabilizer gain
    0.01                  // Tr: Measurement time constant (placeholder)
    //ninebus.exc2.VRmax,
    //ninebus.exc2.VRmax
);


std::shared_ptr<Signal::TurbineGovernorType1> turbineGovernor2 = Signal::TurbineGovernorType1::make("Gen2_TurbineGovernor", CPS::Logger::Level::debug);

turbineGovernor2->setParameters(
    ninebus.gov2.T2,                            // T3: Transient gain time constant
    T4,                                        // T4: Power fraction time constant
    T5,                                        // T5: Reheat time constant
    ninebus.gov2.T3,                            // Tc: Servo time constant
    ninebus.gov2.T1,                            // Ts: Governor time constant
    ninebus.gov2.R,                             // R: Droop
    ninebus.gov2.Vmin,                          // Pmin: Maximum torque
    ninebus.gov2.Vmax,                          // Pmax: Minimum torque
    1.0 //initElecPower_G2/ninebus.gen2.RatedPower   // OmRef: Speed reference
);

gen2DP->addGovernor(turbineGovernor2);

// Generator 3 Initialization
auto gen3DP = DP::Ph3::SynchronGenerator4OrderVBR::make(ninebus.gen3.Name, Logger::Level::debug);

gen3DP->setOperationalParametersPerUnit(
    ninebus.gen3.RatedPower,        // nomPower [VA]
    ninebus.gen3.RatedVoltage,      // nomVolt [V]
    ninebus.nomFreq,                // nomFreq [Hz]
    ninebus.gen3.H,
    ninebus.gen3.Xd,
    ninebus.gen3.Xq,
    ninebus.gen3.Xa,
    ninebus.gen3.XdPrime,
    ninebus.gen3.XqPrime,
    ninebus.gen3.TdoPrime,
    ninebus.gen3.TqoPrime
);

// Add Exciter for Generator 1
gen3DP->addExciter(
    ninebus.exc3.TA,      // Ta: Voltage regulator time constant
    ninebus.exc3.KA,      // Ka: Voltage regulator gain
    ninebus.exc3.TE,      // Te: Exciter time constant
    ninebus.exc3.KE,      // Ke: Exciter field proportional constant
    ninebus.exc3.TF,      // Tf: Stabilizer time constant
    ninebus.exc3.KF,      // Kf: Stabilizer gain
    0.01                  // Tr: Measurement time constant (placeholder)
);


std::shared_ptr<Signal::TurbineGovernorType1> turbineGovernor3 = Signal::TurbineGovernorType1::make("Gen3_TurbineGovernor", CPS::Logger::Level::debug);

turbineGovernor3->setParameters(
    ninebus.gov3.T2,                            // T3: Transient gain time constant
    T4,                                        // T4: Power fraction time constant
    T5,                                        // T5: Reheat time constant
    ninebus.gov3.T3,                            // Tc: Servo time constant
    ninebus.gov3.T1,                            // Ts: Governor time constant
    ninebus.gov3.R,                             // R: Droop
    ninebus.gov3.Vmin,                          // Pmin: Maximum torque
    ninebus.gov3.Vmax,                          // Pmax: Minimum torque
    1.0//initElecPower_G3/ninebus.gen3.RatedPower   // OmRef: Speed reference
);

gen3DP->addGovernor(turbineGovernor3);

// Loads
auto load5DP = DP::Ph3::RXLoad::make(ninebus.load5.Name, Logger::Level::debug);
load5DP->setParameters(
    Math::singlePhasePowerToThreePhase(ninebus.load5.RealPower),
    Math::singlePhasePowerToThreePhase(ninebus.load5.ReactivePower),
    ninebus.load5.BaseVoltage);

auto load6DP = DP::Ph3::RXLoad::make(ninebus.load6.Name, Logger::Level::debug);
load6DP->setParameters(
    Math::singlePhasePowerToThreePhase(ninebus.load6.RealPower),
    Math::singlePhasePowerToThreePhase(ninebus.load6.ReactivePower),
    ninebus.load6.BaseVoltage);

auto load8DP = DP::Ph3::RXLoad::make(ninebus.load8.Name, Logger::Level::debug);
load8DP->setParameters(
    Math::singlePhasePowerToThreePhase(ninebus.load8.RealPower),
    Math::singlePhasePowerToThreePhase(ninebus.load8.ReactivePower),
    ninebus.load8.BaseVoltage);

// Lines
auto line54DP = DP::Ph3::PiLine::make(ninebus.line54.Name, Logger::Level::debug);
line54DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line54.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line54.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line54.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line54.Conductance));

auto line64DP = DP::Ph3::PiLine::make(ninebus.line64.Name, Logger::Level::debug);
line64DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line64.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line64.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line64.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line64.Conductance));

auto line75DP = DP::Ph3::PiLine::make(ninebus.line75.Name, Logger::Level::debug);
line75DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line75.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line75.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line75.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line75.Conductance));

auto line96DP = DP::Ph3::PiLine::make(ninebus.line96.Name, Logger::Level::debug);
line96DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line96.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line96.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line96.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line96.Conductance));

auto line78DP = DP::Ph3::PiLine::make(ninebus.line78.Name, Logger::Level::debug);
line78DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Conductance));

auto line89DP = DP::Ph3::PiLine::make(ninebus.line89.Name, Logger::Level::debug);
line89DP->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Conductance));

// Transformers

auto transf14DP = DP::Ph3::Transformer::make(ninebus.transf14.Name, Logger::Level::debug);
transf14DP->setParameters(
    ninebus.transf14.VoltageLVSide,
    ninebus.transf14.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf14.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf14.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf14.Inductance));

auto transf27DP = DP::Ph3::Transformer::make(ninebus.transf27.Name, Logger::Level::debug);
transf27DP->setParameters(
    ninebus.transf27.VoltageLVSide,
    ninebus.transf27.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf27.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf27.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf27.Inductance));

auto transf39DP = DP::Ph3::Transformer::make(ninebus.transf39.Name, Logger::Level::debug);
transf39DP->setParameters(
    ninebus.transf39.VoltageLVSide,
    ninebus.transf39.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf39.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf39.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf39.Inductance));

// Connect components to nodes
gen1DP->connect({n1DP});
gen2DP->connect({n2DP});
gen3DP->connect({n3DP});

load5DP->connect({n5DP});
load6DP->connect({n6DP});
load8DP->connect({n8DP});

line54DP->connect({n5DP, n4DP});
line64DP->connect({n6DP, n4DP});
line75DP->connect({n7DP, n5DP});
line96DP->connect({n9DP, n6DP});
line78DP->connect({n7DP, n8DP});
line89DP->connect({n8DP, n9DP});

transf14DP->connect({n1DP, n4DP});
transf27DP->connect({n2DP, n7DP});
transf39DP->connect({n3DP, n9DP});

// Create system topology
auto systemDP = SystemTopology(
    ninebus.nomFreq,
    SystemNodeList{n1DP, n2DP, n3DP, n4DP, n5DP, n6DP, n7DP, n8DP, n9DP},
    SystemComponentList{gen1DP, gen2DP, gen3DP, load5DP, load6DP, load8DP,
                        line54DP, line64DP, line75DP, line96DP, line78DP,
                        line89DP, transf14DP, transf27DP, transf39DP});


// Loggin
auto loggerDP = DataLogger::make(simNameDP, Logger::Level::debug);
loggerDP->logAttribute("v1", n1DP->attribute("v"));
loggerDP->logAttribute("v2", n2DP->attribute("v"));
loggerDP->logAttribute("v3", n3DP->attribute("v"));
loggerDP->logAttribute("v4", n4DP->attribute("v"));
loggerDP->logAttribute("v5", n5DP->attribute("v"));
loggerDP->logAttribute("v6", n6DP->attribute("v"));
loggerDP->logAttribute("v7", n7DP->attribute("v"));
loggerDP->logAttribute("v8", n8DP->attribute("v"));
loggerDP->logAttribute("v9", n9DP->attribute("v"));

// log generator's current
for (auto comp : systemDP.mComponents) {
if (std::dynamic_pointer_cast<CPS::DP::Ph3::SynchronGenerator4OrderVBR>(comp)) {
    loggerDP->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
    loggerDP->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
    loggerDP->logAttribute(comp->name() + ".omega", comp->attribute("w_r"));
    loggerDP->logAttribute(comp->name() + ".delta", comp->attribute("delta"));
}
}


// log transfomers voltages & currents
for (auto comp : systemDP.mComponents) {
if (std::dynamic_pointer_cast<CPS::DP::Ph3::Transformer>(comp)) {
    loggerDP->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
    loggerDP->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
}
}

// log Lines voltages & currents
for (auto comp : systemDP.mComponents) {
if (std::dynamic_pointer_cast<CPS::DP::Ph3::PiLine>(comp)) {
    loggerDP->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
    loggerDP->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
}
}

// Simulation setup and run
systemDP.initWithPowerflow(systemPF, Domain::DP);
Simulation simDP(simNameDP, Logger::Level::debug);
simDP.setTimeStep(timeStep);
simDP.setFinalTime(finalTime);
simDP.setDomain(Domain::DP);
simDP.addLogger(loggerDP);
simDP.setSystem(systemDP);
simDP.doSystemMatrixRecomputation(true);
simDP.run();
}

int main(int argc, char *argv[]) {
    // Default simulation parameters
    String simName = "Cosim-IEEE9-4order";
    double finalTime = 5; // seconds (can be adjusted)
    double timeStep = 50e-6;

    // Check for command-line arguments
    if (argc >= 2) {
        simName = argv[1];
    }
    if (argc >= 3) {
        finalTime = atof(argv[2]);
    }
    if (argc >= 4) {
        timeStep = atof(argv[3]);
    }

    // Run the simulation with the given parameters
    IEEE9Bus_Simulation(simName, timeStep, finalTime);

    return 0;
}

