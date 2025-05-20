#define USE_MODEL_4TH_ORDER // Select the generator model to use

#include "../GeneratorFactory.h"
#include "../Examples.h"

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

CPS::CIM::Examples::NineBus::ScenarioConfig ninebus;

void IEEE9Bus_Simulation(String simName, Real timeStep,
                                      Real finalTime) {
  // ----- POWER FLOW FOR INITIALIZATION -----
  String simNamePF = simName + "_PF";
  Logger::setLogDir("logs/" + simNamePF);

  // ----- INITIALIZE COMPONENTS -----

  // Nodes
  auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
  auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);
  auto n3PF = SimNode<Complex>::make("n3", PhaseType::Single);
  auto n4PF = SimNode<Complex>::make("n4", PhaseType::Single);
  auto n5PF = SimNode<Complex>::make("n5", PhaseType::Single);
  auto n6PF = SimNode<Complex>::make("n6", PhaseType::Single);
  auto n7PF = SimNode<Complex>::make("n7", PhaseType::Single);
  auto n8PF = SimNode<Complex>::make("n8", PhaseType::Single);
  auto n9PF = SimNode<Complex>::make("n9", PhaseType::Single);

  auto gen1PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen1.Name, Logger::Level::debug);
  gen1PF->setParameters(ninebus.gen1.RatedPower, ninebus.gen1.RatedVoltage,
                        ninebus.gen1.InitialPower, ninebus.gen1.InitialVoltage,
                        ninebus.gen1.BusType); // Reactive power not given, should be 27e6
  gen1PF->setBaseVoltage(ninebus.gen1.RatedVoltage);

  auto gen2PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen2.Name, Logger::Level::debug);
  gen2PF->setParameters(ninebus.gen2.RatedPower, ninebus.gen2.RatedVoltage,
                        ninebus.gen2.InitialPower, ninebus.gen2.InitialVoltage,
                        ninebus.gen2.BusType); // Reactive power not given, should be 6.7e6
  gen2PF->setBaseVoltage(ninebus.gen2.RatedVoltage);

  auto gen3PF =
      SP::Ph1::SynchronGenerator::make(ninebus.gen3.Name, Logger::Level::debug);
  gen3PF->setParameters(ninebus.gen3.RatedPower, ninebus.gen3.RatedVoltage,
                        ninebus.gen3.InitialPower, ninebus.gen3.InitialVoltage,
                        ninebus.gen3.BusType); // Reactive power not given, should be -10.9e6
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

  auto line45PF =
      SP::Ph1::PiLine::make(ninebus.line45.Name, Logger::Level::debug);
  line45PF->setParameters(ninebus.line45.Resistance, ninebus.line45.Inductance,
                          ninebus.line45.Capacitance,
                          ninebus.line45.Conductance);
  line45PF->setBaseVoltage(ninebus.line45.BaseVoltage);

  auto line46PF =
      SP::Ph1::PiLine::make(ninebus.line46.Name, Logger::Level::debug);
  line46PF->setParameters(ninebus.line46.Resistance, ninebus.line46.Inductance,
                          ninebus.line46.Capacitance,
                          ninebus.line46.Conductance);
  line46PF->setBaseVoltage(ninebus.line46.BaseVoltage);

  auto line57PF =
      SP::Ph1::PiLine::make(ninebus.line57.Name, Logger::Level::debug);
  line57PF->setParameters(ninebus.line57.Resistance, ninebus.line57.Inductance,
                          ninebus.line57.Capacitance,
                          ninebus.line57.Conductance);
  line57PF->setBaseVoltage(ninebus.line57.BaseVoltage);

  auto line69PF =
      SP::Ph1::PiLine::make(ninebus.line69.Name, Logger::Level::debug);
  line69PF->setParameters(ninebus.line69.Resistance, ninebus.line69.Inductance,
                          ninebus.line69.Capacitance,
                          ninebus.line69.Conductance);
  line69PF->setBaseVoltage(ninebus.line69.BaseVoltage);

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
      ninebus.transf14.Resistance, ninebus.transf14.Inductance); // RatioAbs should be 1 from RTDS
  transf14PF->setBaseVoltage(ninebus.transf14.VoltageHVSide);

  // Transformer between bus 2 and bus 7
  auto transf27PF =
      SP::Ph1::Transformer::make(ninebus.transf27.Name, Logger::Level::debug);
  transf27PF->setParameters(
      ninebus.transf27.VoltageLVSide, ninebus.transf27.VoltageHVSide,
      ninebus.transf27.Ratio, 0.0, ninebus.transf27.Resistance,
      ninebus.transf27.Inductance); // RatioAbs should be 1 from RTDS
  transf27PF->setBaseVoltage(ninebus.transf27.VoltageHVSide);

  // Transformer between bus 3 and bus 9
  auto transf39PF =
      SP::Ph1::Transformer::make(ninebus.transf39.Name, Logger::Level::debug);
  transf39PF->setParameters(
      ninebus.transf39.VoltageLVSide, ninebus.transf39.VoltageHVSide,
      ninebus.transf39.Ratio, 0.0, ninebus.transf39.Resistance,
      ninebus.transf39.Inductance); // RatioAbs should be 1 from RTDS
  transf39PF->setBaseVoltage(ninebus.transf39.VoltageHVSide);

  // ----- CONNECT COMPONENTS TO SYSTEM TOPOLOGY -----

  gen1PF->connect({n1PF});
  gen2PF->connect({n2PF});
  gen3PF->connect({n3PF});

  load5PF->connect({n5PF});
  load6PF->connect({n6PF});
  load8PF->connect({n8PF});

  line45PF->connect({n4PF, n5PF});
  line46PF->connect({n4PF, n6PF});
  line57PF->connect({n5PF, n7PF});
  line69PF->connect({n6PF, n9PF});
  line78PF->connect({n7PF, n8PF});
  line89PF->connect({n8PF, n9PF});

  transf14PF->connect({n1PF, n4PF});
  transf27PF->connect({n2PF, n7PF});
  transf39PF->connect({n3PF, n9PF});

  // ----- CREATE SYSTEM TOPOLOGY -----
  auto systemPF = SystemTopology(
      ninebus.nomFreq, SystemNodeList{n1PF, n2PF, n3PF, n4PF, n5PF, n6PF, n7PF, n8PF, n9PF},
      SystemComponentList{gen1PF, gen2PF, gen3PF, load5PF, load6PF, load8PF,
                          line45PF, line46PF, line57PF, line69PF, line78PF,
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
//   simPF.doInitFromNodesAndTerminals(false);
  simPF.addLogger(loggerPF);
  simPF.run();


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


// ----- DYNAMIC SIMULATION -----

String simNameEMT = simName + "_EMT";
Logger::setLogDir("logs/" + simNameEMT);

// Nodes
auto n1EMT = SimNode<Real>::make("n1", PhaseType::ABC);
auto n2EMT = SimNode<Real>::make("n2", PhaseType::ABC);
auto n3EMT = SimNode<Real>::make("n3", PhaseType::ABC);
auto n4EMT = SimNode<Real>::make("n4", PhaseType::ABC);
auto n5EMT = SimNode<Real>::make("n5", PhaseType::ABC);
auto n6EMT = SimNode<Real>::make("n6", PhaseType::ABC);
auto n7EMT = SimNode<Real>::make("n7", PhaseType::ABC);
auto n8EMT = SimNode<Real>::make("n8", PhaseType::ABC);
auto n9EMT = SimNode<Real>::make("n9", PhaseType::ABC);

// Generator 1 Initialization
auto gen1EMT = EMT::Ph3::SynchronGenerator4OrderVBR::make(ninebus.gen1.Name, Logger::Level::debug);
gen1EMT->setOperationalParametersPerUnit(
    ninebus.gen1.RatedPower,        // nomPower [VA]
    ninebus.gen1.RatedVoltage,      // nomVolt [V]
    ninebus.nomFreq,                // nomFreq [Hz]
    ninebus.gen1.H,
    ninebus.gen1.Ld,
    ninebus.gen1.Lq,
    ninebus.gen1.L0,
    ninebus.gen1.LdPrime,
    ninebus.gen1.LqPrime,
    ninebus.gen1.TdoPrime,
    ninebus.gen1.TqoPrime
);
// Get actual active and reactive power of generator's Terminal from Powerflow solution
  Complex initApparentPower_G1 = gen1PF->getApparentPower();
  gen1EMT->setInitialValues(initApparentPower_G1, initApparentPower_G1.real(), n1PF->voltage()(0, 0));
  gen1EMT->setModelAsNortonSource(true);

// Generator 2 Initialization
auto gen2EMT = GeneratorFactory::createGenEMT("4", ninebus.gen2.Name, Logger::Level::debug);
gen2EMT->setOperationalParametersPerUnit(
    ninebus.gen2.RatedPower,
    ninebus.gen2.RatedVoltage,
    ninebus.nomFreq,
    ninebus.gen2.H,
    ninebus.gen2.Ld,
    ninebus.gen2.Lq,
    ninebus.gen2.L0,
    ninebus.gen2.LdPrime,
    ninebus.gen2.LqPrime,
    ninebus.gen2.TdoPrime,
    ninebus.gen2.TqoPrime
);
// Get actual active and reactive power of generator's Terminal from Powerflow solution
  Complex initApparentPower_G2 = gen2PF->getApparentPower();
  gen2EMT->setInitialValues(initApparentPower_G2, initApparentPower_G2.real(), n2PF->voltage()(0, 0));
  gen2EMT->setModelAsNortonSource(true);

  // Generator 3 Initialization
auto gen3EMT = GeneratorFactory::createGenEMT("4", ninebus.gen3.Name, Logger::Level::debug);
gen3EMT->setOperationalParametersPerUnit(
    ninebus.gen3.RatedPower,
    ninebus.gen3.RatedVoltage,
    ninebus.nomFreq,
    ninebus.gen3.H,
    ninebus.gen3.Ld,
    ninebus.gen3.Lq,
    ninebus.gen3.L0,
    ninebus.gen3.LdPrime,
    ninebus.gen3.LqPrime,
    ninebus.gen3.TdoPrime,
    ninebus.gen3.TqoPrime
);
// Get actual active and reactive power of generator's Terminal from Powerflow solution
  Complex initApparentPower_G3 = gen3PF->getApparentPower();
  gen3EMT->setInitialValues(initApparentPower_G3, initApparentPower_G3.real(), n3PF->voltage()(0, 0));
  gen3EMT->setModelAsNortonSource(true);


// Loads
auto load5EMT = EMT::Ph3::RXLoad::make(ninebus.load5.Name, Logger::Level::debug);
load5EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.load5.RealPower),
    Math::singlePhaseParameterToThreePhase(ninebus.load5.ReactivePower),
    ninebus.load5.BaseVoltage);

auto load6EMT = EMT::Ph3::RXLoad::make(ninebus.load6.Name, Logger::Level::debug);
load6EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.load6.RealPower),
    Math::singlePhaseParameterToThreePhase(ninebus.load6.ReactivePower),
    ninebus.load6.BaseVoltage);

auto load8EMT = EMT::Ph3::RXLoad::make(ninebus.load8.Name, Logger::Level::debug);
load8EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.load8.RealPower),
    Math::singlePhaseParameterToThreePhase(ninebus.load8.ReactivePower),
    ninebus.load8.BaseVoltage);

// Lines
auto line45EMT = EMT::Ph3::PiLine::make(ninebus.line45.Name, Logger::Level::debug);
line45EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line45.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line45.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line45.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line45.Conductance));

auto line46EMT = EMT::Ph3::PiLine::make(ninebus.line46.Name, Logger::Level::debug);
line46EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line46.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line46.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line46.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line46.Conductance));

auto line57EMT = EMT::Ph3::PiLine::make(ninebus.line57.Name, Logger::Level::debug);
line57EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line57.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line57.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line57.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line57.Conductance));

auto line69EMT = EMT::Ph3::PiLine::make(ninebus.line69.Name, Logger::Level::debug);
line69EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line69.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line69.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line69.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line69.Conductance));

auto line78EMT = EMT::Ph3::PiLine::make(ninebus.line78.Name, Logger::Level::debug);
line78EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line78.Conductance));

auto line89EMT = EMT::Ph3::PiLine::make(ninebus.line89.Name, Logger::Level::debug);
line89EMT->setParameters(
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Inductance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Capacitance),
    Math::singlePhaseParameterToThreePhase(ninebus.line89.Conductance));

// Transformers
// Check which is the high voltage side and which is the low voltage side
auto transf14EMT = EMT::Ph3::Transformer::make(ninebus.transf14.Name, Logger::Level::debug);
transf14EMT->setParameters(
    ninebus.transf14.VoltageLVSide,
    ninebus.transf14.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf14.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf14.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf14.Inductance));

auto transf27EMT = EMT::Ph3::Transformer::make(ninebus.transf27.Name, Logger::Level::debug);
transf27EMT->setParameters(
    ninebus.transf27.VoltageLVSide,
    ninebus.transf27.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf27.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf27.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf27.Inductance));

auto transf39EMT = EMT::Ph3::Transformer::make(ninebus.transf39.Name, Logger::Level::debug);
transf39EMT->setParameters(
    ninebus.transf39.VoltageLVSide,
    ninebus.transf39.VoltageHVSide,
    ninebus.transf14.RatedPower,
    ninebus.transf39.Ratio,
    0.0,
    Math::singlePhaseParameterToThreePhase(ninebus.transf39.Resistance),
    Math::singlePhaseParameterToThreePhase(ninebus.transf39.Inductance));

// Connect components to nodes
gen1EMT->connect({n1EMT});
gen2EMT->connect({n2EMT});
gen3EMT->connect({n3EMT});

load5EMT->connect({n5EMT});
load6EMT->connect({n6EMT});
load8EMT->connect({n8EMT});

line45EMT->connect({n4EMT, n5EMT});
line46EMT->connect({n4EMT, n6EMT});
line57EMT->connect({n5EMT, n7EMT});
line69EMT->connect({n6EMT, n9EMT});
line78EMT->connect({n7EMT, n8EMT});
line89EMT->connect({n8EMT, n9EMT});

transf14EMT->connect({n1EMT, n4EMT});
transf27EMT->connect({n2EMT, n7EMT});
transf39EMT->connect({n3EMT, n9EMT});

// Create system topology
auto systemEMT = SystemTopology(
    ninebus.nomFreq,  // System frequency in Hz
    SystemNodeList{n1EMT, n2EMT, n3EMT, n4EMT, n5EMT, n6EMT, n7EMT, n8EMT, n9EMT},
    SystemComponentList{gen1EMT, gen2EMT, gen3EMT, load5EMT, load6EMT, load8EMT,
                        line45EMT, line46EMT, line57EMT, line69EMT, line78EMT,
                        line89EMT, transf14EMT, transf27EMT, transf39EMT});


// Loggin
auto loggerEMT = DataLogger::make(simNameEMT, Logger::Level::debug);

// Log node voltages
loggerEMT->logAttribute("v1", n1EMT->attribute("v"));
loggerEMT->logAttribute("v2", n2EMT->attribute("v"));
loggerEMT->logAttribute("v3", n3EMT->attribute("v"));
loggerEMT->logAttribute("v4", n4EMT->attribute("v"));
loggerEMT->logAttribute("v5", n5EMT->attribute("v"));
loggerEMT->logAttribute("v6", n6EMT->attribute("v"));
loggerEMT->logAttribute("v7", n7EMT->attribute("v"));
loggerEMT->logAttribute("v8", n8EMT->attribute("v"));
loggerEMT->logAttribute("v9", n9EMT->attribute("v"));


// Log transformer voltages and currents
loggerEMT->logAttribute("v_transf14", transf14EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_transf14", transf14EMT->attribute("i_intf"));
loggerEMT->logAttribute("v_transf27", transf27EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_transf27", transf27EMT->attribute("i_intf"));
loggerEMT->logAttribute("v_transf39", transf39EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_transf39", transf39EMT->attribute("i_intf"));

// Log interface voltages and currents (abc frame)
loggerEMT->logAttribute("v_gen1", gen1EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_gen1", gen1EMT->attribute("i_intf"));
loggerEMT->logAttribute("v_gen2", gen2EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_gen2", gen2EMT->attribute("i_intf"));
loggerEMT->logAttribute("v_gen3", gen3EMT->attribute("v_intf"));
loggerEMT->logAttribute("i_gen3", gen3EMT->attribute("i_intf"));



// Simulation setup and run
systemEMT.initWithPowerflow(systemPF, Domain::EMT);
Simulation simEMT(simNameEMT, Logger::Level::debug);
simEMT.setTimeStep(timeStep);
simEMT.setFinalTime(finalTime);
simEMT.setDomain(Domain::EMT);
simEMT.addLogger(loggerEMT);
simEMT.setSystem(systemEMT);
simEMT.doSystemMatrixRecomputation(true);

simEMT.run();
}

int main(int argc, char *argv[]) {
    // Default simulation parameters
    String simName = "Cosim-IEEE9";
    double finalTime = 1.0; // seconds (can be adjusted)
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
