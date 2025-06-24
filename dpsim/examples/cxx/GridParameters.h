/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

 #pragma once

 #include <DPsim.h>

 namespace CPS {
 namespace CIM {
 namespace Examples {

 namespace NineBus {
 struct Gen {
   String Name;              // Name of Generator
   Real RatedPower;          // Rated Power of Generator
   Real RatedVoltage;        // Rated Voltage of Generator
   Real InitialPower;        // Initial Power of Generator
   Real InitialPowerReactive; // Initial Reactive Power of Generator
   Real InitialVoltage;      // Initial Voltage of Generator
   PowerflowBusType BusType; // Bus Type of Generator

   //////////////////////////////////////////////////////////////////
   //////////// Generator Parameters (in per unit) //////////////////
   //////////////////////////////////////////////////////////////////

   Real Ra;                  // Armature Resistance (PU)
   // Note:
   // - Ra (from the datasheet) refers to the armature resistance, a physical property of the generator.
   // - Rs (in setStandardParametersPU) represents the stator resistance, which is equivalent to Ra.
   // Therefore, Ra from the datasheet is passed as Rs in the function for initialization.
   Real Xa;

   Real Xd;                  // Synchronous Reactance (PU)
   Real XdPrime;             // Transient Reactance (PU)
   Real XdDoublePrime;       // Sub-Transient Reactance (PU)
   Real Xq;                  // Quadrature-Axis Reactance (PU)
   Real XqPrime;             // Quadrature-Axis Transient Reactance (PU)
   Real XqDoublePrime;       // Quadrature-Axis Sub-Transient Reactance (PU)
   Real Ld;                  // Synchronous Inductance (PU)
   Real LdPrime;             // Transient Inductance (PU)
   Real LdDoublePrime;       // Sub-Transient Inductance (PU)
   Real Lq;                  // Quadrature-Axis Inductance (PU)
   Real LqPrime;             // Quadrature-Axis Transient Inductance (PU)
   Real LqDoublePrime;       // Quadrature-Axis Sub-Transient Inductance (PU)
   Real H;                   // Inertia Constant (s)
   Real D;                   // Damping Coefficient (PU/PU)
   Real TdoPrime;            // Open-Circuit Field Time Constant (s)
   Real TdoDoublePrime;      // Sub-Transient Time Constant (s)
   Real TqoPrime;
   Real TqoDoublePrime;

   Real Taa;                 // Accelerating Time Constant (s)
   Real L0;                  // Leakage Reactance (PU)
 }; // Generator Structure

 struct Exciter {
   // Exciter Data-1 (IEEE Type DC1A)
   Real KA;                  // Voltage regulator gain (PU)
   Real TA;                  // Voltage regulator time constant (s)
   Real VRmin;               // Minimum exciter output limit (PU)
   Real VRmax;               // Maximum exciter output limit (PU)
   Real KE;                  // Exciter field proportional constant (PU)
   Real TE;                  // Exciter field time constant (s)
   Real KF;                  // Stabilizer feedback gain (PU)
   Real TF;                  // Stabilizer feedback time constant (s)

   // Exciter Data-2 (IEEE Type DC1A)
   Real EX1;                 // Saturation point 1 for exciter (PU)
   Real S_EX1;               // Saturation value at EX1 (PU)
   Real EX2;                 // Saturation point 2 for exciter (PU)
   Real S_EX2;               // Saturation value at EX2 (PU)
   }; // Exciter Structure

   struct Governor {
   // Governor Data (TGOV1)
   Real R;                   // Speed droop (PU/PU)
   Real T1;                  // Governor time constant (s)
   Real Vmax;                // Maximum valve position (PU)
   Real Vmin;                // Minimum valve position (PU)
   Real T2;                  // Turbine time constant (s)
   Real T3;                  // Servo time constant (s)
   Real Dt;                  // Turbine damping coefficient (PU)
   }; // Governor Structure

 struct Load {
   String Name;
   Real RealPower;
   Real ReactivePower;
   Real BaseVoltage;
   Real Conductance;
   Real Susceptance;
 }; // Load Structure

 struct Line {
   String Name;
   Real ResistancePU;  // Resistance per unit
   Real ReactancePU;   // Reactance per unit
   Real SusceptancePU; // Susceptance per unit
   Real Resistance;    // Resistance
   Real Inductance;    // Inductance
   Real Capacitance;   // Capacitance
   Real Conductance;   // Conductance
   Real BaseVoltage;   // Base Voltage
 };                    // Line Structure

 // Define the Transformer struct
 struct Transformer {
   String Name;            // Name of the Transformer
   Real VoltageHVSide;     // High Voltage Side Voltage (in volts)
   Real VoltageLVSide;     // Low Voltage Side Voltage (in volts)
   Real RatedPower;        // Rated Power of the Transformer (in VA)
   Real Ratio;             // Turns Ratio (unitless)
   Real WindingConnection; // Winding Connection Type (0 for default)
   Real ResistancePU;      // Resistance per unit
   Real Resistance;        // Resistance (in ohms)
   Real ReactancePU;       // Reactance  per unit
   Real Inductance;        // Inductance (in henrys)
 };                        // Transformer Structure

 struct ScenarioConfig {
   //-----------Network-----------//
   Real Vnom = 230e3;
   Real nomFreq = 60;
   Real nomOmega = nomFreq * 2 * PI;
   Real baseMVA = 100e6;
   Real baseZ = Vnom * Vnom / baseMVA;

   //-----------Generators-----------//
   Gen gen1;
   Gen gen2;
   Gen gen3;
   //----------Exciters-----------//
   Exciter exc1;
   Exciter exc2;
   Exciter exc3;
   //----------Governors-----------//
   Governor gov1;
   Governor gov2;
   Governor gov3;

   //-----------Loads-----------//
   Load load5;
   Load load6;
   Load load8;

   //-----------Lines-----------//
   Line line54;
   Line line64;
   Line line75;
   Line line96;
   Line line78;
   Line line89;

   //-----------Transformers-----------//
   Transformer transf14;
   Transformer transf27;
   Transformer transf39;

   ScenarioConfig() {
     //-----------------Generator 1 (bus1)-----------------//
     gen1.Name = "GEN1";
     gen1.RatedPower = 100e6;    //71.6e6;   // Rated Power of Generator 1
     gen1.RatedVoltage = 16.5e3; // Rated Voltage of Generator 1
     gen1.InitialPower = 71.6e6; // Initial Power of Generator 1
     gen1.InitialPowerReactive = 27e6; // Initial Reactive Power of Generator 1
     gen1.InitialVoltage = 1.04 * gen1.RatedVoltage; // Initial Voltage of Generator 1
     gen1.BusType = PowerflowBusType::VD;
     gen1.Ra = 0.000125;
     gen1.Xa = 0.01460;
     gen1.L0 = gen1.Xa/(2*PI*nomFreq);

     gen1.Xd = 0.146;
     gen1.XdPrime = 0.0608;
     gen1.XdDoublePrime = 0.06;

     gen1.Ld = gen1.Xd/(2*PI*nomFreq);
     gen1.LdPrime = gen1.XdPrime/(2*PI*nomFreq);
     gen1.LdDoublePrime = gen1.XdDoublePrime/(2*PI*nomFreq);

     gen1.Xq = 0.1;
     gen1.XqPrime = 0.0969;
     gen1.XqDoublePrime = 0.06;

     gen1.Lq = gen1.Xq/(2*PI*nomFreq);
     gen1.LqPrime = gen1.XqPrime/(2*PI*nomFreq);
     gen1.LqDoublePrime = gen1.XqDoublePrime/(2*PI*nomFreq);

     gen1.H = 23.64;
     gen1.D = 0.0;
     gen1.TdoPrime = 8.96;
     gen1.TdoDoublePrime = 0.01;
     gen1.TqoPrime = 0.310;
     gen1.TqoDoublePrime = 0.01;
     gen1.Taa = 0.0;

     // Exciter Data-1 (IEEE Type DC1A)
     exc1.KA = 20.0;
     exc1.TA = 0.2;
     exc1.VRmin = -5.0;
     exc1.VRmax = 5.0;
     exc1.KE = 1.0;
     exc1.TE = 0.314;
     exc1.KF = 0.063;
     exc1.TF = 0.35;

     // Exciter Data-2 (IEEE Type DC1A)
     exc1.EX1 = 3.3;
     exc1.S_EX1 = 0.6602;
     exc1.EX2 = 4.5;
     exc1.S_EX2 = 4.2662;

     // Governor Data (TGOV1)
     gov1.R = 0.05;
     gov1.T1 = 0.05;
     gov1.Vmax = 5.00;
     gov1.Vmin = -5.00;
     gov1.T2 = 2.1;
     gov1.T3 = 7.0;
     gov1.Dt = 0.0;

     //-----------------Generator 2 (bus2)-----------------//
     gen2.Name = "GEN2";
     gen2.RatedPower = 100e6; //163e6;
     gen2.RatedVoltage = 18e3;
     gen2.InitialPower = 163e6;
     gen2.InitialPowerReactive = 6.7e6;
     gen2.InitialVoltage = 1.025 * gen2.RatedVoltage;
     gen2.BusType = PowerflowBusType::PV;
     gen2.Ra = 0.000125;
     gen2.Xa = 0.08958;
     gen2.L0 = gen2.Xa/(2*PI*nomFreq);
     gen2.Xd = 0.8958;
     gen2.XdPrime = 0.1198;
     gen2.XdDoublePrime = 0.11;

     gen2.Ld = gen2.Xd/(2*PI*nomFreq);
     gen2.LdPrime = gen2.XdPrime/(2*PI*nomFreq);
     gen2.LdDoublePrime = gen2.XdDoublePrime/(2*PI*nomFreq);

     gen2.Xq = 0.8645;
     gen2.XqPrime = 0.1969;
     gen2.XqDoublePrime = 0.11;

     gen2.Lq = gen2.Xq/(2*PI*nomFreq);
     gen2.LqPrime = gen2.XqPrime/(2*PI*nomFreq);
     gen2.LqDoublePrime = gen2.XqDoublePrime/(2*PI*nomFreq);

     gen2.H = 6.40;
     gen2.D = 0.0;
     gen2.TdoPrime = 6.00;
     gen2.TdoDoublePrime = 0.01;
     gen2.TqoPrime = 0.535;
     gen2.TqoDoublePrime = 0.01;
     gen2.Taa = 0.0;

     // Exciter Data-1 (IEEE Type DC1A)
     exc2.KA = 20.0;
     exc2.TA = 0.2;
     exc2.VRmin = -5.0;
     exc2.VRmax = 5.0;
     exc2.KE = 1.0;
     exc2.TE = 0.314;
     exc2.KF = 0.063;
     exc2.TF = 0.35;

     // Exciter Data-2 (IEEE Type DC1A)
     exc2.EX1 = 3.3;
     exc2.S_EX1 = 0.6602;
     exc2.EX2 = 4.5;
     exc2.S_EX2 = 4.2662;

     // Governor Data (TGOV1)
     gov2.R = 0.05;
     gov2.T1 = 0.05;
     gov2.Vmax = 5.00;
     gov2.Vmin = -5.00;
     gov2.T2 = 2.1;
     gov2.T3 = 7.0;
     gov2.Dt = 0.0;

     //-----------------Generator 3 (bus3)-----------------//
     gen3.Name = "GEN3";
     gen3.RatedPower = 100e6; //85e6;
     gen3.RatedVoltage = 13.8e3;
     gen3.InitialPower = 85e6;
     gen3.InitialPowerReactive = -10.9e6;
     gen3.InitialVoltage = 1.025 * gen3.RatedVoltage;
     gen3.BusType = PowerflowBusType::PV;
     gen3.Ra = 0.000125;
     gen3.Xa = 0.13125;
     gen3.L0 = gen3.Xa/(2*PI*nomFreq);

     gen3.Xd = 1.3125;
     gen3.XdPrime = 0.1813;
     gen3.XdDoublePrime = 0.18;

     gen3.Ld = gen3.Xd/(2*PI*nomFreq);
     gen3.LdPrime = gen3.XdPrime/(2*PI*nomFreq);
     gen3.LdDoublePrime = gen3.XdDoublePrime/(2*PI*nomFreq);

     gen3.Xq = 1.2578;
     gen3.XqPrime = 0.25;
     gen3.XqDoublePrime = 0.18;

     gen3.Lq = gen3.Xq/(2*PI*nomFreq);
     gen3.LqPrime = gen3.XqPrime/(2*PI*nomFreq);
     gen3.LqDoublePrime = gen3.XqDoublePrime/(2*PI*nomFreq);

     gen3.H = 3.01;
     gen3.D = 0.0;
     gen3.TdoPrime = 5.89;
     gen3.TdoDoublePrime = 0.01;
     gen3.TqoPrime = 0.600;
     gen3.TqoDoublePrime = 0.01;
     gen3.Taa = 0.0;

     // Exciter Data-1 (IEEE Type DC1A)
     exc3.KA = 20.0;
     exc3.TA = 0.2;
     exc3.VRmin = -5.0;
     exc3.VRmax = 5.0;
     exc3.KE = 1.0;
     exc3.TE = 0.314;
     exc3.KF = 0.063;
     exc3.TF = 0.35;

     // Exciter Data-2 (IEEE Type DC1A)
     exc3.EX1 = 3.3;
     exc3.S_EX1 = 0.6602;
     exc3.EX2 = 4.5;
     exc3.S_EX2 = 4.2662;

     // Governor Data (TGOV1)
     gov3.R = 0.05;
     gov3.T1 = 0.05;
     gov3.Vmax = 5.00;
     gov3.Vmin = -5.00;
     gov3.T2 = 2.1;
     gov3.T3 = 7.0;
     gov3.Dt = 0.0;

     //-----------------Load 5 (bus5)----------------------//
     load5.Name = "LOAD5";
     load5.RealPower = 125e6;
     load5.ReactivePower = 50e6;
     load5.BaseVoltage = 230e3;
     load5.Conductance = load5.RealPower / std::pow(load5.BaseVoltage, 2);
     load5.Susceptance = -load5.ReactivePower / std::pow(load5.BaseVoltage, 2);

     //-----------------Load 6 (bus6)----------------------//
     load6.Name = "LOAD6";
     load6.RealPower = 90e6;
     load6.ReactivePower = 30e6;
     load6.BaseVoltage = 230e3;
     load6.Conductance = load6.RealPower / std::pow(load6.BaseVoltage, 2);
     load6.Susceptance = -load6.ReactivePower / std::pow(load6.BaseVoltage, 2);

     //-----------------Load 8 (bus8)----------------------//
     load8.Name = "LOAD8";
     load8.RealPower = 100e6;
     load8.ReactivePower = 35e6;
     load8.BaseVoltage = 230e3;
     load8.Conductance = load8.RealPower / std::pow(load8.BaseVoltage, 2);
     load8.Susceptance = -load8.ReactivePower / std::pow(load8.BaseVoltage, 2);

     //-----------------Line 4-5---------------------------//
     line54.Name = "LINE54";
     line54.ResistancePU = 0.01005;//0.01;
     line54.ReactancePU = 0.08521;//0.085;
     line54.SusceptancePU = 1.0/5.688921;//0.176;
     line54.Conductance = 0.0;
     line54.Resistance = line54.ResistancePU * baseZ;
     line54.Inductance = line54.ReactancePU * baseZ / nomOmega;
     line54.Capacitance = line54.SusceptancePU / baseZ / nomOmega;
     line54.BaseVoltage = 230e3;

     //-----------------Line 4-6---------------------------//
     line64.Name = "LINE64";
     line64.ResistancePU = 0.017083;//0.017;
     line64.ReactancePU = 0.092216;//0.092;
     line64.SusceptancePU = 1.0/6.336801;//0.158;
     line64.Conductance = 0.0;
     line64.Resistance = line64.ResistancePU * baseZ;
     line64.Inductance = line64.ReactancePU * baseZ / nomOmega;
     line64.Capacitance = line64.SusceptancePU / baseZ / nomOmega;
     line64.BaseVoltage = 230e3;

     //-----------------Line 5-7---------------------------//
     line75.Name = "LINE75";
     line75.ResistancePU = 0.032533;//0.032;
     line75.ReactancePU = 0.162281;//0.161;
     line75.SusceptancePU = 1.0/3.28151;//0.306;
     line75.Conductance = 0.0;
     line75.Resistance = line75.ResistancePU * baseZ;
     line75.Inductance = line75.ReactancePU * baseZ / nomOmega;
     line75.Capacitance = line75.SusceptancePU / baseZ / nomOmega;
     line75.BaseVoltage = 230e3;

     //-----------------Line 6-9---------------------------//
     line96.Name = "LINE96";
     line96.ResistancePU = 0.039806;//0.039;
     line96.ReactancePU = 0.171651;//0.17;
     line96.SusceptancePU = 1.0/2.807618;//0.358;
     line96.Conductance = 0.0;
     line96.Resistance = line96.ResistancePU * baseZ;
     line96.Inductance = line96.ReactancePU * baseZ / nomOmega;
     line96.Capacitance = line96.SusceptancePU / baseZ / nomOmega;
     line96.BaseVoltage = 230e3;

     //-----------------Line 7-8---------------------------//
     line78.Name = "LINE78";
     line78.ResistancePU = 0.00853;//0.0085;
     line78.ReactancePU = 0.072127;//0.072;
     line78.SusceptancePU = 1.0/6.717421;//0.149;
     line78.Conductance = 0.0;
     line78.Resistance = line78.ResistancePU * baseZ;
     line78.Inductance = line78.ReactancePU * baseZ / nomOmega;
     line78.Capacitance = line78.SusceptancePU / baseZ / nomOmega;
     line78.BaseVoltage = 230e3;

     //-----------------Line 8-9---------------------------//
     line89.Name = "LINE89";
     line89.ResistancePU = 0.011984;//0.0119;
     line89.ReactancePU = 0.10115;//0.1008;
     line89.SusceptancePU = 1.0/4.793121;//0.209;
     line89.Conductance = 0.0;
     line89.Resistance = line89.ResistancePU * baseZ;
     line89.Inductance = line89.ReactancePU * baseZ / nomOmega;
     line89.Capacitance = line89.SusceptancePU / baseZ / nomOmega;
     line89.BaseVoltage = 230e3;

     //-----------------Step-up Transformer 1-4---------------------------//
     transf14.Name = "TR14";
     transf14.VoltageLVSide = 16.5e3; // Low voltage side
     transf14.VoltageHVSide = 230e3;  // High voltage side
     transf14.RatedPower = 100e6;     // Power rating (100 MVA)
     transf14.Ratio = transf14.VoltageLVSide / transf14.VoltageHVSide;
     transf14.WindingConnection = 0; // Assuming default connection
     transf14.ResistancePU = 0;      // Resistance per unit
     transf14.Resistance = transf14.ResistancePU * baseZ; // Resistance
     transf14.ReactancePU = 0.0576;                       // Reactance per unit
     transf14.Inductance = transf14.ReactancePU * baseZ / nomOmega;
     ; // Inductance

     //-----------------Step-up Transformer 2-7---------------------------//
     transf27.Name = "TR27";
     transf27.VoltageLVSide = 18e3;
     transf27.VoltageHVSide = 230e3;
     transf27.RatedPower = 100e6;
     transf27.Ratio = transf27.VoltageLVSide / transf27.VoltageHVSide;
     transf27.WindingConnection = 0;
     transf27.ResistancePU = 0;
     transf27.Resistance = transf27.ResistancePU * baseZ; // Resistance
     transf27.ReactancePU = 0.0625;
     transf27.Inductance = transf27.ReactancePU * baseZ / nomOmega;
     ;

     //-----------------Step-up Transformer 3-9---------------------------//
     transf39.Name = "TR39";
     transf39.VoltageLVSide = 13.8e3;
     transf39.VoltageHVSide = 230e3;
     transf39.RatedPower = 100e6;
     transf39.Ratio = transf39.VoltageLVSide / transf39.VoltageHVSide;
     transf39.WindingConnection = 0;
     transf39.ResistancePU = 0;
     transf39.Resistance = transf39.ResistancePU * baseZ; // Resistance
     transf39.ReactancePU = 0.0586;
     transf39.Inductance = transf39.ReactancePU * baseZ / nomOmega;
     ;
   }
 };
 } // namespace NineBus
 namespace Components {
namespace SynchronousGeneratorKundur {
// P. Kundur, "Power System Stability and Control", Example 3.2, pp. 102
// and Example 3.5, pp. 134f.
struct MachineParameters {
  // Thermal generating unit, 3600r/min, 2-pole
  Real nomPower = 555e6;
  Real nomVoltage = 24e3; // Phase-to-Phase RMS
  Real nomFreq = 60;
  Real nomFieldCurr = 1300;
  Int poleNum = 2;
  Real H = 3.7;

  // Define machine parameters in per unit
  // Fundamental parameters
  Real Rs = 0.003;
  Real Ll = 0.15;
  Real Lmd = 1.6599;
  Real Lmq = 1.61;
  Real Rfd = 0.0006;
  Real Llfd = 0.1648;
  Real Rkd = 0.0284;
  Real Llkd = 0.1713;
  Real Rkq1 = 0.0062;
  Real Llkq1 = 0.7252;
  Real Rkq2 = 0.0237;
  Real Llkq2 = 0.125;

  // Operational parameters
  Real Td0_t = 8.0669;
  Real Td0_s = 0.0300;
  Real Td_t = 1.3368;
  Real Td_s = 0.0230;
  Real Ld_t = 0.2999;
  Real Ld_s = 0.2299;
  Real Tq0_t = 0.9991;
  Real Tq0_s = 0.0700;
  Real Lq_t = 0.6500;
  Real Lq_s = 0.2500;
  Real Ld = 1.8099;
  Real Lq = 1.7600;
};
} // namespace SynchronousGeneratorKundur

namespace GovernorKundur {
struct Parameters {
  // Turbine model parameters (tandem compound single reheat steam turbine, fossil-fuelled)
  // from P. Kundur, "Power System Stability and Control", 1994, p. 427
  Real Ta_t = 0.3; // T_CH
  Real Fa = 0.3;   // F_HP
  Real Tb = 7.0;   // T_RH
  Real Fb = 0.3;   // F_IP
  Real Tc = 0.5;   // T_CO
  Real Fc = 0.4;   // F_LP

  // Governor parameters (mechanical-hydraulic control)
  // from P. Kundur, "Power System Stability and Control", 1994, p. 437
  Real Kg = 20; // 5% droop
  Real Tsr = 0.1;
  Real Tsm = 0.3;
};
} // namespace GovernorKundur

namespace ExcitationSystemEremia {
struct Parameters {
  // Excitation system parameters (IEEE Type DC1A)
  // from M. Eremia, "Handbook of Electrical Power System Dynamics", 2013, p.96 and 106
  // voltage-regulator
  Real Ka = 46;
  Real Ta = 0.06;
  // exciter
  Real Ke = -0.0435;
  Real Te = 0.46;
  // stabilizing feedback
  Real Kf = 0.1;
  Real Tf = 1;
  // voltage transducer
  Real Tr = 0.02;
};
} // namespace ExcitationSystemEremia
} // namespace Components
} // namespace Examples
} // namespace CIM
} // namespace CPS
