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
   Real LoadStep;
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
   Real LoadStep = 90e6 * 0.33; // Load 6 = 90e6;

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
   Line line45;
   Line line46;
   Line line57;
   Line line69;
   Line line78;
   Line line89;

   //-----------Transformers-----------//
   Transformer transf14;
   Transformer transf27;
   Transformer transf39;

   ScenarioConfig() {
     //-----------------Generator 1 (bus1)-----------------//
     gen1.Name = "SynGen1";
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
     gen2.Name = "SynGen2";
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

     // // Exciter Data-1 (IEEE Type DC1A)
     // gen2.KA = 20.0;
     // gen2.TA = 0.2;
     // gen2.VRmin = -5.0;
     // gen2.VRmax = 5.0;
     // gen2.KE = 1.0;
     // gen2.TE = 0.314;
     // gen2.KF = 0.063;
     // gen2.TF = 0.35;

     // // Exciter Data-2 (IEEE Type DC1A)
     // gen2.EX1 = 3.3;
     // gen2.S_EX1 = 0.6602;
     // gen2.EX2 = 4.5;
     // gen2.S_EX2 = 4.2662;

     // // Governor Data (TGOV1)
     // gen2.R = 0.05;
     // gen2.T1 = 0.05;
     // gen2.Vmax = 5.00;
     // gen2.Vmin = -5.00;
     // gen2.T2 = 2.1;
     // gen2.T3 = 7.0;
     // gen2.Dt = 0.0;

     //-----------------Generator 3 (bus3)-----------------//
     gen3.Name = "SynGen3";
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

     // // Exciter Data-1 (IEEE Type DC1A)
     // gen3.KA = 20.0;
     // gen3.TA = 0.2;
     // gen3.VRmin = -5.0;
     // gen3.VRmax = 5.0;
     // gen3.KE = 1.0;
     // gen3.TE = 0.314;
     // gen3.KF = 0.063;
     // gen3.TF = 0.35;

     // // Exciter Data-2 (IEEE Type DC1A)
     // gen3.EX1 = 3.3;
     // gen3.S_EX1 = 0.6602;
     // gen3.EX2 = 4.5;
     // gen3.S_EX2 = 4.2662;

     // // Governor Data (TGOV1)
     // gen3.R = 0.05;
     // gen3.T1 = 0.05;
     // gen3.Vmax = 5.00;
     // gen3.Vmin = -5.00;
     // gen3.T2 = 2.1;
     // gen3.T3 = 7.0;
     // gen3.Dt = 0.0;

     //-----------------Load 5 (bus5)----------------------//
     load5.Name = "Load5";
     load5.RealPower = 125e6;
     load5.ReactivePower = 50e6;
     load5.BaseVoltage = 230e3;
     load5.Conductance = load5.RealPower / std::pow(load5.BaseVoltage, 2);
     load5.Susceptance = -load5.ReactivePower / std::pow(load5.BaseVoltage, 2);

     //-----------------Load 6 (bus6)----------------------//
     load6.Name = "Load6";
     load6.RealPower = 90e6;
     load6.ReactivePower = 30e6;
     load6.BaseVoltage = 230e3;
     load6.Conductance = load6.RealPower / std::pow(load6.BaseVoltage, 2);
     load6.Susceptance = -load6.ReactivePower / std::pow(load6.BaseVoltage, 2);


     //-----------------Load 8 (bus8)----------------------//
     load8.Name = "Load8";
     load8.RealPower = 100e6;
     load8.ReactivePower = 35e6;
     load8.BaseVoltage = 230e3;
     load8.Conductance = load8.RealPower / std::pow(load8.BaseVoltage, 2);
     load8.Susceptance = -load8.ReactivePower / std::pow(load8.BaseVoltage, 2);

     //-----------------Line 4-5---------------------------//
     line45.Name = "Line45";
     line45.ResistancePU = 0.01005;//0.01;
     line45.ReactancePU = 0.08521;//0.085;
     line45.SusceptancePU = 1.0/5.688921;//0.176;
     line45.Conductance = 0.0;
     line45.Resistance = line45.ResistancePU * baseZ;
     line45.Inductance = line45.ReactancePU * baseZ / nomOmega;
     line45.Capacitance = line45.SusceptancePU / baseZ / nomOmega;
     line45.BaseVoltage = 230e3;

     //-----------------Line 4-6---------------------------//
     line46.Name = "Line46";
     line46.ResistancePU = 0.017083;//0.017;
     line46.ReactancePU = 0.092216;//0.092;
     line46.SusceptancePU = 1.0/6.336801;//0.158;
     line46.Conductance = 0.0;
     line46.Resistance = line46.ResistancePU * baseZ;
     line46.Inductance = line46.ReactancePU * baseZ / nomOmega;
     line46.Capacitance = line46.SusceptancePU / baseZ / nomOmega;
     line46.BaseVoltage = 230e3;

     //-----------------Line 5-7---------------------------//
     line57.Name = "Line57";
     line57.ResistancePU = 0.032533;//0.032;
     line57.ReactancePU = 0.162281;//0.161;
     line57.SusceptancePU = 1.0/3.28151;//0.306;
     line57.Conductance = 0.0;
     line57.Resistance = line57.ResistancePU * baseZ;
     line57.Inductance = line57.ReactancePU * baseZ / nomOmega;
     line57.Capacitance = line57.SusceptancePU / baseZ / nomOmega;
     line57.BaseVoltage = 230e3;

     //-----------------Line 6-9---------------------------//
     line69.Name = "Line69";
     line69.ResistancePU = 0.039806;//0.039;
     line69.ReactancePU = 0.171651;//0.17;
     line69.SusceptancePU = 1.0/2.807618;//0.358;
     line69.Conductance = 0.0;
     line69.Resistance = line69.ResistancePU * baseZ;
     line69.Inductance = line69.ReactancePU * baseZ / nomOmega;
     line69.Capacitance = line69.SusceptancePU / baseZ / nomOmega;
     line69.BaseVoltage = 230e3;

     //-----------------Line 7-8---------------------------//
     line78.Name = "Line78";
     line78.ResistancePU = 0.00853;//0.0085;
     line78.ReactancePU = 0.072127;//0.072;
     line78.SusceptancePU = 1.0/6.717421;//0.149;
     line78.Conductance = 0.0;
     line78.Resistance = line78.ResistancePU * baseZ;
     line78.Inductance = line78.ReactancePU * baseZ / nomOmega;
     line78.Capacitance = line78.SusceptancePU / baseZ / nomOmega;
     line78.BaseVoltage = 230e3;

     //-----------------Line 8-9---------------------------//
     line89.Name = "Line89";
     line89.ResistancePU = 0.011984;//0.0119;
     line89.ReactancePU = 0.10115;//0.1008;
     line89.SusceptancePU = 1.0/4.793121;//0.209;
     line89.Conductance = 0.0;
     line89.Resistance = line89.ResistancePU * baseZ;
     line89.Inductance = line89.ReactancePU * baseZ / nomOmega;
     line89.Capacitance = line89.SusceptancePU / baseZ / nomOmega;
     line89.BaseVoltage = 230e3;

     //-----------------Step-up Transformer 1-4---------------------------//
     transf14.Name = "Transformer14";
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
     transf27.Name = "Transformer27";
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
     transf39.Name = "Transformer39";
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

  struct Cosim_9bus {
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
  
  Cosim_9bus() {
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
     load6.LoadStep = 90e6 * 0.33; // Load 6 = 90e6;

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
   };
 }; //struct Cosim_9bus
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

namespace TurbineGovernor {
struct TurbineGovernorPSAT1 {
  // Turbine Governor type 1
  // Taken from from PSAT - example d_014_pss_l14

  // Reference speed (p.u.)
  Real OmegaRef = 1.0;
  // Pilot valve droop (p.u.)
  Real R = 0.02;
  // Maximum Torque (p.u.)
  Real Tmax = 1.2;
  // Minimim Torque (p.u.)
  Real Tmin = 0.3;
  // Governor time constant (s)
  Real Ts = 0.1;
  // Servo time constant (s)
  Real Tc = 0.45;
  // Transient gain time constant (s)
  Real T3 = 0.0;
  // Power fraction time constant (s)
  Real T4 = 12.0;
  // Reheat time constant (s)
  Real T5 = 50.0;
};

struct TurbineGovernorPSAT2 {
  // Turbine Governor type 1
  // Taken from PSAT - example d_anderson_farmer

  // Reference speed (p.u.)
  Real OmegaRef = 1.0;
  // Pilot valve droop (p.u.)
  Real R = 0.04;
  // Maximum Torque (p.u.)
  Real Tmax = 100;
  // Minimim Torque (p.u.)
  Real Tmin = 0.0;
  // Governor time constant (s)
  Real Ts = 20;
  // Servo time constant (s)
  Real Tc = 0.2;
  // Transient gain time constant (s)
  Real T3 = 0.2;
  // Power fraction time constant (s)
  Real T4 = 0.2;
  // Reheat time constant (s)
  Real T5 = 0.2;
  };
} // namespace TurbineGovernor
} // namespace Components

namespace Grids {
namespace CIGREHVEuropean {
struct LineParameters {
  // 220 kV
  Real Vnom = 220e3;
  Real nomFreq = 50;
  Real nomOmega = nomFreq * 2 * PI;
  //PiLine parameters (given in [p.u/km] with Z_base= 529ohm)
  Real lineResistancePerKm = 1.35e-4 * 529;
  Real lineReactancePerKm = 8.22e-4 * 529;
  Real lineSusceptancePerKm = 1.38e-3 / 529;
  Real lineConductancePerKm = 0;
};
} // namespace CIGREHVEuropean
namespace CIGREHVAmerican {
struct LineParameters {
  // 230 kV
  Real Vnom = 230e3;
  Real nomFreq = 60;
  Real nomOmega = nomFreq * 2 * PI;
  //PiLine parameters (given in [p.u/km] with Z_base= 529ohm)
  Real lineResistancePerKm = 1.27e-4 * 529;
  Real lineReactancePerKm = 9.05e-4 * 529;
  Real lineSusceptancePerKm = 1.81e-3 / 529;
  Real lineConductancePerKm = 0;
};
} // namespace CIGREHVAmerican

namespace KundurExample1 {
// P. Kundur, "Power System Stability and Control", Example 13.2, pp. 864-869.
struct Network {
  Real nomVoltage = 400e3;
};

struct Gen {
  Real nomPower = 2220e6;
  Real nomVoltage = 400e3;
  Real H = 3.5;
  Real XpdPU = 0.3;
  Real RsPU = 0;
  Real D = 1.0;
};
struct Line1 {
  // Vnom = 400kV
  Real lineResistance = 0.0721;
  Real lineReactance = 36.0360;
  Real lineSusceptance = 0;
  Real lineConductance = 0;
};

struct Transf1 {
  Real nomVoltageHV = 400e3;
  Real nomVoltageMV = 400e3;
  Real transformerResistance = 0;      // referred to HV side
  Real transformerReactance = 10.8108; // referred to HV side
};
} // namespace KundurExample1

namespace SMIB {
struct ScenarioConfig {
  //-----------Network-----------//
  Real Vnom = 230e3;
  Real nomFreq = 60;
  Real nomOmega = nomFreq * 2 * PI;
  //-----------Generator-----------//
  Real nomPower = 500e6;
  Real nomPhPhVoltRMS = 22e3;
  Real H = 5;
  Real Xpd = 0.31;
  Real Rs = 0.003 * 0;
  Real D = 1.5;
  // Initialization parameters
  Real initMechPower = 300e6;
  Real initActivePower = 300e6;
  Real setPointVoltage = nomPhPhVoltRMS + 0.05 * nomPhPhVoltRMS;
  //-----------Transformer-----------//
  Real t_ratio = Vnom / nomPhPhVoltRMS;
  //-----------Transmission Line-----------//
  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  Real lineLeng = 100;
  Real lineResistance = lineCIGREHV.lineResistancePerKm * lineLeng;
  Real lineInductance = lineCIGREHV.lineReactancePerKm / nomOmega * lineLeng;
  Real lineCapacitance = lineCIGREHV.lineSusceptancePerKm / nomOmega * lineLeng;
  Real lineConductance = lineCIGREHV.lineConductancePerKm * lineLeng;
};

struct ScenarioConfig2 {
  //Scenario used to validate reduced order SG VBR models against PSAT (in SP domain)

  // General grid parameters
  Real VnomMV = 24e3;
  Real VnomHV = 230e3;
  Real nomFreq = 60;
  Real ratio = VnomMV / VnomHV;
  Real nomOmega = nomFreq * 2 * PI;

  // Generator parameters
  Real setPointActivePower = 300e6;
  Real setPointVoltage = 1.05 * VnomMV;

  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  Real lineLength = 100;
  Real lineResistance =
      lineCIGREHV.lineResistancePerKm * lineLength * std::pow(ratio, 2);
  Real lineInductance = lineCIGREHV.lineReactancePerKm * lineLength *
                        std::pow(ratio, 2) / nomOmega;
  Real lineCapacitance = lineCIGREHV.lineSusceptancePerKm * lineLength /
                         std::pow(ratio, 2) / nomOmega;
  Real lineConductance = 1e-15;

  // In PSAT SwitchClosed is equal to 1e-3 p.u.
  Real SwitchClosed = 0.001 * (24 * 24 / 555);
  Real SwitchOpen = 1e6;
};

struct ScenarioConfig3 {
  //Scenario used to validate reduced order SG VBR models in the DP and EMT domain against the SP domain

  // General grid parameters
  Real VnomMV = 24e3;
  Real nomFreq = 60;
  Real nomOmega = nomFreq * 2 * PI;

  //-----------Generator-----------//
  Real setPointActivePower = 300e6;
  Real mechPower = 300e6;
  Real initActivePower = 300e6;
  Real initReactivePower = 0;
  Real initVoltAngle = -PI / 2;
  Complex initComplexElectricalPower =
      Complex(initActivePower, initReactivePower);
  Complex initTerminalVolt =
      VnomMV * Complex(cos(initVoltAngle), sin(initVoltAngle));

  //
  Real SwitchClosed = 0.1;
  Real SwitchOpen = 1e6;
};

namespace ReducedOrderSynchronGenerator {
namespace Scenario4 {
//Scenario used to compare DP against SP accuracy in Martin's thesis

struct Config {
  // default configuration of scenario
  // adjustable using applyCommandLineArgsOptions
  String sgType = "4";
  Real startTimeFault = 30.0;
  Real endTimeFault = 30.1;
};

struct GridParams {
  // General grid parameters
  Real VnomMV = 24e3;
  Real VnomHV = 230e3;
  Real nomFreq = 60;
  Real ratio = VnomMV / VnomHV;
  Real nomOmega = nomFreq * 2 * PI;

  // Generator parameters
  Real setPointActivePower = 300e6;
  Real setPointVoltage = 1.05 * VnomMV;

  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  Real lineLength = 100;
  Real lineResistance =
      lineCIGREHV.lineResistancePerKm * lineLength * std::pow(ratio, 2);
  Real lineInductance = lineCIGREHV.lineReactancePerKm * lineLength *
                        std::pow(ratio, 2) / nomOmega;
  Real lineCapacitance = lineCIGREHV.lineSusceptancePerKm * lineLength /
                         std::pow(ratio, 2) / nomOmega;
  Real lineConductance = 8e-2;

  // In PSAT SwitchClosed is equal to 1e-3 p.u.
  Real SwitchClosed = 0.1;
  Real SwitchOpen = 1e6;
};
} // namespace Scenario4

namespace Scenario5 {
// SMIB scenario with RX trafo and load step as event
struct Config {
  // default configuration of scenario
  // adjustable using applyCommandLineArgsOptions
  String sgType = "4";
  Real startTimeFault = 1.0;
  Real endTimeFault = 1.1;
};

struct GridParams {

  // General grid parameters
  Real VnomMV = 24e3;
  Real VnomHV = 230e3;
  Real nomFreq = 60;
  Real ratio = VnomMV / VnomHV;
  Real nomOmega = nomFreq * 2 * PI;

  // Generator parameters
  Real setPointActivePower = 300e6;
  Real setPointVoltage = 1.05 * VnomMV;

  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  Real lineLength = 100;
  Real lineResistance = lineCIGREHV.lineResistancePerKm * lineLength;
  Real lineInductance = lineCIGREHV.lineReactancePerKm * lineLength / nomOmega;
  Real lineCapacitance =
      lineCIGREHV.lineSusceptancePerKm * lineLength / nomOmega;
  Real lineConductance = 1.0491e-05; // Psnub 0.1% of 555MW

  // Switch for load step
  Real SwitchClosed = 529;      // 100 MW load step
  Real SwitchOpen = 9.1840e+07; // corresponds to 1e6 Ohms at MV level of 24kV
};

struct Transf1 {
  Real nomVoltageHV = 230e3;
  Real nomVoltageMV = 24e3;
  Real transformerResistance = 0;     // referred to HV side
  Real transformerReactance = 5.2900; // referred to HV side
  Real transformerNominalPower = 555e6;
};
} // namespace Scenario5

namespace Scenario6 {
// SMIB scenario with ideal trafo and load step as event

struct Config {
  // default configuration of scenario
  // adjustable using applyCommandLineArgsOptions
  String sgType = "4";
  Real loadStepEventTime = 10.0;

  // parameters of iterative model
  Real maxIter = 25;
  Real tolerance = 1e-10;
};

struct GridParams {
  // General grid parameters
  Real VnomMV = 24e3;
  Real VnomHV = 230e3;
  Real nomFreq = 60;
  Real ratio = VnomMV / VnomHV;
  Real nomOmega = nomFreq * 2 * PI;

  // Generator parameters
  Real setPointActivePower = 300e6;
  Real setPointVoltage = 1.05 * VnomMV;

  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  Real lineLength = 100;
  Real lineResistance =
      lineCIGREHV.lineResistancePerKm * lineLength * std::pow(ratio, 2);
  Real lineInductance = lineCIGREHV.lineReactancePerKm * lineLength *
                        std::pow(ratio, 2) / nomOmega;
  Real lineCapacitance = lineCIGREHV.lineSusceptancePerKm * lineLength /
                         std::pow(ratio, 2) / nomOmega;
  Real lineConductance = 0.0048; // Psnub 0.5% of 555MW

  // Load step
  Real loadStepActivePower = 100e6;
};
} // namespace Scenario6

} // namespace ReducedOrderSynchronGenerator
} // namespace SMIB

namespace ThreeBus {
struct ScenarioConfig {

  //-----------Network-----------//
  Real Vnom = 230e3;
  Real nomFreq = 60;
  Real nomOmega = nomFreq * 2 * PI;

  //-----------Generator 1 (bus1)-----------//
  Real nomPower_G1 = 300e6;
  Real nomPhPhVoltRMS_G1 = 25e3;
  Real nomFreq_G1 = 60;
  Real H_G1 = 6;
  Real Xpd_G1 = 0.3;      //in p.u
  Real Rs_G1 = 0.003 * 0; //in p.u
  Real D_G1 = 1.5;        //in p.u
  // Initialization parameters
  Real initActivePower_G1 = 270e6;
  Real initMechPower_G1 = 270e6;
  Real setPointVoltage_G1 = nomPhPhVoltRMS_G1 + 0.05 * nomPhPhVoltRMS_G1;

  //-----------Generator 2 (bus2)-----------//
  Real nomPower_G2 = 50e6;
  Real nomPhPhVoltRMS_G2 = 13.8e3;
  Real nomFreq_G2 = 60;
  Real H_G2 = 2;
  Real Xpd_G2 = 0.1;      //in p.u
  Real Rs_G2 = 0.003 * 0; //in p.u
  Real D_G2 = 1.5;        //in p.u
  // Initialization parameters
  Real initActivePower_G2 = 45e6;
  Real initMechPower_G2 = 45e6;
  Real setPointVoltage_G2 = nomPhPhVoltRMS_G2 - 0.05 * nomPhPhVoltRMS_G2;

  //-----------Transformers-----------//
  Real t1_ratio = Vnom / nomPhPhVoltRMS_G1;
  Real t2_ratio = Vnom / nomPhPhVoltRMS_G2;

  //-----------Load (bus3)-----------
  Real activePower_L = 310e6;
  Real reactivePower_L = 150e6;

  // -----------Transmission Lines-----------//
  // CIGREHVAmerican (230 kV)
  Grids::CIGREHVAmerican::LineParameters lineCIGREHV;
  //line 1-2 (180km)
  Real lineResistance12 = lineCIGREHV.lineResistancePerKm * 180;
  Real lineInductance12 = lineCIGREHV.lineReactancePerKm / nomOmega * 180;
  Real lineCapacitance12 = lineCIGREHV.lineSusceptancePerKm / nomOmega * 180;
  Real lineConductance12 = lineCIGREHV.lineConductancePerKm * 180;
  //line 1-3 (150km)
  Real lineResistance13 = lineCIGREHV.lineResistancePerKm * 150;
  Real lineInductance13 = lineCIGREHV.lineReactancePerKm / nomOmega * 150;
  Real lineCapacitance13 = lineCIGREHV.lineSusceptancePerKm / nomOmega * 150;
  Real lineConductance13 = lineCIGREHV.lineConductancePerKm * 150;
  //line 2-3 (80km)
  Real lineResistance23 = lineCIGREHV.lineResistancePerKm * 80;
  Real lineInductance23 = lineCIGREHV.lineReactancePerKm / nomOmega * 80;
  Real lineCapacitance23 = lineCIGREHV.lineSusceptancePerKm / nomOmega * 80;
  Real lineConductance23 = lineCIGREHV.lineConductancePerKm * 80;
};
} // namespace ThreeBus

namespace SGIB {

struct ScenarioConfig {
  Real systemFrequency = 50;
  Real systemNominalVoltage = 20e3;

  // Line parameters (R/X = 1)
  Real length = 5;
  Real lineResistance = 0.5 * length;
  Real lineInductance = 0.5 / 314 * length;
  Real lineCapacitance = 50e-6 / 314 * length;

  // PV controller parameters
  Real scalingKp = 1;
  Real scalingKi = 0.1;

  Real KpPLL = 0.25 * scalingKp;
  Real KiPLL = 2 * scalingKi;
  Real KpPowerCtrl = 0.001 * scalingKp;
  Real KiPowerCtrl = 0.08 * scalingKi;
  Real KpCurrCtrl = 0.3 * scalingKp;
  Real KiCurrCtrl = 10 * scalingKi;
  Real OmegaCutoff = 2 * PI * systemFrequency;

  // Initial state values
  Real thetaPLLInit = 0; // only for debug
  Real phiPLLInit = 0;   // only for debug
  Real phi_dInit = 0;
  Real phi_qInit = 0;
  Real gamma_dInit = 0;
  Real gamma_qInit = 0;

  // Nominal generated power values of PV
  Real pvNominalVoltage = 1500.;
  Real pvNominalActivePower = 100e3;
  Real pvNominalReactivePower = 50e3;

  // PV filter parameters
  Real Lf = 0.002;
  Real Cf = 789.3e-6;
  Real Rf = 0.1;
  Real Rc = 0.1;

  // PV connection transformer parameters
  Real transformerNominalPower = 5e6;
  Real transformerInductance = 0.928e-3;

  // Further parameters
  Real systemOmega = 2 * PI * systemFrequency;
};
} // namespace SGIB

namespace CIGREMV {

struct ScenarioConfig {
  Real systemFrequency = 50;
  Real systemNominalVoltage = 20e3;
  Real penetrationLevel = 1;
  Real totalLoad =
      4319.1e3; // calculated total load in CIGRE MV left feeder (see CIM data)

  // parameters of one PV unit
  Real pvUnitNominalVoltage = 1500.;
  Real pvUnitNominalPower = 50e3;
  Real pvUnitPowerFactor = 1;

  // calculate PV units per plant to reach penetration level
  Int numberPVUnits = Int(totalLoad * penetrationLevel / pvUnitNominalPower);
  Int numberPVPlants = 9;
  Int numberPVUnitsPerPlant = numberPVUnits / numberPVPlants;

  // PV controller parameters
  Real scalingKp = 1;
  Real scalingKi = 0.001;

  Real KpPLL = 0.25 * scalingKp;
  Real KiPLL = 2 * scalingKi;
  Real KpPowerCtrl = 0.001 * scalingKp;
  Real KiPowerCtrl = 0.08 * scalingKi;
  Real KpCurrCtrl = 0.3 * scalingKp;
  Real KiCurrCtrl = 10 * scalingKi;
  Real OmegaCutoff = 2 * PI * systemFrequency;

  // PV filter parameters
  Real Lf = 0.002;
  Real Cf = 789.3e-6;
  Real Rf = 0.1;
  Real Rc = 0.1;

  // PV connection transformer parameters
  Real transformerNominalPower = 5e6;
  Real transformerInductance = 0.928e-3;

  // Further parameters
  Real systemOmega = 2 * PI * systemFrequency;

  // Initial state values (use known values with scaled control params)
  Real thetaPLLInit = 314.168313 - systemOmega;
  Real phiPLLInit = 8e-06;
  Real pInit = 450000.716605;
  Real qInit = -0.577218;
  Real phi_dInit = 3854197405 * scalingKi;
  Real phi_qInit = -3737 * scalingKi;
  Real gamma_dInit = 128892668 * scalingKi;
  Real gamma_qInit = 23068682 * scalingKi;
};

void addInvertersToCIGREMV(SystemTopology &system,
                           CIGREMV::ScenarioConfig scenario, Domain domain) {
  Real pvActivePower =
      scenario.pvUnitNominalPower * scenario.numberPVUnitsPerPlant;
  Real pvReactivePower =
      sqrt(std::pow(pvActivePower / scenario.pvUnitPowerFactor, 2) -
           std::pow(pvActivePower, 2));

  // add PVs to network topology
  for (Int n = 3; n <= 11; ++n) {
    // TODO: cast to BaseAverageVoltageSourceInverter and move set functions out of case distinction
    if (domain == Domain::SP) {
      SimNode<Complex>::Ptr connectionNode =
          system.node<CPS::SimNode<Complex>>("N" + std::to_string(n));
      auto pv = SP::Ph1::AvVoltageSourceInverterDQ::make(
          "pv_" + connectionNode->name(), "pv_" + connectionNode->name(),
          Logger::Level::debug, true);
      pv->setParameters(scenario.systemOmega, scenario.pvUnitNominalVoltage,
                        pvActivePower, pvReactivePower);
      pv->setControllerParameters(scenario.KpPLL, scenario.KiPLL,
                                  scenario.KpPowerCtrl, scenario.KiPowerCtrl,
                                  scenario.KpCurrCtrl, scenario.KiCurrCtrl,
                                  scenario.OmegaCutoff);
      pv->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                              scenario.Rc);
      pv->setTransformerParameters(
          scenario.systemNominalVoltage, scenario.pvUnitNominalVoltage,
          scenario.transformerNominalPower,
          scenario.systemNominalVoltage / scenario.pvUnitNominalVoltage, 0, 0,
          scenario.transformerInductance);
      pv->setInitialStateValues(scenario.pInit, scenario.qInit,
                                scenario.phi_dInit, scenario.phi_qInit,
                                scenario.gamma_dInit, scenario.gamma_qInit);
      system.addComponent(pv);
      system.connectComponentToNodes<Complex>(pv, {connectionNode});
    } else if (domain == Domain::DP) {
      SimNode<Complex>::Ptr connectionNode =
          system.node<CPS::SimNode<Complex>>("N" + std::to_string(n));
      auto pv = DP::Ph1::AvVoltageSourceInverterDQ::make(
          "pv_" + connectionNode->name(), "pv_" + connectionNode->name(),
          Logger::Level::debug, true);
      pv->setParameters(scenario.systemOmega, scenario.pvUnitNominalVoltage,
                        pvActivePower, pvReactivePower);
      pv->setControllerParameters(scenario.KpPLL, scenario.KiPLL,
                                  scenario.KpPowerCtrl, scenario.KiPowerCtrl,
                                  scenario.KpCurrCtrl, scenario.KiCurrCtrl,
                                  scenario.OmegaCutoff);
      pv->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                              scenario.Rc);
      pv->setTransformerParameters(
          scenario.systemNominalVoltage, scenario.pvUnitNominalVoltage,
          scenario.transformerNominalPower,
          scenario.systemNominalVoltage / scenario.pvUnitNominalVoltage, 0, 0,
          scenario.transformerInductance);
      pv->setInitialStateValues(scenario.pInit, scenario.qInit,
                                scenario.phi_dInit, scenario.phi_qInit,
                                scenario.gamma_dInit, scenario.gamma_qInit);
      system.addComponent(pv);
      system.connectComponentToNodes<Complex>(pv, {connectionNode});
    } else if (domain == Domain::EMT) {
      SimNode<Real>::Ptr connectionNode =
          system.node<CPS::SimNode<Real>>("N" + std::to_string(n));
      auto pv = EMT::Ph3::AvVoltageSourceInverterDQ::make(
          "pv_" + connectionNode->name(), "pv_" + connectionNode->name(),
          Logger::Level::debug, true);
      pv->setParameters(scenario.systemOmega, scenario.pvUnitNominalVoltage,
                        pvActivePower, pvReactivePower);
      pv->setControllerParameters(scenario.KpPLL, scenario.KiPLL,
                                  scenario.KpPowerCtrl, scenario.KiPowerCtrl,
                                  scenario.KpCurrCtrl, scenario.KiCurrCtrl,
                                  scenario.OmegaCutoff);
      pv->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                              scenario.Rc);
      pv->setTransformerParameters(
          scenario.systemNominalVoltage, scenario.pvUnitNominalVoltage,
          scenario.transformerNominalPower,
          scenario.systemNominalVoltage / scenario.pvUnitNominalVoltage, 0, 0,
          scenario.transformerInductance, scenario.systemOmega);
      pv->setInitialStateValues(scenario.pInit, scenario.qInit,
                                scenario.phi_dInit, scenario.phi_qInit,
                                scenario.gamma_dInit, scenario.gamma_qInit);
      system.addComponent(pv);
      system.connectComponentToNodes<Real>(pv, {connectionNode});
    }
  }
}

void logPVAttributes(DPsim::DataLogger::Ptr logger,
                     CPS::TopologicalPowerComp::Ptr pv) {

  // power controller
  std::vector<String> inputNames = {pv->name() + "_powerctrl_input_pref",
                                    pv->name() + "_powerctrl_input_qref",
                                    pv->name() + "_powerctrl_input_vcd",
                                    pv->name() + "_powerctrl_input_vcq",
                                    pv->name() + "_powerctrl_input_ircd",
                                    pv->name() + "_powerctrl_input_ircq"};
  logger->logAttribute(inputNames, pv->attribute("powerctrl_inputs"));
  std::vector<String> stateNames = {pv->name() + "_powerctrl_state_p",
                                    pv->name() + "_powerctrl_state_q",
                                    pv->name() + "_powerctrl_state_phid",
                                    pv->name() + "_powerctrl_state_phiq",
                                    pv->name() + "_powerctrl_state_gammad",
                                    pv->name() + "_powerctrl_state_gammaq"};
  logger->logAttribute(stateNames, pv->attribute("powerctrl_states"));
  std::vector<String> outputNames = {pv->name() + "_powerctrl_output_vsd",
                                     pv->name() + "_powerctrl_output_vsq"};
  logger->logAttribute(outputNames, pv->attribute("powerctrl_outputs"));

  // interface variables
  logger->logAttribute(pv->name() + "_v_intf", pv->attribute("v_intf"));
  logger->logAttribute(pv->name() + "_i_intf", pv->attribute("i_intf"));

  // additional variables
  logger->logAttribute(pv->name() + "_pll_output", pv->attribute("pll_output"));
  logger->logAttribute(pv->name() + "_vsref", pv->attribute("Vsref"));
  logger->logAttribute(pv->name() + "_vs", pv->attribute("Vs"));
}

} // namespace CIGREMV
} // namespace Grids

namespace Events {
std::shared_ptr<DPsim::SwitchEvent> createEventAddPowerConsumption(
    String nodeName, Real eventTime, Real additionalActivePower,
    SystemTopology &system, Domain domain, DPsim::DataLogger::Ptr logger) {

  // TODO: use base classes ph1
  if (domain == CPS::Domain::DP) {
    auto loadSwitch = DP::Ph1::Switch::make("Load_Add_Switch_" + nodeName,
                                            Logger::Level::debug);
    auto connectionNode = system.node<CPS::SimNode<Complex>>(nodeName);
    Real resistance = std::abs(connectionNode->initialSingleVoltage()) *
                      std::abs(connectionNode->initialSingleVoltage()) /
                      additionalActivePower;
    loadSwitch->setParameters(1e9, resistance);
    loadSwitch->open();
    system.addComponent(loadSwitch);
    system.connectComponentToNodes<Complex>(
        loadSwitch, {CPS::SimNode<Complex>::GND, connectionNode});
    logger->logAttribute("switchedload_i", loadSwitch->attribute("i_intf"));
    return DPsim::SwitchEvent::make(eventTime, loadSwitch, true);
  } else if (domain == CPS::Domain::SP) {
    auto loadSwitch = SP::Ph1::Switch::make("Load_Add_Switch_" + nodeName,
                                            Logger::Level::debug);
    auto connectionNode = system.node<CPS::SimNode<Complex>>(nodeName);
    Real resistance = std::abs(connectionNode->initialSingleVoltage()) *
                      std::abs(connectionNode->initialSingleVoltage()) /
                      additionalActivePower;
    loadSwitch->setParameters(1e9, resistance);
    loadSwitch->open();
    system.addComponent(loadSwitch);
    system.connectComponentToNodes<Complex>(
        loadSwitch, {CPS::SimNode<Complex>::GND, connectionNode});
    logger->logAttribute("switchedload_i", loadSwitch->attribute("i_intf"));
    return DPsim::SwitchEvent::make(eventTime, loadSwitch, true);
  } else {
    return nullptr;
  }
}

std::shared_ptr<DPsim::SwitchEvent3Ph> createEventAddPowerConsumption3Ph(
    String nodeName, Real eventTime, Real additionalActivePower,
    SystemTopology &system, Domain domain, DPsim::DataLogger::Ptr logger) {

  // TODO: use base classes ph3
  if (domain == CPS::Domain::EMT) {
    auto loadSwitch = EMT::Ph3::Switch::make("Load_Add_Switch_" + nodeName,
                                             Logger::Level::debug);
    auto connectionNode = system.node<CPS::SimNode<Real>>(nodeName);
    Real resistance = std::abs(connectionNode->initialSingleVoltage()) *
                      std::abs(connectionNode->initialSingleVoltage()) /
                      additionalActivePower;
    loadSwitch->setParameters(Matrix::Identity(3, 3) * 1e9,
                              Matrix::Identity(3, 3) * resistance);
    loadSwitch->openSwitch();
    system.addComponent(loadSwitch);
    system.connectComponentToNodes<Real>(
        loadSwitch,
        {CPS::SimNode<Real>::GND, system.node<CPS::SimNode<Real>>(nodeName)});
    logger->logAttribute("switchedload_i", loadSwitch->attribute("i_intf"));
    return DPsim::SwitchEvent3Ph::make(eventTime, loadSwitch, true);
  } else {
    return nullptr;
  }
}
} // namespace Events
} // namespace Examples
} // namespace CIM
} // namespace CPS
