%% Parameters for Three-Phase Grid-Tied Inverter Example
%
% This example shows how to control the voltage in a grid-tied inverter
% system. The Voltage regulator subsystem implements the PI-based control
% strategy. The three-phase inverter is connected to the grid via a
% Circuit Breaker. The Circuit Breaker is open at the beginning of the
% simulation to allow synchronization. At time 0.15 seconds, the Circuit
% breaker closes, and the inverter is connected to the grid. The Scopes
% subsystem contains scopes that allow you to see the simulation results.
% The inverter is implemented using IGBTs. To speed up simulation, or for
% real-time deployment, the IGBTs can be replaced with Averaged Switches.
% In this way the gate signals can be averaged over a specified period or
% replaced with modulation waveforms.

% Copyright 2020-2023 The MathWorks, Inc.

%% Control Parameters
Ts   = 1e-6;    % Fundamental sample time       [s]
fsw  = 6000;    % Inverter switching frequency [Hz]
Tsc  = 5e-4;    % Control sample time           [s]

Kp_vd = 0.36;   % Proportional term d-axis voltage controller
Ki_vd = 850;    % Integral term d-axis voltage controller
Kp_vq = 0.36;   % Proportional term q-axis voltage controller
Ki_vq = 850;    % Integral term q-axis voltage controller
