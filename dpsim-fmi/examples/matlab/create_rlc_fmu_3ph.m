function outputDir = create_rlc_fmu_3ph(outputDir, exportFMU)
%CREATE_RLC_FMU_3PH Three-phase EMT RLC benchmark and FMI 3.0 ME export.
% Run: create_rlc_fmu_3ph
% Or:  create_rlc_fmu_3ph('/Volumes/SSDMatthias/rlc_3ph_test')
% To generate models/references without exporting: create_rlc_fmu_3ph([], false)
%
% Requires MATLAB R2025a+, Simulink, the standalone FMU exporter and a
% supported C compiler. No Simscape or Specialized Power Systems is needed.
%
% Three identical, uncoupled series RLC branches, connected in grounded wye:
%   phase a -- R_a -- L_a -- C_a -- neutral/ground
%   phase b -- R_b -- L_b -- C_b -- neutral/ground
%   phase c -- R_c -- L_c -- C_c -- neutral/ground
%
% The source stays outside the FMU, as in the scalar benchmark.
% FMU voltage inputs: u_a, u_b, u_c [V], positive minus negative terminal.
% FMU current outputs: i_a, i_b, i_c [A], entering the positive terminal.
% Diagnostics: v_C_a, v_C_b, v_C_c [V], capacitor voltage in the same direction.
% MATLAB coordinates: x = [i_a; i_b; i_c; v_C_a; v_C_b; v_C_c].
% The generated FMU state order/names must still be read from its metadata.

if nargin < 2, exportFMU = true; end
if nargin < 1 || isempty(outputDir)
    stem = fullfile(pwd, ['rlc_3ph_fmu_test_' ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
    outputDir = stem;
    attempt = 1;
    while isfolder(outputDir)
        outputDir = sprintf('%s_%d', stem, attempt);
        attempt = attempt + 1;
    end
end
assert(~verLessThan('matlab', '25.1'), ...
    'FMI 3.0 Model Exchange export requires MATLAB R2025a or later.');
if exportFMU
    assert(exist('exportToFMU', 'file') ~= 0, ...
        'The standalone FMU export function exportToFMU is not installed.');
end
if ~isfolder(outputDir), mkdir(outputDir); end
[ok, info] = fileattrib(outputDir);
assert(ok, 'Cannot access output directory.');
outputDir = info.Name;
previousDir = pwd;
restoreDir = onCleanup(@() cd(previousDir)); %#ok<NASGU>
cd(outputDir);

plant = 'dpsim_rlc_3ph_me';
harness = 'dpsim_rlc_3ph_reference';
assert(~bdIsLoaded(plant) && ~bdIsLoaded(harness), ...
    'Close dpsim_rlc_3ph_me and dpsim_rlc_3ph_reference before rerunning.');
assert(~isfile([plant '.slx']) && ~isfile([harness '.slx']) && ...
    ~isfile([plant '.fmu']), ...
    'Benchmark files already exist here. Select a fresh output directory.');

% Same per-phase values as the scalar benchmark. Values are SI units.
I3 = eye(3);
R = 10 * I3;            % ohm; no mutual coupling
L = 20e-3 * I3;         % H
Ccap = 100e-6 * I3;     % F
VrmsPhase = 230;        % phase-to-neutral RMS voltage, not line-to-line
frequency = 50;         % Hz
Vpeak = sqrt(2) * VrmsPhase;
phaseRadians = [0, -2*pi/3, 2*pi/3]; % positive abc sequence
h = 10e-6;
stopTime = 0.1;
x0 = zeros(6, 1);
u0 = (Vpeak * sin(phaseRadians)).';

% L di_abc/dt = u_abc - R i_abc - v_C_abc; Ccap dv_C_abc/dt = i_abc.
A = [-(L \ R), -(L \ I3); Ccap \ I3, zeros(3)];
B = [L \ I3; zeros(3)];
Cout = eye(6);
Dout = zeros(6, 3);
Ccurrent = [I3, zeros(3)];
Dcurrent = zeros(3);
voltageNames = {'u_a', 'u_b', 'u_c'};
currentNames = {'i_a', 'i_b', 'i_c'};
outputNames = [currentNames, {'v_C_a', 'v_C_b', 'v_C_c'}];

load_system('simulink');
new_system(plant);
configureSolver(plant, stopTime, h);
add_block('simulink/Signal Routing/Mux', [plant '/Voltages_abc'], ...
    'Inputs', '3', 'Position', [125 80 130 220]);
for p = 1:3
    name = voltageNames{p};
    y = 90 + (p-1)*60;
    add_block('simulink/Sources/In1', [plant '/' name], ...
        'Port', num2str(p), 'PortDimensions', '1', ...
        'OutDataTypeStr', 'double', 'SampleTime', '0', ...
        'Position', [35 y 65 y+20]);
    add_line(plant, [name '/1'], ['Voltages_abc/' num2str(p)]);
end
add_block('simulink/Continuous/State-Space', [plant '/RLC_abc'], ...
    'A', mat2str(A, 17), 'B', mat2str(B, 17), ...
    'C', mat2str(Cout, 17), 'D', mat2str(Dout, 17), ...
    'X0', mat2str(x0, 17), 'Position', [200 110 330 190]);
add_block('simulink/Signal Routing/Demux', [plant '/Outputs_abc'], ...
    'Outputs', '6', 'Position', [390 30 395 310]);
add_line(plant, 'Voltages_abc/1', 'RLC_abc/1');
add_line(plant, 'RLC_abc/1', 'Outputs_abc/1');
for p = 1:6
    y = 35 + (p-1)*50;
    add_block('simulink/Sinks/Out1', [plant '/' outputNames{p}], ...
        'Port', num2str(p), 'Position', [465 y 495 y+20]);
    add_line(plant, ['Outputs_abc/' num2str(p)], [outputNames{p} '/1']);
end
save_system(plant, fullfile(outputDir, [plant '.slx']));

% MATLAB reference: balanced sinusoidal voltage source and the same dynamics.
new_system(harness);
configureSolver(harness, stopTime, h);
add_block('simulink/Signal Routing/Mux', [harness '/Voltages_abc'], ...
    'Inputs', '3', 'Position', [170 60 175 240]);
for p = 1:3
    sourceName = ['Source_' char('a' + p-1)];
    y = 65 + (p-1)*70;
    add_block('simulink/Sources/Sine Wave', [harness '/' sourceName], ...
        'Amplitude', num2str(Vpeak, 17), ...
        'Frequency', num2str(2*pi*frequency, 17), ...
        'Phase', num2str(phaseRadians(p), 17), ...
        'Bias', '0', 'SampleTime', '0', 'Position', [35 y 95 y+35]);
    add_line(harness, [sourceName '/1'], ['Voltages_abc/' num2str(p)]);
end
add_block([plant '/RLC_abc'], [harness '/RLC_abc'], ...
    'Position', [245 105 375 185]);
add_block('simulink/Sinks/To Workspace', [harness '/Record'], ...
    'VariableName', 'rlc_abc_output', 'SaveFormat', 'Timeseries', ...
    'MaxDataPoints', 'inf', 'Decimation', '1', ...
    'Position', [455 120 565 165]);
add_line(harness, 'Voltages_abc/1', 'RLC_abc/1');
add_line(harness, 'RLC_abc/1', 'Record/1');
save_system(harness, fullfile(outputDir, [harness '.slx']));

result = sim(harness, 'ReturnWorkspaceOutputs', 'on');
signal = result.get('rlc_abc_output');
referenceTime = signal.Time(:);
referenceOutput = reshape(signal.Data, [], 6);
referenceVoltage = phaseVoltages(referenceTime, Vpeak, frequency, phaseRadians);
columnNames = {'time_s', 'u_a_V', 'u_b_V', 'u_c_V', ...
    'i_a_A', 'i_b_A', 'i_c_A', 'v_C_a_V', 'v_C_b_V', 'v_C_c_V'};
writetable(array2table([referenceTime, referenceVoltage, referenceOutput], ...
    'VariableNames', columnNames), fullfile(outputDir, 'reference_simulink.csv'));

% SSN-style trapezoidal reference; the nonzero B/C starting voltages matter.
t = (0:round(stopTime/h)).' * h;
u = phaseVoltages(t, Vpeak, frequency, phaseRadians);
x = zeros(6, numel(t));
x(:, 1) = x0;
left = eye(6) - h*A/2;
right = eye(6) + h*A/2;
for k = 2:numel(t)
    x(:, k) = left \ (right*x(:, k-1) + h*B*(u(k-1, :) + u(k, :)).'/2);
end
y = (Cout*x + Dout*u.').';
writetable(array2table([t, u, y], 'VariableNames', columnNames), ...
    fullfile(outputDir, 'reference_trapezoidal.csv'));

% Balance checks test the phase wiring and initial-input handling.
assert(max(abs(sum(u, 2))) < 1e-9, 'Three-phase source is not balanced.');
assert(max(abs(sum(y(:, 1:3), 2))) < 1e-8, 'Phase currents are not balanced.');
assert(max(abs(sum(referenceOutput(:, 1:3), 2))) < 1e-6, ...
    'Simulink phase currents are not balanced.');

save(fullfile(outputDir, 'rlc_3ph_parameters.mat'), 'R', 'L', 'Ccap', ...
    'A', 'B', 'Cout', 'Dout', 'Ccurrent', 'Dcurrent', 'x0', 'u0', ...
    'h', 'stopTime', 'VrmsPhase', 'frequency', 'phaseRadians', ...
    'voltageNames', 'currentNames', 'outputNames');
writematrix(A, fullfile(outputDir, 'A.csv'));
writematrix(B, fullfile(outputDir, 'B.csv'));
writematrix(Ccurrent, fullfile(outputDir, 'C_current.csv'));
writematrix(Dcurrent, fullfile(outputDir, 'D_current.csv'));
writematrix(Cout, fullfile(outputDir, 'C_all_outputs.csv'));
writematrix(Dout, fullfile(outputDir, 'D_all_outputs.csv'));

mapping = struct('fmi_version', '3.0', 'interface', 'ModelExchange', ...
    'domain', 'EMT', 'frame', 'abc', 'connection', 'grounded_wye', ...
    'voltage_inputs', {voltageNames}, 'current_outputs', {currentNames}, ...
    'diagnostic_outputs', {{'v_C_a', 'v_C_b', 'v_C_c'}}, ...
    'terminal_order', {{'negative_abc', 'positive_abc'}}, ...
    'initial_voltage_V', u0.', 'initial_current_A', x0(1:3).', ...
    'matlab_state_order', {{'i_a', 'i_b', 'i_c', 'v_C_a', 'v_C_b', 'v_C_c'}}, ...
    'step_s', h, 'stop_time_s', stopTime, ...
    'phase_to_neutral_rms_V', VrmsPhase, 'frequency_Hz', frequency, ...
    'phase_angles_rad', phaseRadians);
writeText(fullfile(outputDir, 'port_mapping.json'), jsonencode(mapping));

figure('Name', 'Three-phase RLC benchmark');
subplot(3, 1, 1);
plot(t, u); grid on; ylabel('Voltage [V]'); legend('a', 'b', 'c');
title('Grounded-wye RLC: 230 V phase-to-neutral RMS, 50 Hz');
subplot(3, 1, 2);
plot(referenceTime, referenceOutput(:, 1:3), t, y(:, 1:3), '--');
grid on; ylabel('Current [A]');
legend('Simulink a', 'Simulink b', 'Simulink c', 'Trap a', 'Trap b', 'Trap c');
subplot(3, 1, 3);
plot(referenceTime, referenceOutput(:, 4:6), t, y(:, 4:6), '--');
grid on; xlabel('Time [s]'); ylabel('Capacitor voltage [V]');
savefig(fullfile(outputDir, 'rlc_3ph_reference.fig'));

fprintf('Created six-state three-phase models and references in:\n%s\n', outputDir);
fprintf('Initial input [u_a u_b u_c] = [%.9g %.9g %.9g] V\n', u0);
if exportFMU
    % Model Exchange: the importing DPsim SSN/MNA solver will integrate states.
    exportToFMU(plant, 'FMIVersion', '3.0', 'FMUType', 'ME', ...
        'SaveDirectory', outputDir);
    fprintf('Exported: %s\n', fullfile(outputDir, [plant '.fmu']));
else
    fprintf('FMU export skipped. To export the open model later:\n');
    fprintf('exportToFMU(''%s'', ''FMIVersion'', ''3.0'', ''FMUType'', ''ME'', ''SaveDirectory'', ''%s'');\n', ...
        plant, strrep(outputDir, '''', ''''''));
end
end

function configureSolver(model, stopTime, h)
set_param(model, 'StartTime', '0', 'StopTime', num2str(stopTime, 17), ...
    'SolverType', 'Variable-step', 'Solver', 'ode45', ...
    'MaxStep', num2str(h, 17), 'RelTol', '1e-9', 'AbsTol', '1e-10');
end

function values = phaseVoltages(time, amplitude, frequency, phase)
values = amplitude * sin(2*pi*frequency*time(:) + phase);
end

function writeText(path, text)
file = fopen(path, 'w');
assert(file ~= -1, 'Cannot write %s.', path);
closeFile = onCleanup(@() fclose(file)); %#ok<NASGU>
fprintf(file, '%s\n', text);
end
