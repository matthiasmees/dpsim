% SPDX-License-Identifier: MPL-2.0
% Native Simulink state-space reference for the dual-interface DPsim test grid.
scriptPath = mfilename('fullpath');
caseDir = fileparts(fileparts(scriptPath));
repoDir = fileparts(fileparts(fileparts(fileparts(caseDir))));
buildDir = fullfile(repoDir, 'build-fmi-cosim', 'dual-interface-rl', 'native-simulink');
if ~exist(buildDir, 'dir')
    mkdir(buildDir);
end

frequency = 50;
phasePeak = sqrt(2) * 230;
sourceR = 2;
sourceL = 0.01;
loadR = 100;
branchR = 8;
branchL = 0.02;
sagStart = 0.04;
sagEnd = 0.06;
stopTime = 0.08;
step = 2.5e-6;
t = (0:step:stopTime)';
envelope = min(1, max(0, t / 0.01));
sag = ones(size(t));
sag(t >= sagStart & t < sagEnd) = 0.75;
sourceVoltage = zeros(numel(t), 3);
for phase = 1:3
    sourceVoltage(:, phase) = envelope .* sag .* phasePeak .* ...
        sin(2*pi*frequency*t - (phase-1)*2*pi/3);
end
sourceVoltageTs = timeseries(sourceVoltage, t);

% Per-phase states are [grid current; imported branch current].
Aphase = [-(sourceR + loadR)/sourceL, loadR/sourceL; ...
          loadR/branchL, -(loadR + branchR)/branchL];
Bphase = [1/sourceL; 0];
A = kron(eye(3), Aphase);
B = kron(eye(3), Bphase);
C = zeros(9, 6);
for phase = 1:3
    stateRows = 2*phase-1:2*phase;
    C(phase, stateRows) = [loadR, -loadR];
    C(3+phase, stateRows) = [0, 1];
    C(6+phase, stateRows) = [1, 0];
end
D = zeros(9, 3);
x0 = zeros(6, 1);

model = 'DualInterfaceRL_Native';
modelFile = fullfile(buildDir, [model '.slx']);
if bdIsLoaded(model)
    close_system(model, 0);
end
new_system(model);
open_system(model);
add_block('simulink/Sources/From Workspace', [model '/Grid voltage'], ...
    'VariableName', 'sourceVoltageTs', 'Interpolate', 'on', ...
    'Position', [40 90 170 130]);
add_block('simulink/Continuous/State-Space', [model '/R-L grid and load'], ...
    'A', 'A', 'B', 'B', 'C', 'C', 'D', 'D', 'X0', 'x0', ...
    'Position', [240 80 400 140]);
add_block('simulink/Sinks/To Workspace', [model '/Waveforms'], ...
    'VariableName', 'ySimulink', 'SaveFormat', 'Structure With Time', ...
    'Position', [480 90 610 130]);
add_line(model, 'Grid voltage/1', 'R-L grid and load/1');
add_line(model, 'R-L grid and load/1', 'Waveforms/1');
set_param(model, 'StopTime', num2str(stopTime, '%.12g'), ...
    'SolverType', 'Fixed-step', 'Solver', 'ode4', ...
    'FixedStep', num2str(step, '%.12g'), 'ReturnWorkspaceOutputs', 'on', ...
    'SaveTime', 'off', 'SaveOutput', 'off');
save_system(model, modelFile);

timer = tic;
simOut = sim(model, 'ReturnWorkspaceOutputs', 'on');
simulationWall = toc(timer);
waveforms = simOut.get('ySimulink');
y = waveforms.signals.values;
out = table(waveforms.time, y(:,1), y(:,2), y(:,3), ...
    y(:,4), y(:,5), y(:,6), y(:,7), y(:,8), y(:,9), ...
    'VariableNames', {'time','pcc_v_a','pcc_v_b','pcc_v_c', ...
    'fmu_i_a','fmu_i_b','fmu_i_c','grid_i_a','grid_i_b','grid_i_c'});
timer = tic;
writetable(out, fullfile(buildDir, 'native_reference.csv'));
csvWriteWall = toc(timer);
info = struct('simulation_wall_s', simulationWall, ...
    'csv_write_wall_s', csvWriteWall, 'step_s', step, ...
    'stop_time_s', stopTime, 'samples', height(out), ...
    'solver', 'Simulink fixed-step ode4', 'source_R_ohm', sourceR, ...
    'source_L_H', sourceL, 'load_R_ohm', loadR, ...
    'branch_R_ohm', branchR, 'branch_L_H', branchL, ...
    'sag_start_s', sagStart, 'sag_end_s', sagEnd, 'sag_factor', 0.75);
fid = fopen(fullfile(buildDir, 'timing.json'), 'w');
fwrite(fid, jsonencode(info, 'PrettyPrint', true));
fclose(fid);
close_system(model, 0);
disp(info);
