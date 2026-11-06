function run_ieee39_closed_loop(manifestPath, outputRoot, fmuPath, libraryPath, modelPath)
% Build and run short validation or matched IEEE39 timestep sweep cases.
cases = jsondecode(fileread(manifestPath));
if isempty(cases), error('The case manifest is empty.'); end
matlabDir = fileparts(mfilename('fullpath'));
fmuFolder = fileparts(fmuPath);
addpath(fmuFolder);
dt0 = double(cases(1).step_s);
stop0 = double(cases(1).stop_time_s);
step0 = double(cases(1).load_step_s);

if ~isfile(fullfile(matlabDir, ['ieee39_network_sfun.' mexext]))
    error('Build ieee39_network_sfun.%s with build_ieee39_network_sfun.sh first.', mexext);
end
build_ieee39_closed_loop(fmuPath, libraryPath, modelPath, dt0, stop0, step0);
mdl = 'IEEE39_ClosedLoop';
load_system(modelPath);

for index = 1:numel(cases)
    entry = cases(index);
    dt = double(entry.step_s);
    stopTime = double(entry.stop_time_s);
    loadStepTime = double(entry.load_step_s);
    fprintf('Closed-loop IEEE39 %d/%d: %s\n', index, numel(cases), entry.label);
    simInput = Simulink.SimulationInput(mdl);
    simInput = simInput.setVariable('dt', dt);
    simInput = simInput.setVariable('stopTime', stopTime);
    simInput = simInput.setVariable('loadStepTime', loadStepTime);
    simInput = simInput.setModelParameter('FixedStep', num2str(dt, '%.17g'), ...
        'StopTime', num2str(stopTime, '%.17g'), 'SimulationMode', 'normal');
    started = tic;
    result = sim(simInput);
    elapsed = toc(started);
    bus = result.get('bus_voltage');
    current = result.get('inverter_current');
    timing = result.get('network_timing');
    if isempty(bus) || isempty(current) || isempty(timing)
        error('Missing logged outputs for %s.', entry.label);
    end
    t = bus.Time(:);
    voltage = squeeze(bus.Data);
    iabc = squeeze(current.Data);
    timingData = squeeze(timing.Data);
    if size(voltage, 1) ~= numel(t) && size(voltage, 2) == numel(t), voltage = voltage.'; end
    if size(iabc, 1) ~= numel(t) && size(iabc, 2) == numel(t), iabc = iabc.'; end
    if size(timingData, 1) ~= numel(t) && size(timingData, 2) == numel(t), timingData = timingData.'; end
    if size(voltage,2) ~= 119 || size(iabc,2) ~= 3 || size(timingData,2) ~= 2
        error('Unexpected log dimensions: V=%s, I=%s, timing=%s.', ...
            mat2str(size(voltage)), mat2str(size(iabc)), mat2str(size(timingData)));
    end
    if abs(t(1)) > 1e-12 || abs(t(end) - stopTime) > dt * 1e-6 || ...
            any(abs(diff(t) - dt) > max(1e-11, dt * 1e-7))
        error('Simulink did not use the requested time grid for %s.', entry.label);
    end
    outDir = fullfile(outputRoot, entry.sweep_dir, entry.case_dir);
    if ~isfolder(outDir), mkdir(outDir); end
    save(fullfile(outDir, 'simulink_closed_loop.mat'), ...
        't', 'voltage', 'iabc', 'timingData', '-v7.3');
    caseTiming = struct('label', entry.label, ...
        'network', 'DPsim EMT IEEE39 MNA (same C++ topology/parameters)', ...
        'inverter', 'FMI 3.0 Co-Simulation FMU Import block', ...
        'closed_loop', true, 'communication_step_s', dt, ...
        'stop_time_s', stopTime, 'load_step_s', loadStepTime, ...
        'samples', numel(t), 'simulink_wall_time_s', elapsed, ...
        'network_initialization_s', max(timingData(:,2)), ...
        'network_step_mean_s', mean(timingData(:,1)), ...
        'network_step_max_s', max(timingData(:,1)), ...
        'network_step_total_s', sum(timingData(:,1)));
    fid = fopen(fullfile(outDir, 'simulink_closed_loop_timing.json'), 'w');
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, jsonencode(caseTiming, 'PrettyPrint', true)); fwrite(fid, newline);
    clear cleanup;
    % Compact trace for matching DPsim output comparison and notebook use.
    compact = table(t, voltage(:,88), voltage(:,89), voltage(:,90), ...
        voltage(:,115), voltage(:,116), voltage(:,117), iabc(:,1), ...
        iabc(:,2), iabc(:,3), timingData(:,1), ...
        'VariableNames', {'time','bus30_v_a','bus30_v_b','bus30_v_c', ...
        'bus39_v_a','bus39_v_b','bus39_v_c','i_a','i_b','i_c', ...
        'network_step_wall_s'});
    writetable(compact, fullfile(outDir, 'simulink_closed_loop_pcc.csv'));
end
save_system(mdl, modelPath);
close_system(mdl, 0);
fprintf('Closed-loop results saved under %s\n', outputRoot);
end
