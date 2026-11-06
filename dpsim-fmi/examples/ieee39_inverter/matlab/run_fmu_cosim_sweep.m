function run_fmu_cosim_sweep(manifestPath, outputRoot, modelPath)
% Run the same FMI 3 CS FMU beside the native inverter in one Simulink model.
% Both branches receive the paired DPsim CS terminal-voltage trace.
cases = jsondecode(fileread(manifestPath));
fmuFolder = fileparts(cases(1).fmu_path);
addpath(fmuFolder);
mdl = 'IEEE39_Simulink_FMU_Cosim_Validation';
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
if isfile(modelPath)
    delete(modelPath);
end
buildHarness(mdl, cases(1).fmu_path);
save_system(mdl, modelPath);

for index = 1:numel(cases)
    entry = cases(index);
    fprintf('Simulink native + FMU CS %d/%d: %s\n', ...
        index, numel(cases), entry.label);
    T = readtable(entry.input_csv, 'VariableNamingRule', 'preserve');
    loggedTime = T.time;
    vabc = [T.('inverter.v_0'), T.('inverter.v_1'), T.('inverter.v_2')];
    dt = double(entry.step_s);
    gridTime = (0:numel(loggedTime)-1)' * dt;
    if any(abs(loggedTime-gridTime) > max(5.1e-6, 1e-9*max(1,dt)))
        error('DPsim voltage trace is not a rounded uniform grid: %s', entry.label);
    end
    stopTime = double(entry.stop_time_s);
    if abs(gridTime(end)-stopTime) > dt/2
        if abs(stopTime-gridTime(end)-dt) > max(5.1e-6, 1e-9*max(1,dt))
            error('DPsim trace does not reach the requested stop time: %s', entry.label);
        end
        % DPsim omits the terminal log row for a few steps. Add one linearly
        % extrapolated PCC sample so both Simulink branches run to stopTime.
        vabc(end+1,:) = 2*vabc(end,:) - vabc(end-1,:);
        gridTime(end+1,1) = stopTime;
    elseif abs(gridTime(end)-stopTime) > max(5.1e-6, 1e-9*max(1,dt))
        error('DPsim trace endpoint is not aligned: %s', entry.label);
    end

    i0 = localPassiveCurrent(vabc(1,:));
    simInput = Simulink.SimulationInput(mdl);
    simInput = simInput.setVariable('vabc', timeseries(vabc, gridTime));
    simInput = simInput.setVariable('i0', i0);
    simInput = simInput.setVariable('dt', dt);
    simInput = simInput.setModelParameter('FixedStep', ...
        num2str(dt, '%.17g'), 'StopTime', num2str(stopTime, '%.17g'), ...
        'SimulationMode', 'normal');
    started = tic;
    simOutput = sim(simInput);
    elapsed = toc(started);
    native = simOutput.get('i_native');
    cs = simOutput.get('i_fmu_cs');
    nativeData = squeeze(native.Data);
    csData = squeeze(cs.Data);
    if size(nativeData, 2) ~= 3 && size(nativeData, 1) == 3
        nativeData = nativeData.';
    end
    if size(csData, 2) ~= 3 && size(csData, 1) == 3
        csData = csData.';
    end
    if size(nativeData, 2) ~= 3 || size(csData, 2) ~= 3 || ...
            size(nativeData, 1) ~= size(csData, 1)
        error('Unexpected Simulink outputs for %s: native=%s FMU-CS=%s', ...
            entry.label, mat2str(size(nativeData)), mat2str(size(csData)));
    end
    t = native.Time;
    if numel(t) ~= size(nativeData,1) || any(abs(cs.Time-t) > dt*1e-5)
        error('Native and FMI CS output time grids differ for %s', entry.label);
    end
    result = table(t, vabc(1:numel(t),1), vabc(1:numel(t),2), ...
        vabc(1:numel(t),3), nativeData(:,1), nativeData(:,2), ...
        nativeData(:,3), csData(:,1), csData(:,2), csData(:,3), ...
        'VariableNames', {'time','v_a','v_b','v_c', ...
        'native_i_a','native_i_b','native_i_c', ...
        'fmu_cs_i_a','fmu_cs_i_b','fmu_cs_i_c'});
    outDir = fullfile(outputRoot, entry.sweep_dir, entry.case_dir);
    if ~isfolder(outDir), mkdir(outDir); end
    writetable(result, fullfile(outDir, 'simulink_cosim.csv'));
    timing = struct('label', entry.label, 'host', 'Simulink FMU Import block', ...
        'fmi_mode', 'FMI 3.0 Co-Simulation', 'communication_step_s', dt, ...
        'stop_time_s', stopTime, 'sample_count', height(result), ...
        'simulink_run_wall_time_s', elapsed, 'native_and_fmu_run_together', true, ...
        'input_voltage_source', 'paired DPsim CS inverter terminal trace');
    fid = fopen(fullfile(outDir, 'simulink_cosim_timing.json'), 'w');
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, jsonencode(timing, 'PrettyPrint', true));
    fwrite(fid, newline);
    clear cleanup;
end
save_system(mdl, modelPath);
close_system(mdl, 0);
fprintf('Saved Simulink harness: %s\n', modelPath);
end

function buildHarness(mdl, fmuPath)
assignin('base', 'dt', 50e-6);
assignin('base', 'stopTime', 0.05);
assignin('base', 'i0', [0;0;0]);
assignin('base', 'vabc', timeseries(zeros(2,3), [0;0.05]));
new_system(mdl);
set_param(mdl, 'Solver', 'ode4', 'SolverType', 'Fixed-step', ...
    'FixedStep', 'dt', 'StopTime', 'stopTime', 'SaveTime', 'off', ...
    'SaveOutput', 'off', 'ReturnWorkspaceOutputs', 'on');

add_block('simulink/Sources/From Workspace', [mdl '/PCC voltage'], ...
    'VariableName', 'vabc', 'Interpolate', 'on', 'SampleTime', 'dt', ...
    'Position', [35 185 150 215]);
add_block('simulink/Signal Routing/Demux', [mdl '/Split PCC voltage'], ...
    'Outputs', '3', 'Position', [185 170 190 230]);
add_block('simulink/User-Defined Functions/MATLAB Function', ...
    [mdl '/Native inverter current controller'], ...
    'Position', [295 135 505 240]);
chart = sfroot().find('-isa', 'Stateflow.EMChart', ...
    'Path', [mdl '/Native inverter current controller']);
chart.Script = sprintf([ ...
    'function [dia,dib,dic] = fcn(va,vb,vc,ia0,ib0,ic0)\n' ...
    '%%#codegen\n' ...
    'P=125e6; Q=80.881e6; Vnom=345e3; Ipu=1.20; tau=5e-3;\n' ...
    'a=(2/3)*(va-0.5*vb-0.5*vc);\n' ...
    'b=(2/3)*(sqrt(3)/2)*(vb-vc); v2=a*a+b*b;\n' ...
    'if v2<1, iaref=0; ibref=0; icref=0; else\n' ...
    ' k=2/(3*v2); ialpha=k*(P*a+Q*b); ibeta=k*(P*b-Q*a);\n' ...
    ' imax=Ipu*sqrt(2)*hypot(P,Q)/(sqrt(3)*Vnom); imag=hypot(ialpha,ibeta);\n' ...
    ' if imag>imax, ialpha=ialpha*imax/imag; ibeta=ibeta*imax/imag; end\n' ...
    ' iaref=ialpha; ibref=-0.5*ialpha+sqrt(3)/2*ibeta;\n' ...
    ' icref=-0.5*ialpha-sqrt(3)/2*ibeta; end\n' ...
    ' dia=(-iaref-ia0)/tau; dib=(-ibref-ib0)/tau; dic=(-icref-ic0)/tau;\n' ...
    'end\n']);

for phase = 1:3
    y = 80 + 55*phase;
    add_block('simulink/Continuous/Integrator', ...
        [mdl '/Native current state ' char('a'+phase-1)], ...
        'InitialCondition', sprintf('i0(%d)', phase), ...
        'Position', [565 y 600 y+30]);
end
add_block('simulink/Signal Routing/Mux', [mdl '/Native current abc'], ...
    'Inputs', '3', 'Position', [665 205 670 265]);
add_block('simulink/Signal Routing/Mux', [mdl '/FMI CS current abc'], ...
    'Inputs', '3', 'Position', [665 390 670 450]);
add_block('simulink/Sinks/To Workspace', [mdl '/Log native current'], ...
    'VariableName', 'i_native', 'SaveFormat', 'Timeseries', ...
    'Position', [725 215 835 245]);
add_block('simulink/Sinks/To Workspace', [mdl '/Log FMI CS current'], ...
    'VariableName', 'i_fmu_cs', 'SaveFormat', 'Timeseries', ...
    'Position', [725 400 835 430]);

addpath(fileparts(fmuPath));
add_block('simulink_extras/FMU Import/FMU', [mdl '/IEEE39 averaged inverter FMU CS'], ...
    'FMUName', 'IEEE39AvgInverter.fmu', 'FMUMode', 'CoSimulation', ...
    'FMUInputMapping', 'Flat', 'FMUOutputMapping', 'Flat', ...
    'Position', [350 340 500 435]);
set_param(mdl, 'SimulationCommand', 'update');

from = [mdl '/PCC voltage']; split = [mdl '/Split PCC voltage'];
nativeCtl = [mdl '/Native inverter current controller'];
fmu = [mdl '/IEEE39 averaged inverter FMU CS'];
add_line(mdl, 'PCC voltage/1', 'Split PCC voltage/1');
for phase = 1:3
    add_line(mdl, sprintf('Split PCC voltage/%d',phase), ...
        sprintf('Native inverter current controller/%d',phase));
    add_line(mdl, sprintf('Split PCC voltage/%d',phase), ...
        sprintf('IEEE39 averaged inverter FMU CS/%d',phase));
    state = sprintf('Native current state %c', char('a'+phase-1));
    add_line(mdl, [state '/1'], ...
        sprintf('Native inverter current controller/%d',phase+3));
    add_line(mdl, sprintf('Native inverter current controller/%d',phase), ...
        [state '/1']);
    add_line(mdl, [state '/1'], sprintf('Native current abc/%d',phase));
    add_line(mdl, sprintf('IEEE39 averaged inverter FMU CS/%d',phase), ...
        sprintf('FMI CS current abc/%d',phase));
end
add_line(mdl, 'Native current abc/1', 'Log native current/1');
add_line(mdl, 'FMI CS current abc/1', 'Log FMI CS current/1');
end

function current = localPassiveCurrent(vabc)
P = 125e6; Q = 80.881e6; Vnom = 345e3; Ipu = 1.20;
alpha = (2/3)*(vabc(1)-0.5*vabc(2)-0.5*vabc(3));
beta = (2/3)*(sqrt(3)/2)*(vabc(2)-vabc(3));
v2 = alpha*alpha + beta*beta;
if v2 < 1, current = [0;0;0]; return; end
k = 2/(3*v2); ialpha = k*(P*alpha+Q*beta); ibeta = k*(P*beta-Q*alpha);
imax = Ipu*sqrt(2)*hypot(P,Q)/(sqrt(3)*Vnom); imag = hypot(ialpha,ibeta);
if imag > imax, ialpha=ialpha*imax/imag; ibeta=ibeta*imax/imag; end
current = [-ialpha; 0.5*ialpha-sqrt(3)/2*ibeta; ...
    0.5*ialpha+sqrt(3)/2*ibeta];
end
