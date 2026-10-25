% Regenerate the physical Simulink traces used by run_paper_validation.py.
% Run this script from MATLAB. It copies each checked-in model to build-fmi,
% uses a fixed 1 us Simscape local solver step, and writes CSVs in the existing
% build/fmi-validation/matlab-reference-attempt result directory. The checked-
% in source models and FMU assets are never modified.

scriptFolder = fileparts(mfilename('fullpath'));
caseFolder = fileparts(fileparts(scriptFolder));
repositoryRoot = fileparts(fileparts(fileparts(caseFolder)));
out = fullfile(repositoryRoot, 'build-fmi/fmi-validation/matlab-reference-attempt');
if ~isfolder(out)
    mkdir(out);
end
addpath(fullfile(caseFolder, 'sources'), scriptFolder);
Simulink.fileGenControl('set', 'CacheFolder', fullfile(out, 'cache'), ...
    'CodeGenFolder', fullfile(out, 'codegen'), 'createDir', true);
maxNumCompThreads(1);
h = 1e-6;
models = {
    fullfile(scriptFolder, 'balanced_load', ...
        'MixedReference_balanced_load_PhysicalEvent.slx'), 'balanced_load';
    fullfile(scriptFolder, 'unbalanced_load', ...
        'MixedReference_unbalanced_load_PhysicalEvent.slx'), 'unbalanced_load';
    fullfile(scriptFolder, 'voltage_sag', ...
        'MixedReference_voltage_sag_GridCorrection.slx'), 'voltage_sag'
};
labels = {'time','i_grid_a','v_pcc_a','i_grid_b','v_pcc_b','i_grid_c', ...
    'v_pcc_c','v_d_pu','v_q_pu','pi_d','pi_q','pll_frequency','pll_theta', ...
    'pll_magnitude'};

for k = 1:size(models, 1)
    sourceModel = models{k, 1};
    tag = models{k, 2};
    [~, modelName] = fileparts(sourceModel);
    localModel = fullfile(out, [modelName '.slx']);
    copyfile(sourceModel, localModel, 'f');
    load_system(localModel);
    set_param(modelName, 'FixedStep', num2str(h, 17), 'StopTime', '0.2');
    solverConfig = find_system(modelName, 'SearchDepth', 1, ...
        'LookUnderMasks', 'all', 'RegExp', 'on', 'Name', 'Solver.*Configuration');
    assert(~isempty(solverConfig), ['No Simscape Solver Configuration block in ' tag]);
    set_param(solverConfig{1}, 'LocalSolverSampleTime', num2str(h, 17));
    fprintf('Simulating %s at 1 us through 0.2 s...\n', tag);
    result = sim(modelName, 'ReturnWorkspaceOutputs', 'on', 'SaveOutput', 'on', ...
        'OutputSaveName', 'yout', 'SaveFormat', 'Dataset', 'LimitDataPoints', 'off');
    assert(result.yout.numElements == 13, ['Unexpected signals in ' tag]);
    t = (0:200000)' * h;
    data = t;
    for j = 1:result.yout.numElements
        signal = result.yout.getElement(j);
        method = 'linear';
        if j >= 7
            method = 'previous';
        end
        values = squeeze(signal.Values.Data);
        data(:, end+1) = interp1(signal.Values.Time, values, t, method, 'extrap'); %#ok<SAGROW>
    end
    filename = fullfile(out, [tag '_reference_trap_1us.csv']);
    writetable(array2table(data, 'VariableNames', labels), filename);
    fprintf('Wrote %s\n', filename);
    close_system(modelName, 0);
end
