function export_linear_rl_fmus(outputDirectory)
% Export the same continuous RL model as FMI 2.0 Model Exchange and CS FMUs.
% Requires MATLAB/Simulink and the Simulink FMU Exporter (exportToFMU).
arguments
    outputDirectory (1,1) string = fullfile(pwd, "build-fmi", "linear_rl_stability", "fmus")
end

assert(exist('exportToFMU', 'file') ~= 0, ...
    'The Simulink FMU Exporter (exportToFMU) is not installed.');
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
sourceDirectory = fileparts(fileparts(mfilename('fullpath')));
if ~isfolder(fullfile(sourceDirectory, 'sources'))
    mkdir(fullfile(sourceDirectory, 'sources'));
end
modelName = 'LinearRLStability';
modelFile = fullfile(sourceDirectory, 'sources', [modelName '.slx']);
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
new_system(modelName);
add_block('simulink/Ports & Subsystems/In1', [modelName '/u'], ...
    'Position', [30 95 60 115]);
add_block('simulink/Continuous/Integrator', [modelName '/Current'], ...
    'InitialCondition', '0', 'Position', [330 90 360 120]);
add_block('simulink/Math Operations/Gain', [modelName '/Rf'], ...
    'Gain', '0.1', 'Position', [110 150 155 180]);
add_block('simulink/Math Operations/Sum', [modelName '/VoltageDrop'], ...
    'Inputs', '+-', 'Position', [205 85 230 125]);
add_block('simulink/Math Operations/Gain', [modelName '/InverseL'], ...
    'Gain', '1/20e-6', 'Position', [260 85 300 125]);
add_block('simulink/Ports & Subsystems/Out1', [modelName '/i'], ...
    'Position', [420 95 450 115]);
add_line(modelName, 'u/1', 'VoltageDrop/1');
add_line(modelName, 'Current/1', 'Rf/1');
add_line(modelName, 'Rf/1', 'VoltageDrop/2');
add_line(modelName, 'VoltageDrop/1', 'InverseL/1');
add_line(modelName, 'InverseL/1', 'Current/1');
add_line(modelName, 'Current/1', 'i/1');
set_param(modelName, 'SolverType', 'Variable-step', 'Solver', 'ode15s', ...
    'StopTime', '0.1', 'SaveFormat', 'Dataset', 'ReturnWorkspaceOutputs', 'off');
save_system(modelName, modelFile);

exportToFMU(modelName, 'FMIVersion', '2.0', 'FMUType', 'ME', ...
    'FMUName', 'LinearRL_ME', 'SaveDirectory', outputDirectory, ...
    'ExportedInputNames', {'u'}, 'ExportedOutputNames', {'i'}, ...
    'SaveSourceCodeToFMU', 'on');
exportToFMU(modelName, 'FMIVersion', '2.0', 'FMUType', 'CS', ...
    'FMUName', 'LinearRL_CS', 'SaveDirectory', outputDirectory, ...
    'ExportedInputNames', {'u'}, 'ExportedOutputNames', {'i'}, ...
    'SaveSourceCodeToFMU', 'on');
assert(isfile(fullfile(outputDirectory, 'LinearRL_ME.fmu')));
assert(isfile(fullfile(outputDirectory, 'LinearRL_CS.fmu')));
copyfile(modelFile, fullfile(sourceDirectory, 'sources', [modelName '.slx']), 'f');
for mode = ["ME", "CS"]
    archive = fullfile(outputDirectory, "LinearRL_" + mode + ".fmu");
    copyfile(archive, fullfile(sourceDirectory, 'sources', "LinearRL_" + mode + ".fmu"), 'f');
    metadata = fullfile(sourceDirectory, 'sources', "modelDescription_" + mode + ".xml");
    entries = unzip(archive, tempdir);
    cleanup = onCleanup(@() delete(fullfile(tempdir, 'modelDescription.xml')));
    copyfile(fullfile(tempdir, 'modelDescription.xml'), metadata, 'f');
    clear cleanup
end
fprintf('Exported FMI 2.0 ME and CS archives to %s\n', outputDirectory);
close_system(modelName, 0);
end
