% Export the repository's supported Simulink source models as FMI 3.0 ME FMUs.
% Run from Fedora/Linux with MATLAB R2025a and the standalone FMU exporter:
%   matlab -batch "run('/home/matthiasmees/dpsim/dpsim-fmi/export_all_simulink_fmus.m')"
%
% Linux FMUs and extracted metadata are written beside their source models as
% <model>.x86_64-linux.fmu and modelDescription.x86_64-linux.xml. The existing
% unsuffixed macOS FMUs and modelDescription.xml files remain untouched.
%
% This is an explicit allowlist. Reference harnesses and inverter validation
% models under matlab/ are not standalone FMU source models.

fmiRoot = fileparts(mfilename('fullpath'));
stagingRoot = tempname;
mkdir(stagingRoot);
removeStagingRoot = onCleanup(@() rmdir(stagingRoot, 's')); %#ok<NASGU>

% MATLAB prepends its bundled C++ runtime to LD_LIBRARY_PATH. On current
% Fedora releases this can make the system GCC fail before compilation because
% GCC needs newer GLIBCXX/CXXABI symbols. Keep MATLAB running, but let compiler
% subprocesses resolve the system C++ runtime.
originalLdLibraryPath = getenv('LD_LIBRARY_PATH');
restoreLdLibraryPath = onCleanup(@() setenv( ...
    'LD_LIBRARY_PATH', originalLdLibraryPath)); %#ok<NASGU>
matlabRuntimeDir = fullfile(matlabroot, 'sys', 'os', 'glnxa64');
libraryPaths = strsplit(originalLdLibraryPath, pathsep);
libraryPaths(strcmp(libraryPaths, matlabRuntimeDir)) = [];
setenv('LD_LIBRARY_PATH', strjoin(libraryPaths, pathsep));

previousDir = pwd;
restoreDir = onCleanup(@() cd(previousDir)); %#ok<NASGU>
cd(stagingRoot);
models = {
    'examples/linear_scalar/sources/dpsim_rlc_me.slx'
    'examples/linear_three_phase/sources/dpsim_rlc_3ph_me.slx'
    'examples/synchronous_generator/sources/dpsim_sg_3ph_me.slx'
    'examples/two_generators/sources/dpsim_sg_3ph_second.slx'
    'examples/inverter/sources/GridTiedInverterPortedV3.slx'
};

assert(~verLessThan('matlab', '25.1'), ...
    'This batch export requires MATLAB R2025a or later.');
assert(strcmp(computer('arch'), 'glnxa64'), ...
    'This script creates the x86_64 Linux FMU variants and must run on Linux.');
assert(exist('exportToFMU', 'file') ~= 0, ...
    ['The standalone FMU exporter is unavailable. Install/enable the ' ...
     'MathWorks FMU Exporter, then retry.']);

for k = 1:numel(models)
    modelFile = fullfile(fmiRoot, models{k});
    assert(isfile(modelFile), 'Missing Simulink model: %s', modelFile);
end

failed = {};
for k = 1:numel(models)
    modelFile = fullfile(fmiRoot, models{k});
    [sourceDir, modelName] = fileparts(modelFile);
    outputDir = tempname(stagingRoot);
    mkdir(outputDir);
    addpath(sourceDir);
    removeSourcePath = onCleanup(@() rmpath(sourceDir)); %#ok<NASGU>
    wasLoaded = bdIsLoaded(modelName);
    fprintf('\n[%d/%d] Exporting %s\n', k, numel(models), modelFile);

    try
        load_system(modelFile);
        exportToFMU(modelName, ...
            'FMIVersion', '3.0', ...
            'FMUType', 'ME', ...
            'SaveDirectory', outputDir);

        fmuFile = fullfile(outputDir, [modelName '.fmu']);
        assert(isfile(fmuFile) && dir(fmuFile).bytes > 0, ...
            'Exporter returned without creating a nonempty FMU: %s', fmuFile);
        installedFmu = fullfile(sourceDir, ...
            [modelName '.x86_64-linux.fmu']);
        copyfile(fmuFile, installedFmu, 'f');
        extractDir = tempname(outputDir);
        mkdir(extractDir);
        removeExtractDir = onCleanup(@() rmdir(extractDir, 's')); %#ok<NASGU>
        unzip(fmuFile, extractDir);
        copyfile(fullfile(extractDir, 'modelDescription.xml'), ...
            fullfile(sourceDir, 'modelDescription.x86_64-linux.xml'), 'f');
        clear removeExtractDir;
        assert(isfile(fullfile(sourceDir, ...
            'modelDescription.x86_64-linux.xml')), ...
            'FMU is missing modelDescription.xml: %s', installedFmu);
        fprintf('Created %s (%g bytes)\n', ...
            installedFmu, dir(installedFmu).bytes);
    catch exception
        fprintf(2, 'FAILED %s: %s\n', modelFile, exception.message);
        failed{end + 1} = modelFile; %#ok<SAGROW>
    end

    if ~wasLoaded && bdIsLoaded(modelName)
        close_system(modelName, 0);
    end
    clear removeSourcePath;
end

if ~isempty(failed)
    fprintf(2, '\nFMU export failed for %d of %d model(s):\n', ...
        numel(failed), numel(models));
    fprintf(2, '  %s\n', failed{:});
    error('dpsim-fmi:BatchExportFailed', ...
        'One or more FMU exports failed; see the per-model diagnostics above.');
end

% The two-generator example intentionally reuses the first generator export
% under the established first_generator name.
firstGeneratorSource = fullfile(fmiRoot, ...
    'examples/synchronous_generator/sources/dpsim_sg_3ph_me.x86_64-linux.fmu');
firstGeneratorTarget = fullfile(fmiRoot, ...
    'examples/two_generators/sources/first_generator.x86_64-linux.fmu');
copyfile(firstGeneratorSource, firstGeneratorTarget, 'f');

fprintf('\nExported all %d FMUs successfully.\n', numel(models));
