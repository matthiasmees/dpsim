% Export the same Simulink model as a macOS FMI 3.0 Co-Simulation FMU.
% This produces a sibling artifact for tools that support FMI doStep().
% DPsim's current importer in this branch is Model Exchange only.

fmiRoot = fileparts(mfilename('fullpath'));
modelFile = fullfile(fmiRoot, 'examples', 'inverter', 'sources', ...
    'GridTiedInverterPortedV3.slx');
sourceDir = fileparts(modelFile);
modelName = 'GridTiedInverterPortedV3';
assert(isfile(modelFile), 'Missing Simulink model: %s', modelFile);
assert(~verLessThan('matlab', '25.1'), 'Requires MATLAB R2025a or newer.');
assert(ismac && any(strcmp(computer('arch'), {'maca64', 'maci64'})), ...
    'Run this exporter on macOS arm64.');
assert(exist('exportToFMU', 'file') ~= 0, ...
    'Install or enable the MathWorks FMU Exporter first.');

staging = tempname;
mkdir(staging);
cleanup = onCleanup(@() rmdir(staging, 's')); %#ok<NASGU>
addpath(sourceDir);
pathCleanup = onCleanup(@() rmpath(sourceDir)); %#ok<NASGU>
load_system(modelFile);
exportToFMU(modelName, 'FMIVersion', '3.0', 'FMUType', 'CS', ...
    'SaveDirectory', staging);

fmu = fullfile(staging, [modelName '.fmu']);
assert(isfile(fmu) && dir(fmu).bytes > 0, 'Exporter did not create an FMU.');
target = fullfile(sourceDir, [modelName '.cs.fmu']);
copyfile(fmu, target, 'f');

inspection = tempname;
mkdir(inspection);
inspectionCleanup = onCleanup(@() rmdir(inspection, 's')); %#ok<NASGU>
unzip(target, inspection);
assert(isfile(fullfile(inspection, 'modelDescription.xml')), ...
    'FMU is missing modelDescription.xml.');
darwinBinary = dir(fullfile(inspection, 'binaries', 'aarch64-darwin', '*.dylib'));
assert(~isempty(darwinBinary), ...
    'FMU is missing binaries/aarch64-darwin/*.dylib.');
fprintf('Created macOS Co-Simulation FMU: %s\n', target);
fprintf('Binary: %s (%d bytes)\n', ...
    fullfile('binaries', 'aarch64-darwin', darwinBinary(1).name), ...
    darwinBinary(1).bytes);
close_system(modelName, 0);
