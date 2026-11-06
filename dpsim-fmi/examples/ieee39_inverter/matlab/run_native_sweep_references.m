function run_native_sweep_references(manifestPath)
% Run matched native Simulink replay references for an input/output manifest.
cases = jsondecode(fileread(manifestPath));
for index = 1:numel(cases)
    entry = cases(index);
    fprintf('Native reference %d/%d: %s\n', index, numel(cases), entry.label);
    % Compile/cache Accelerator once for this communication step before timing.
    run_native_replay(char(entry.input_csv), char(entry.output_csv), 0.05, false);
    run_native_replay(char(entry.input_csv), char(entry.output_csv), ...
        double(entry.stop_time_s), true);
end
end
