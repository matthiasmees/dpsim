function run_smallgrid_sweep(mode, manifestPath, outputRoot, fmuPath, ~, modelPath)
% Build and run fresh, closed-loop Simulink solves for every manifest case.
cases = jsondecode(fileread(manifestPath));
if isempty(cases), error('Case manifest is empty.'); end
first = cases(1);
build_smallgrid_model(mode, fmuPath, modelPath, double(first.step_s), ...
    double(first.stop_time_s), double(first.load_step_s));
 [~,mdl]=fileparts(modelPath);
load_system(modelPath);
for k=1:numel(cases)
    entry=cases(k); dt=double(entry.step_s); stopTime=double(entry.stop_time_s);
    loadStepTime=double(entry.load_step_s);
    fprintf('Fresh closed-loop Simulink %s %d/%d: %s\n',mode,k,numel(cases),entry.label);
    simInput=Simulink.SimulationInput(mdl);
    simInput=simInput.setVariable('dt',dt);
    simInput=simInput.setVariable('stopTime',stopTime);
    simInput=simInput.setVariable('loadStepTime',loadStepTime);
    simInput=simInput.setModelParameter('FixedStep',num2str(dt,'%.17g'), ...
        'StopTime',num2str(stopTime,'%.17g'),'SimulationMode','normal');
    started=tic; result=sim(simInput); elapsed=toc(started);
    voltage=result.get('pcc_voltage'); current=result.get('inverter_current');
    if isempty(voltage)||isempty(current), error('Missing fresh run outputs for %s.',entry.label); end
    t=voltage.Time(:); vabc=squeeze(voltage.Data); iabc=squeeze(current.Data);
    if size(vabc,1)~=numel(t)&&size(vabc,2)==numel(t), vabc=vabc.'; end
    if size(iabc,1)~=numel(t)&&size(iabc,2)==numel(t), iabc=iabc.'; end
    if size(vabc,2)~=3||size(iabc,2)~=3
        error('Unexpected signal dimensions V=%s I=%s.',mat2str(size(vabc)),mat2str(size(iabc)));
    end
    if any(~isfinite(vabc),'all')||any(~isfinite(iabc),'all')
        error('Non-finite voltage or current in fresh Simulink run %s.',entry.label);
    end
    if abs(t(1))>1e-12||abs(t(end)-stopTime)>dt*1e-6||any(abs(diff(t)-dt)>max(1e-11,dt*1e-7))
        error('Simulink did not use the requested fixed-step grid for %s.',entry.label);
    end
    outDir=fullfile(outputRoot,entry.sweep_dir); if ~isfolder(outDir), mkdir(outDir); end
    save(fullfile(outDir,['simulink_' mode '.mat']),'t','vabc','iabc','-v7.3');
    trace=table(t,vabc(:,1),vabc(:,2),vabc(:,3),iabc(:,1),iabc(:,2),iabc(:,3), ...
        'VariableNames',{'time_s','v_a_V','v_b_V','v_c_V','i_a_A','i_b_A','i_c_A'});
    writetable(trace,fullfile(outDir,['simulink_' mode '.csv']));
    timer=struct('mode',mode,'label',entry.label, ...
        'network','Simulink ODE state equations for three-phase voltage source, PI feeder and PCC load', ...
        'inverter',ternary(strcmp(mode,'cs'),'FMI 3.0 Co-Simulation FMU Import block','native averaged inverter current state, solved with the network'), ...
        'closed_loop',true,'fresh_run',true,'network_solver','Simulink fixed-step ode14x, two Newton iterations', ...
        'step_s',dt,'stop_time_s',stopTime,'load_step_s',loadStepTime,'samples',numel(t), ...
        'simulink_wall_time_s',elapsed);
    fid=fopen(fullfile(outDir,['simulink_' mode '_timing.json']),'w'); cleanup=onCleanup(@()fclose(fid));
    fwrite(fid,jsonencode(timer,'PrettyPrint',true)); fwrite(fid,newline); clear cleanup;
end
save_system(mdl,modelPath); close_system(mdl,0);
fid=fopen(fullfile(outputRoot,[mode '_complete.marker']),'w');
fwrite(fid,'complete'); fclose(fid);
fprintf('Saved fresh closed-loop Simulink %s outputs under %s\n',mode,outputRoot);
end

function y=ternary(test,a,b)
if test, y=a; else, y=b; end
end
