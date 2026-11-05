function outputDir = create_sg_grid_fmu_3ph(outputDir, exportFMU)
%CREATE_SG_GRID_FMU_3PH Nonlinear synchronous-generator EMT benchmark.
% create_sg_grid_fmu_3ph([], false) builds and simulates without FMU export.
% create_sg_grid_fmu_3ph('/absolute/fresh/directory', true) also exports the
% ten-state controlled generator as an FMI 3.0 Model Exchange FMU.
% Requires MATLAB R2025a+, Simulink; no Simscape/Optimization Toolbox.
% FMU export additionally requires the standalone exporter and C compiler.
if nargin<2, exportFMU=false; end
if nargin<1 || isempty(outputDir)
    outputDir=fullfile('/Volumes/SSDMatthias/development/dpsim/build-macos-check',['sg_3ph_fmu_test_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
    stem=outputDir; k=1;
    while isfolder(outputDir), outputDir=sprintf('%s_%d',stem,k); k=k+1; end
end
sourceDir=fileparts(mfilename('fullpath'));
kernelDir=fullfile(sourceDir,'nonlinear_3ph');
addpath(kernelDir);
assert(~verLessThan('matlab','25.1'),'Requires MATLAB R2025a or later.');
if exportFMU
    assert(exist('exportToFMU','file')~=0,'Standalone FMU exporter is unavailable.');
end
if ~isfolder(outputDir), mkdir(outputDir); end
[ok,info]=fileattrib(outputDir); assert(ok,'Cannot access output directory.');
outputDir=info.Name;
oldDir=pwd; restore=onCleanup(@() cd(oldDir)); %#ok<NASGU>
cd(outputDir);
plant='dpsim_sg_3ph_me'; harness='dpsim_sg_3ph_reference';
assert(~bdIsLoaded(plant) && ~bdIsLoaded(harness),'Close existing SG benchmark models before rerunning.');
assert(~isfile([plant '.slx']) && ~isfile([harness '.slx']), ...
    'Benchmark models already exist. Choose a fresh output directory.');
b=sg_benchmark_parameters;
save(fullfile(outputDir,'sg_3ph_parameters.mat'),'b');
load_system('simulink');
new_system(plant); configure(plant,b);
% Keep machine, AVR and governor visibly separate. All dynamic states live in
% continuous Integrator blocks, so FMI exports actual continuous states.
add_block('simulink/Signal Routing/Mux',[plant '/Terminal_voltage_abc'], ...
    'Inputs','3','Position',[105 45 110 140]);
inputNames={'u_a','u_b','u_c','vref_pu','pmref_pu'};
for k=1:5
    y=25+45*k;
    add_block('simulink/Sources/In1',[plant '/' inputNames{k}], ...
        'Port',num2str(k),'PortDimensions','1','OutDataTypeStr','double', ...
        'SampleTime','0','Position',[25 y 55 y+18]);
    if k<=3, wire(plant,[inputNames{k} '/1'],['Terminal_voltage_abc/' num2str(k)]); end
end
addSubsystem(plant,'Machine',kernelDir,'sg_machine_rhs',b.machine,b.xMachine0, ...
    {'v_abc','efd_pu','pm_pu'},{'i_gen','diagnostic'},[650 40 845 190]);
addSubsystem(plant,'AVR',kernelDir,'sg_avr_rhs',b.avr,b.xAVR0, ...
    {'v_abc','vref_pu'},{'efd_pu','diagnostic'},[220 30 440 140]);
addSubsystem(plant,'Governor',kernelDir,'sg_governor_rhs',b.governor,b.pm0, ...
    {'speed_pu','pmref_pu'},{'diagnostic'},[220 230 440 350]);
% Subsystem ports: state first, followed by kernel's non-derivative outputs.
wire(plant,'Terminal_voltage_abc/1','Machine/1');
wire(plant,'Terminal_voltage_abc/1','AVR/1');
wire(plant,'vref_pu/1','AVR/2');
wire(plant,'AVR/2','Machine/2');
wire(plant,'Governor/1','Machine/3');
add_block('simulink/Signal Routing/Demux',[plant '/Machine_states'], ...
    'Outputs','6','Position',[920 20 925 210]);
wire(plant,'Machine/1','Machine_states/1');
wire(plant,'Machine_states/6','Governor/1');
wire(plant,'pmref_pu/1','Governor/2');
add_block('simulink/Signal Routing/Mux',[plant '/Record_vector'], ...
    'Inputs','7','Position',[1030 70 1035 430]);
% output vector = current3, machine diagnostics5, AVR diagnostics3,
% governor state1, governor diagnostics2, AVR state3, machine state6 = 23.
recordSources={'Machine/2','Machine/3','AVR/3','Governor/1', ...
    'Governor/2','AVR/1','Machine/1'};
for k=1:numel(recordSources), wire(plant,recordSources{k},['Record_vector/' num2str(k)]); end
outputNames={'i_a','i_b','i_c','speed_pu','theta_e_rad','field_current_A', ...
    'p_gap_W','p_gen_W','v_measured_pu','avr_raw_pu','avr_limited_pu', ...
    'pm_pu','governor_raw_pu','governor_limited_pu','v2_filter_pu', ...
    'avr_integral_pu','efd_pu','psi_d_Wb','psi_q_Wb','psi_0_Wb', ...
    'psi_f_Wb','theta_state_rad','speed_state_pu'};
add_block('simulink/Signal Routing/Demux',[plant '/Named_outputs'], ...
    'Outputs',num2str(numel(outputNames)),'Position',[1130 25 1135 1170]);
wire(plant,'Record_vector/1','Named_outputs/1');
for k=1:numel(outputNames)
    y=20+48*(k-1);
    add_block('simulink/Sinks/Out1',[plant '/' outputNames{k}], ...
        'Port',num2str(k),'Position',[1200 y 1230 y+18]);
    wire(plant,['Named_outputs/' num2str(k)],[outputNames{k} '/1']);
end
save_system(plant,fullfile(outputDir,[plant '.slx']));

new_system(harness); configure(harness,b);
% A copied subsystem retains ordinary blocks, no external model dependency.
add_block('built-in/Subsystem',[harness '/Controlled_generator'], ...
    'Position',[460 80 665 450]);
Simulink.BlockDiagram.copyContentsToSubsystem(plant,[harness '/Controlled_generator']);
addSubsystem(harness,'Network',kernelDir,'sg_network_rhs',b.network,b.xNetwork0, ...
    {'i_gen','v_grid','g_load'},{'diagnostic'},[760 80 960 240]);
add_block('simulink/Signal Routing/Demux',[harness '/Network_states'], ...
    'Outputs','[3 3]','Position',[1030 45 1035 180]);
wire(harness,'Network/1','Network_states/1');
add_block('simulink/Signal Routing/Demux',[harness '/Bus_voltages'], ...
    'Outputs','3','Position',[310 35 315 150]);
wire(harness,'Network_states/1','Bus_voltages/1');
add_block('simulink/Signal Routing/Mux',[harness '/Generator_current'], ...
    'Inputs','3','Position',[705 50 710 160]);
for k=1:3
    wire(harness,['Bus_voltages/' num2str(k)],['Controlled_generator/' num2str(k)]);
    wire(harness,['Controlled_generator/' num2str(k)],['Generator_current/' num2str(k)]);
end
wire(harness,'Generator_current/1','Network/1');
add_block('simulink/Signal Routing/Mux',[harness '/Grid_voltages'], ...
    'Inputs','3','Position',[180 530 185 670]);
for k=1:3
    name=['Grid_' char('a'+k-1)]; y=510+60*k;
    add_block('simulink/Sources/Sine Wave',[harness '/' name], ...
        'Amplitude',num2str(sqrt(2)*b.Vphase,17),'Frequency',num2str(b.wb,17), ...
        'Phase',num2str(b.phase(k)+pi/2,17),'SampleTime','0', ...
        'Position',[35 y 105 y+35]);
    wire(harness,[name '/1'],['Grid_voltages/' num2str(k)]);
end
wire(harness,'Grid_voltages/1','Network/2');
add_block('simulink/Sources/Step',[harness '/Load_step'], ...
    'Time',num2str(b.loadStepTime,17),'Before',num2str(b.G0,17), ...
    'After',num2str(b.G0+b.Gstep,17),'SampleTime','0', ...
    'Position',[250 570 330 610]);
wire(harness,'Load_step/1','Network/3');
add_block('simulink/Sources/Step',[harness '/Voltage_reference_step'], ...
    'Time',num2str(b.vrefStepTime,17),'Before','1', ...
    'After',num2str(1+b.vrefStep,17),'SampleTime','0','Position',[250 250 340 290]);
wire(harness,'Voltage_reference_step/1','Controlled_generator/4');
add_block('simulink/Sources/Constant',[harness '/Mechanical_reference'], ...
    'Value',num2str(b.pm0,17),'Position',[250 370 340 410]);
wire(harness,'Mechanical_reference/1','Controlled_generator/5');
add_block('simulink/Signal Routing/Mux',[harness '/Generator_outputs'], ...
    'Inputs',num2str(numel(outputNames)),'Position',[1090 250 1095 800]);
for k=1:numel(outputNames)
    wire(harness,['Controlled_generator/' num2str(k)],['Generator_outputs/' num2str(k)]);
end
add_block('simulink/Signal Routing/Mux',[harness '/All_measurements'], ...
    'Inputs','5','Position',[1200 70 1205 310]);
records={'Generator_outputs/1','Network/1','Load_step/1', ...
    'Voltage_reference_step/1','Network/2'};
for k=1:numel(records), wire(harness,records{k},['All_measurements/' num2str(k)]); end
add_block('simulink/Sinks/To Workspace',[harness '/Record'], ...
    'VariableName','sg_reference','SaveFormat','Timeseries', ...
    'MaxDataPoints','inf','Decimation','1','Position',[1280 120 1390 160]);
wire(harness,'All_measurements/1','Record/1');
% Save the model before simulation so failures leave an inspectable artifact.
save_system(harness,fullfile(outputDir,[harness '.slx']));
fprintf('Built nonlinear grid. Simulating 16 continuous states for %.1f s...\n',b.stopTime);
result=sim(harness,'ReturnWorkspaceOutputs','on');
signal=result.get('sg_reference');
columns=[outputNames,{'v_bus_a_V','v_bus_b_V','v_bus_c_V', ...
    'i_line_a_A','i_line_b_A','i_line_c_A','g_load_S','vref_pu', ...
    'p_load_W','p_line_bus_W','p_line_loss_W'}];
time=signal.Time(:);
% Matrix-valued Simulink signals store time in the LAST dimension, whereas
% vector-valued timeseries commonly store it first. Never reshape blindly.
if signal.IsTimeFirst
    data=reshape(signal.Data,numel(time),[]);
else
    data=reshape(signal.Data,[],numel(time)).';
end
assert(size(data,2)==numel(columns),'Unexpected recorded signal dimensions.');
% Remove solver event duplicate timestamps by retaining the right-hand sample.
[time,index]=unique(time,'last'); data=data(index,:);
reference=array2table([time,data],'VariableNames',[{'time_s'},columns]);
writetable(reference,fullfile(outputDir,'reference_simulink.csv'));
% Store a uniform sample grid, preserving precisely scheduled input steps.
t=(0:round(b.stopTime/b.sampleStep)).'*b.sampleStep;
uniform=interp1(time,data,t,'linear');
uniform(:,strcmp(columns,'g_load_S'))=b.G0+b.Gstep*(t>=b.loadStepTime);
uniform(:,strcmp(columns,'vref_pu'))=1+b.vrefStep*(t>=b.vrefStepTime);
writetable(array2table([t,uniform],'VariableNames',[{'time_s'},columns]), ...
    fullfile(outputDir,'reference_uniform.csv'));
% Independent wiring/integration check: same equations in a MATLAB ODE solver,
% with exact restarts at the two known discontinuity times.
referenceODE=sg_reference_ode(b,t);
writetable(array2table([t,referenceODE],'VariableNames',[{'time_s'},columns]), ...
    fullfile(outputDir,'reference_matlab_ode15s.csv'));
metrics=sg_validate_reference(b,t,uniform,referenceODE,outputDir);
save(fullfile(outputDir,'sg_3ph_reference.mat'),'b','t','uniform','referenceODE','columns','metrics');
mapping=struct('fmi_version','3.0','interface','ModelExchange','domain','EMT', ...
    'frame','abc','voltage_inputs',{{'u_a','u_b','u_c'}}, ...
    'current_outputs',{{'i_a','i_b','i_c'}},'current_sign','out_of_generator_into_bus', ...
    'control_inputs',{{'vref_pu','pmref_pu'}}, ...
    'initial_voltages_V',b.xNetwork0(1:3).','initial_vref_pu',1, ...
    'initial_pmref_pu',b.pm0,'generator_states',10,'network_states',6, ...
    'load_step_time_s',b.loadStepTime,'load_conductance_step_S',b.Gstep, ...
    'vref_step_time_s',b.vrefStepTime,'vref_step_pu',b.vrefStep, ...
    'sample_step_s',b.sampleStep,'stop_time_s',b.stopTime, ...
    'requires_nonlinear_importer',true);
writeText(fullfile(outputDir,'port_mapping.json'),jsonencode(mapping));
sg_plot_reference(b,t,uniform,columns,outputDir);
if exportFMU
    exportToFMU(plant,'FMIVersion','3.0','FMUType','ME','SaveDirectory',outputDir);
end
fprintf('Completed nonlinear three-phase benchmark:\n%s\n',outputDir);
fprintf('Initial Pgen = %.3f kW, Pgrid = %.3f kW, field = %.3f A\n', ...
    b.initialPgen/1000,b.initialPgrid/1000,b.initialFieldCurrent);
fprintf('This nonlinear FMU requires a future nonlinear DPsim adapter; LinearSSN must reject it.\n');
end

function configure(model,b)
set_param(model,'StartTime','0','StopTime',num2str(b.stopTime,17), ...
    'SolverType','Variable-step','Solver','ode15s', ...
    'MaxStep',num2str(b.sampleStep,17),'RelTol','1e-7','AbsTol','1e-9');
end

function addSubsystem(model,name,kernelDir,kernel,p,x0,inputs,outputs,position)
sys=[model '/' name];
add_block('built-in/Subsystem',sys,'Position',position);
add_block('simulink/Continuous/Integrator',[sys '/States'], ...
    'InitialCondition',mat2str(x0,17),'Position',[340 180 380 220]);
add_block('simulink/Sources/Constant',[sys '/Parameters_SI'], ...
    'Value',mat2str(p,17),'Position',[40 280 150 320]);
add_block('simulink/User-Defined Functions/MATLAB Function',[sys '/Equations'], ...
    'Position',[205 50 365 140]);
chart=find(sfroot,'-isa','Stateflow.EMChart','Path',[sys '/Equations']);
chart.Script=fileread(fullfile(kernelDir,[kernel '.m']));
wire(sys,'States/1','Equations/1');
for k=1:numel(inputs)
    y=30+55*k;
    add_block('simulink/Sources/In1',[sys '/' inputs{k}], ...
        'Port',num2str(k),'Position',[30 y 60 y+18]);
    wire(sys,[inputs{k} '/1'],['Equations/' num2str(k+1)]);
end
wire(sys,'Parameters_SI/1',['Equations/' num2str(numel(inputs)+2)]);
wire(sys,'Equations/1','States/1');
add_block('simulink/Sinks/Out1',[sys '/state'], ...
    'Port','1','Position',[455 195 485 215]);
wire(sys,'States/1','state/1');
for k=1:numel(outputs)
    y=35+55*k;
    add_block('simulink/Sinks/Out1',[sys '/' outputs{k}], ...
        'Port',num2str(k+1),'Position',[455 y 485 y+18]);
    wire(sys,['Equations/' num2str(k+1)],[outputs{k} '/1']);
end
end

function wire(model,source,destination)
add_line(model,source,destination,'autorouting','on');
end

function writeText(path,text)
f=fopen(path,'w'); assert(f~=-1,'Cannot write output file.');
c=onCleanup(@() fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',text);
end
