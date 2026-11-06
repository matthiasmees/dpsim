function run_post_powerflow_simulink_sweep(repoRoot, outRoot, fmuPath, operatingPointPath, fullSweep)
% Run fresh native and FMI-CS Simulink closed-loop cases after an independent
% one-bus AC power-flow initialization. Outputs are saved case by case.
if nargin < 4, error('repoRoot, outRoot, fmuPath and operatingPointPath are required.'); end
if nargin < 5, fullSweep=true; end
addpath(fullfile(repoRoot,'dpsim-fmi','examples','small_inverter_grid','matlab'));
if ~isfolder(outRoot),mkdir(outRoot);end
op=jsondecode(fileread(operatingPointPath));
stepsUs=[25 50 100 500 1000];
stopTime=5; loadStepTime=2.5;
writeJson(fullfile(outRoot,'operating_point.json'),op);

% Short fresh checks run first, before the full sweep.
for modeCell={'native','cs'}
    runOne(modeCell{1},25e-6,0.05,0.025,'smoke',op,repoRoot,outRoot,fmuPath);
end
if ~fullSweep
    fprintf('Both 50 ms post-PF Simulink smoke runs completed.\n');
    return;
end

for modeCell={'native','cs'}
    mode=modeCell{1};
    for dtUs=stepsUs
        runOne(mode,dtUs*1e-6,stopTime,loadStepTime,sprintf('dt_%dus',dtUs), ...
            op,repoRoot,outRoot,fmuPath);
    end
end
writeJson(fullfile(outRoot,'complete.json'),struct('complete',true, ...
    'steps_us',stepsUs,'stop_time_s',stopTime,'load_step_s',loadStepTime, ...
    'modes',{{'native','cs'}},'smoke_stop_s',0.05));
fprintf('Completed fresh Simulink native and CS post-PF sweep under %s\n',outRoot);
end

function runOne(mode,dt,stopTime,loadStepTime,caseName,op,repoRoot,outRoot,fmuPath)
caseDir=fullfile(outRoot,caseName,mode);
if ~isfolder(caseDir),mkdir(caseDir);end
modelDirExisting=fullfile(repoRoot,'build-fmi','small-inverter-grid','models');
sourceModel=fullfile(modelDirExisting,sprintf('SmallGrid_%s.slx',mode));
if ~isfile(sourceModel),error('Existing closed-loop Simulink model is missing: %s',sourceModel);end
load_system(sourceModel); [~,mdl]=fileparts(sourceModel);
assignin('base','dt',dt); assignin('base','stopTime',stopTime);
assignin('base','loadStepTime',loadStepTime);
set_param([mdl '/PI line series current'],'InitialCondition',mat2str(op.i_line_abc_A(:),16));
set_param([mdl '/PCC voltage'],'InitialCondition',mat2str(op.v_pcc_abc_V(:),16));
if strcmp(mode,'native')
    set_param([mdl '/Inverter passive current state'],'InitialCondition', ...
        mat2str(op.i_inverter_abc_A(:),16));
    chart=sfroot().find('-isa','Stateflow.EMChart','Path', ...
        [mdl '/Closed-loop PI line and inverter equations']);
    chart.Script=nativeEquations(op.source_phase_peak_V,op.source_phase_rad);
else
    chart=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/Closed-loop PI line equations']);
    chart.Script=networkEquations(op.source_phase_peak_V,op.source_phase_rad);
    addpath(fileparts(fmuPath));
end
set_param(mdl,'Solver','ode14x','SolverType','Fixed-step', ...
    'NumberNewtonIterations','2','FixedStep',num2str(dt,'%.17g'), ...
    'StopTime',num2str(stopTime,'%.17g'),'SaveTime','off','SaveOutput','off', ...
    'ReturnWorkspaceOutputs','on','SignalLogging','off');
set_param(mdl,'SimulationCommand','update');
input=Simulink.SimulationInput(mdl);
input=input.setVariable('dt',dt).setVariable('stopTime',stopTime) ...
    .setVariable('loadStepTime',loadStepTime);
input=input.setModelParameter('FixedStep',num2str(dt,'%.17g'), ...
    'StopTime',num2str(stopTime,'%.17g'),'SimulationMode','normal');
timer=tic; result=sim(input); wall=toc(timer);
voltage=result.get('pcc_voltage'); current=result.get('inverter_current');
if isempty(voltage)||isempty(current),error('Missing voltage/current output: %s %s',mode,caseName);end
t=voltage.Time(:); vabc=squeeze(voltage.Data); iabc=squeeze(current.Data);
if size(vabc,1)~=numel(t)&&size(vabc,2)==numel(t),vabc=vabc.';end
if size(iabc,1)~=numel(t)&&size(iabc,2)==numel(t),iabc=iabc.';end
if size(vabc,2)~=3||size(iabc,2)~=3,error('Bad signal dimensions in %s %s.',mode,caseName);end
if any(~isfinite(vabc),'all')||any(~isfinite(iabc),'all'),error('Non-finite signal in %s %s.',mode,caseName);end
expected=round(stopTime/dt)+1;
if numel(t)~=expected||abs(t(1))>1e-12||abs(t(end)-stopTime)>dt*1e-6|| ...
        any(abs(diff(t)-dt)>max(1e-11,dt*1e-7))
    error('Unexpected fixed-step time grid for %s %s.',mode,caseName);
end
trace=table(t,vabc(:,1),vabc(:,2),vabc(:,3),iabc(:,1),iabc(:,2),iabc(:,3), ...
    'VariableNames',{'time_s','v_a_V','v_b_V','v_c_V','i_a_A','i_b_A','i_c_A'});
tracePath=fullfile(caseDir,sprintf('simulink_%s.csv',mode));writetable(trace,tracePath);
timing=struct('mode',mode,'case',caseName,'step_s',dt,'stop_time_s',stopTime, ...
    'load_step_s',loadStepTime,'samples',numel(t),'fresh_run',true,'closed_loop',true, ...
    'network','Independent analytical AC PF initial state; Simulink ode14x closed-loop PI-line network', ...
    'inverter',ternary(strcmp(mode,'native'),'native averaged inverter with periodic current state', ...
    'FMI 3.0 CS FMU with matching independent AC PF starts'), ...
    'outer_solver','ode14x, two Newton iterations','simulink_wall_time_s',wall, ...
    'trace',tracePath);
writeJson(fullfile(caseDir,sprintf('simulink_%s_timing.json',mode)),timing);
fprintf('%-6s %-9s %g us: %d samples, %.3f s wall\n',mode,caseName,dt*1e6,numel(t),wall);
clear result voltage current t vabc iabc trace input;
close_system(mdl,0);
end

function text=nativeEquations(sourceVpk,sourcePhase)
text=sprintf(['function [diLine,dvPcc,diInv] = fcn(t,iLine,vPcc,iInv,loadStepTime)\n' ...
    '%%#codegen\n' ...
    'R=0.30; L=2.0e-2; C=5.0e-7; f=60; w=2*pi*f; tau=5e-3;\n' ...
    'Vpk=%.17g; phi=[%.17g;%.17g;%.17g]; vsrc=Vpk*cos(w*t+phi);\n' ...
    'P=125e6; Q=80.881e6; Vnom=345e3; Ipu=1.20;\n' ...
    'alpha=(2/3)*(vPcc(1)-0.5*vPcc(2)-0.5*vPcc(3)); beta=(2/3)*(sqrt(3)/2)*(vPcc(2)-vPcc(3)); v2=alpha*alpha+beta*beta;\n' ...
    'if v2<1, iRef=zeros(3,1); else\n' ...
    'k=2/(3*v2); ia=k*(P*alpha+Q*beta); ib=k*(P*beta-Q*alpha); imax=Ipu*sqrt(2)*hypot(P,Q)/(sqrt(3)*Vnom); imag=hypot(ia,ib);\n' ...
    'if imag>imax, ia=ia*imax/imag; ib=ib*imax/imag; end\n' ...
    'iRef=[-ia;0.5*ia-sqrt(3)/2*ib;0.5*ia+sqrt(3)/2*ib]; end\n' ...
    'iLoad=zeros(3,1); if t>=loadStepTime, iLoad=vPcc/(362215.5^2/(50e6)); end\n' ...
    'diLine=(vsrc-vPcc-R*iLine)/L; dvPcc=(iLine-iInv-iLoad)/C; diInv=(iRef-iInv)/tau;\n' ...
    'end\n'],sourceVpk,sourcePhase(1),sourcePhase(2),sourcePhase(3));
end

function text=networkEquations(sourceVpk,sourcePhase)
text=sprintf(['function [diLine,dvPcc] = fcn(t,iLine,vPcc,iInv,loadStepTime)\n' ...
    '%%#codegen\n' ...
    'R=0.30; L=2.0e-2; C=5.0e-7; f=60; w=2*pi*f;\n' ...
    'Vpk=%.17g; phi=[%.17g;%.17g;%.17g]; vsrc=Vpk*cos(w*t+phi);\n' ...
    'iLoad=zeros(3,1); if t>=loadStepTime, iLoad=vPcc/(362215.5^2/(50e6)); end\n' ...
    'diLine=(vsrc-vPcc-R*iLine)/L; dvPcc=(iLine-iInv-iLoad)/C;\n' ...
    'end\n'],sourceVpk,sourcePhase(1),sourcePhase(2),sourcePhase(3));
end

function writeJson(path,value)
fid=fopen(path,'w');if fid<0,error('Cannot write %s',path);end
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(value,'PrettyPrint',true));fwrite(fid,newline);clear cleanup;
end

function y=ternary(test,a,b)
if test,y=a;else,y=b;end
end
