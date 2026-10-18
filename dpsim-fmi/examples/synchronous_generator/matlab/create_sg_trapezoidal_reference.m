function create_sg_trapezoidal_reference(outputDir, steps)
%CREATE_SG_TRAPEZOIDAL_REFERENCE Full-state implicit trapezoidal + Newton.
% All outputs go to the supplied directory. Scheduled inputs are held per
% interval, with exact endpoints at 1/2 s; no half-step anticipation of events.
if nargin<2, steps=[50e-6,25e-6]; end
source=fileparts(mfilename('fullpath'));
addpath(fullfile(source,'nonlinear_3ph'));
if ~isfolder(outputDir), mkdir(outputDir); end
b=sg_benchmark_parameters;
columns={'i_a','i_b','i_c','speed_pu','theta_e_rad','field_current_A', ...
'p_gap_W','p_gen_W','v_measured_pu','avr_raw_pu','avr_limited_pu', ...
'pm_pu','governor_raw_pu','governor_limited_pu','v2_filter_pu', ...
'avr_integral_pu','efd_pu','psi_d_Wb','psi_q_Wb','psi_0_Wb', ...
'psi_f_Wb','theta_state_rad','speed_state_pu','v_bus_a_V','v_bus_b_V', ...
'v_bus_c_V','i_line_a_A','i_line_b_A','i_line_c_A','g_load_S','vref_pu', ...
'p_load_W','p_line_bus_W','p_line_loss_W'};
for h=steps
    assert(all(abs([b.loadStepTime,b.vrefStepTime,b.stopTime]/h-round( ...
        [b.loadStepTime,b.vrefStepTime,b.stopTime]/h))<1e-8),'Events must align to h.');
    t=(0:round(b.stopTime/h)).'*h;
    states=zeros(numel(t),16); states(1,:)=b.x0.';
    iterations=zeros(numel(t),1); residuals=zeros(numel(t),1);
    tic;
    for k=2:numel(t)
        midpoint=(t(k-1)+t(k))/2;
        g=b.G0+b.Gstep*(midpoint>=b.loadStepTime);
        vref=1+b.vrefStep*(midpoint>=b.vrefStepTime);
        old=states(k-1,:).';
        f0=sg_full_rhs_reference(t(k-1),old,b,g,vref);
        trial=old+h*f0;
        converged=false;
        for it=1:15
            f=sg_full_rhs_reference(t(k),trial,b,g,vref);
            r=trial-old-h*(f0+f)/2;
            scale=1e-11+1e-10*max([ones(16,1),abs(old),abs(trial)],[],2);
            normResidual=max(abs(r)./scale);
            if normResidual<=1, converged=true; break; end
            J=zeros(16);
            for j=1:16
                point=complex(trial); point(j)=point(j)+1i*1e-20;
                J(:,j)=imag(sg_full_rhs_reference(t(k),point,b,g,vref))/1e-20;
            end
            correction=(eye(16)-h*J/2)\r;
            trial=trial-correction;
        end
        assert(converged && all(isfinite(trial)),'Trapezoidal Newton failed at t=%.9g',t(k));
        states(k,:)=trial.'; iterations(k)=it; residuals(k)=normResidual;
        if mod(k-1,round(.5/h))==0
            fprintf('Trapezoidal h=%.9g, t=%.2f s, elapsed %.1f s\n',h,t(k),toc);
        end
    end
    y=zeros(numel(t),34);
    for k=1:numel(t)
        x=states(k,:).';
        g=b.G0+b.Gstep*(t(k)>=b.loadStepTime-1e-12);
        vref=1+b.vrefStep*(t(k)>=b.vrefStepTime-1e-12);
        [~,igen,dm]=sg_machine_rhs(x(1:6),x(11:13),x(9),x(10),b.machine);
        [~,~,da]=sg_avr_rhs(x(7:9),x(11:13),vref,b.avr);
        [~,dg]=sg_governor_rhs(x(10),x(6),b.pm0,b.governor);
        grid=sqrt(2)*b.Vphase*cos(b.wb*t(k)+b.phase);
        [~,dn]=sg_network_rhs(x(11:16),igen,grid,g,b.network);
        y(k,:)=[igen;dm;da;x(10);dg;x(7:9);x(1:6);x(11:16);g;vref;dn].';
    end
    suffix=sprintf('%gus',h*1e6);
    writetable(array2table([t,y,iterations,residuals], ...
        'VariableNames',[{'time_s'},columns,{'newton_iterations','scaled_residual'}]), ...
        fullfile(outputDir,['reference_trapezoidal_' suffix '.csv']));
    save(fullfile(outputDir,['reference_trapezoidal_' suffix '.mat']), ...
        'b','h','t','states','y','iterations','residuals','columns');
    fprintf('Completed h=%.9g; maximum Newton iterations %d\n',h,max(iterations));
end
% Parameters/ports for the C++ network example, without machine equations.
config=struct('frequency_Hz',b.frequency,'phase_rms_V',b.Vphase, ...
    'stop_time_s',b.stopTime,'step_s',steps(1),'network',b.network, ...
    'initial_network_state',b.xNetwork0.','initial_generator_current', ...
    (b.xNetwork0(4:6)+b.G0*b.xNetwork0(1:3)+ ...
    b.network(1)*(-b.wb*sqrt(2)*b.Vphase*sin(b.beta+b.phase))).', ...
    'G0_S',b.G0,'Gstep_S',b.Gstep,'load_step_s',b.loadStepTime, ...
    'vref_step_s',b.vrefStepTime,'vref_step_pu',b.vrefStep,'pmref_pu',b.pm0, ...
    'voltage_inputs',{{'u_a','u_b','u_c'}},'current_outputs',{{'i_a','i_b','i_c'}}, ...
    'control_names',{{'vref_pu','pmref_pu'}},'control_initial_values',[1,b.pm0], ...
    'current_sign_to_passive',-1,'power_base_VA',b.Sbase);
config.bases=struct('phase_voltage_V',sqrt(2)*b.Vphase, ...
    'phase_current_A',sqrt(2)*b.Sbase/(3*b.Vphase),'power_W',b.Sbase, ...
    'field_current_A',b.machine(12)/b.machine(7), ...
    'stator_flux_Wb',sqrt(3)*b.Vphase/b.wb, ...
    'field_flux_Wb',b.machine(5)*b.machine(12)/b.machine(7), ...
    'angle_rad',1,'per_unit',1,'conductance_S',b.G0);
f=fopen(fullfile(outputDir,'benchmark.json'),'w'); assert(f~=-1);
c=onCleanup(@() fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(config));
end
