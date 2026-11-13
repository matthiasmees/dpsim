function metrics=sg_validate_reference(b,t,y,reference,outputDir)
% Signal order is defined alongside the model's named ports in the builder.
assert(all(isfinite(y),'all') && all(isfinite(reference),'all'),'Nonfinite reference.');
error=y-reference;
metrics.max_current_error_A=max(abs(error(:,1:3)),[],'all');
metrics.max_bus_voltage_error_V=max(abs(error(:,24:26)),[],'all');
metrics.max_speed_error_pu=max(abs(error(:,4)));
metrics.max_efd_error_pu=max(abs(error(:,17)));
% Interpolation of the adaptive Simulink output is included in these bounds.
assert(metrics.max_current_error_A<.03,'Simulink/MATLAB current mismatch.');
assert(metrics.max_bus_voltage_error_V<.1,'Simulink/MATLAB voltage mismatch.');
assert(metrics.max_speed_error_pu<2e-5,'Simulink/MATLAB speed mismatch.');
assert(metrics.max_efd_error_pu<2e-4,'Simulink/MATLAB exciter mismatch.');
metrics.max_current_sum_A=max(abs(sum(y(:,1:3),2)));
metrics.max_voltage_sum_V=max(abs(sum(y(:,24:26),2)));
assert(metrics.max_current_sum_A<1e-6 && metrics.max_voltage_sum_V<1e-6, ...
    'Balanced abc wiring check failed.');
metrics.speed_min_pu=min(y(:,4)); metrics.speed_max_pu=max(y(:,4));
assert(metrics.speed_min_pu>.9 && metrics.speed_max_pu<1.1,'Machine lost synchronous operating range.');
assert(all(y(:,17)>=-1e-7 & y(:,17)<=b.avr(7)+1e-7),'Excitation escaped limits.');
assert(all(y(:,12)>=-1e-7 & y(:,12)<=b.governor(4)+1e-7),'Prime mover escaped limits.');
% Before disturbances the initialized machine stays on its periodic orbit.
pre=t>.1 & t<.9;
metrics.predisturbance_speed_error_pu=max(abs(y(pre,4)-1));
metrics.predisturbance_voltage_error_pu=max(abs(y(pre,9)-1));
assert(metrics.predisturbance_speed_error_pu<1e-6,'Initialization created a speed transient.');
assert(metrics.predisturbance_voltage_error_pu<1e-6,'Initialization created a voltage transient.');
metrics.peak_speed_deviation_Hz=max(abs((y(:,4)-1)*b.frequency));
metrics.voltage_min_pu=min(y(:,9)); metrics.voltage_max_pu=max(y(:,9));
metrics.final_voltage_pu=y(end,9); metrics.final_pm_pu=y(end,12);
metrics.max_avr_clipping_pu=max(abs(y(:,10)-y(:,11)));
metrics.max_governor_clipping_pu=max(abs(y(:,13)-y(:,14)));
% Check full instantaneous energy conservation at off-equilibrium samples.
% Magnetic energy uses the symmetric power-invariant flux matrix; kinetic
% energy H*Sbase*speed^2. This detects Park scaling and torque sign errors.
residual=zeros(51,1);
indices=round(linspace(1,numel(t),numel(residual)));
for k=1:numel(indices)
    j=indices(k); row=y(j,:);
    xm=row(18:23).'; xa=row(15:17).'; pm=row(12);
    xn=row(24:29).'; v=xn(1:3);
    [dm,igen,d]=sg_machine_rhs(xm,v,xa(3),pm,b.machine);
    grid=sqrt(2)*b.Vphase*cos(b.wb*t(j)+b.phase);
    [dn,nd]=sg_network_rhs(xn,igen,grid,row(30),b.network);
    a=xm(5)+b.phase;
    T=sqrt(2/3)*[cos(a).';-sin(a).';ones(1,3)/sqrt(2)];
    passive=-T*igen;
    dWmag=passive.'*dm(1:3)+d(3)*dm(4);
    dWkin=2*b.machine(10)*b.Sbase*xm(6)*dm(6);
    dWnetwork=b.network(1)*v.'*dn(1:3)+b.network(3)*xn(4:6).'*dn(4:6);
    fieldPower=b.machine(12)*xa(3)*d(3);
    losses=b.machine(1)*(passive.'*passive)+b.machine(7)*d(3)^2+nd(3);
    gridPower=grid.'*xn(4:6);
    damping=b.machine(11)*(xm(6)-1)*b.Sbase;
    balance=pm*b.Sbase+fieldPower-nd(1)-gridPower-losses-damping;
    residual(k)=dWmag+dWkin+dWnetwork-balance;
end
metrics.max_energy_balance_residual_W=max(abs(residual));
assert(metrics.max_energy_balance_residual_W<1e-7,'Energy/torque/power convention error.');
file=fopen(fullfile(outputDir,'validation.json'),'w');
assert(file~=-1,'Cannot save validation report.');
cleanup=onCleanup(@() fclose(file)); %#ok<NASGU>
fprintf(file,'%s\n',jsonencode(metrics));
disp(metrics);
end
