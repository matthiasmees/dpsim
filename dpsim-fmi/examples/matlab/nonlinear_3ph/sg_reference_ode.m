function y = sg_reference_ode(b,t)
% Independently wire the same component kernels in MATLAB, no Simulink blocks.
% Controls/load are held on each interval; restart exactly at each step.
boundaries=[0,b.loadStepTime,b.vrefStepTime,b.stopTime];
x=b.x0; state=zeros(numel(t),16);
opts=odeset('RelTol',2e-9,'AbsTol',1e-10,'MaxStep',b.sampleStep/2);
for segment=1:numel(boundaries)-1
    start=boundaries(segment); stop=boundaries(segment+1);
    g=b.G0+b.Gstep*(start>=b.loadStepTime);
    vref=1+b.vrefStep*(start>=b.vrefStepTime);
    sol=ode15s(@(time,x) rhs(time,x,b,g,vref),[start,stop],x,opts);
    mask=t>=start & t<=stop;
    state(mask,:)=deval(sol,t(mask)).';
    x=sol.y(:,end);
end
y=zeros(numel(t),34);
for k=1:numel(t)
    x=state(k,:).';
    g=b.G0+b.Gstep*(t(k)>=b.loadStepTime);
    vref=1+b.vrefStep*(t(k)>=b.vrefStepTime);
    [~,igen,dm]=sg_machine_rhs(x(1:6),x(11:13),x(9),x(10),b.machine);
    [~,~,da]=sg_avr_rhs(x(7:9),x(11:13),vref,b.avr);
    [~,dg]=sg_governor_rhs(x(10),x(6),b.pm0,b.governor);
    grid=sqrt(2)*b.Vphase*cos(b.wb*t(k)+b.phase);
    [~,dn]=sg_network_rhs(x(11:16),igen,grid,g,b.network);
    y(k,:)=[igen;dm;da;x(10);dg;x(7:9);x(1:6);x(11:16);g;vref;dn].';
end
end

function dx=rhs(t,x,b,g,vref)
v=x(11:13);
[dm,igen]=sg_machine_rhs(x(1:6),v,x(9),x(10),b.machine);
da=sg_avr_rhs(x(7:9),v,vref,b.avr);
dg=sg_governor_rhs(x(10),x(6),b.pm0,b.governor);
grid=sqrt(2)*b.Vphase*cos(b.wb*t+b.phase);
dn=sg_network_rhs(x(11:16),igen,grid,g,b.network);
dx=[dm;da;dg;dn];
end
