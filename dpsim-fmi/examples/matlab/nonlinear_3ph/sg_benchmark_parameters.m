function b = sg_benchmark_parameters
% Example SI parameters, deliberately not a fit to a commercial generator.
b.frequency=50; b.wb=2*pi*b.frequency; b.Vphase=230; b.Sbase=15000;
b.phase=[0; -2*pi/3; 2*pi/3];
b.machine=[.2, .060, .045, .010, 1.5, .27, 1.5, b.wb, b.Sbase, 1.5, 2, 0];
b.machine(12)=b.machine(7)*sqrt(3)*b.Vphase/(b.wb*b.machine(6));
b.avr=[b.Vphase, .02, 6, 12, 10, 0, 5, .05];
b.governor=[.05, .2, 0, 1.2];
b.network=[40e-6, .3, 5e-3];
b.Pload0=6000; b.PloadStep=6000; b.Pgen0=8000;
b.loadStepTime=1; b.vrefStepTime=2; b.vrefStep=.02;
b.stopTime=4; b.sampleStep=50e-6;
b.G0=b.Pload0/(3*b.Vphase^2);
b.Gstep=b.PloadStep/(3*b.Vphase^2);
% Solve a balanced sinusoidal operating point without Optimization Toolbox.
zline=b.network(2)+1i*b.wb*b.network(3);
powerAt=@(beta) 3*real((b.Vphase*exp(1i*beta))*conj( ...
    (b.G0+1i*b.wb*b.network(1))*b.Vphase*exp(1i*beta)+ ...
    (b.Vphase*exp(1i*beta)-b.Vphase)/zline));
b.beta=fzero(@(beta) powerAt(beta)-b.Pgen0,[-.5,.5]);
vphasor=b.Vphase*exp(1i*b.beta);
ilinephasor=(vphasor-b.Vphase)/zline;
igenphasor=(b.G0+1i*b.wb*b.network(1))*vphasor+ilinephasor;
vabc=sqrt(2)*real(vphasor*exp(1i*b.phase));
iabc=-sqrt(2)*real(igenphasor*exp(1i*b.phase));
b.delta=fzero(@(delta) dResidual(delta,vabc,iabc,b),[-1,1.3]);
theta=b.delta-pi/2;
a=theta+b.phase;
T=sqrt(2/3)*[cos(a).'; -sin(a).'; ones(1,3)/sqrt(2)];
v=T*vabc; i=T*iabc;
Rs=b.machine(1); Ld=b.machine(2); Lq=b.machine(3); Lf=b.machine(5); M=b.machine(6);
ifield=(v(2)-Rs*i(2)-b.wb*Ld*i(1))/(b.wb*M);
psi_d=Ld*i(1)+M*ifield; psi_q=Lq*i(2); psi_f=M*i(1)+Lf*ifield;
efd0=b.machine(7)*ifield/b.machine(12);
pm0=-b.wb*(psi_d*i(2)-psi_q*i(1))/b.Sbase;
assert(efd0>b.avr(6) && efd0<b.avr(7) && ifield>0,'Invalid field initialization.');
assert(pm0>b.governor(3) && pm0<b.governor(4),'Invalid prime-mover initialization.');
b.xMachine0=[psi_d; psi_q; 0; psi_f; theta; 1];
b.xAVR0=[1; efd0; efd0]; b.pm0=pm0;
b.xNetwork0=[vabc; sqrt(2)*real(ilinephasor*exp(1i*b.phase))];
b.x0=[b.xMachine0; b.xAVR0; b.pm0; b.xNetwork0];
[dx,igen,~]=sg_machine_rhs(b.xMachine0,vabc,efd0,pm0,b.machine);
assert(norm(dx([1:4,6]))<1e-8,'Machine is not initialized at equilibrium.');
assert(norm(igen+iabc)<1e-9,'Park-transform current sign is inconsistent.');
b.initialPgen= vabc.'*igen;
b.initialPgrid=3*real(vphasor*conj(ilinephasor));
b.initialFieldCurrent=ifield;
end

function r=dResidual(delta,vabc,iabc,b)
a=delta-pi/2+b.phase;
T=sqrt(2/3)*[cos(a).';-sin(a).';ones(1,3)/sqrt(2)];
v=T*vabc; i=T*iabc;
r=v(1)-b.machine(1)*i(1)+b.wb*b.machine(3)*i(2);
end
