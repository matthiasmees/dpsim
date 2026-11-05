function dx=sg_full_rhs_reference(t,x,b,g,vref)
% Independent combined equations for the sixteen-state trapezoidal reference.
% Unlike FMU/MNA condensation, Newton here solves one full state residual.
% Nonconjugating products and real-branch tests permit complex-step Jacobians.
p=b.machine; theta=x(5); omega=p(8)*x(6);
a=theta+b.phase;
T=sqrt(2/3)*[cos(a).';-sin(a).';ones(1,3)/sqrt(2)];
id=(p(5)*x(1)-p(6)*x(4))/(p(2)*p(5)-p(6)^2);
ifield=(p(2)*x(4)-p(6)*x(1))/(p(2)*p(5)-p(6)^2);
iq=x(2)/p(3); i0=x(3)/p(4);
v=x(11:13); vdq=T*v;
igen=-T.'*[id;iq;i0];
pGap=-omega*(x(1)*iq-x(2)*id);
dx=zeros(16,1,'like',x);
dx(1)=vdq(1)-p(1)*id+omega*x(2);
dx(2)=vdq(2)-p(1)*iq-omega*x(1);
dx(3)=vdq(3)-p(1)*i0;
dx(4)=p(12)*x(9)-p(7)*ifield;
dx(5)=omega;
dx(6)=(x(10)-pGap/p(9)-p(11)*(x(6)-1))/(2*p(10)*x(6));
if real(x(7))>1e-12, vm=sqrt(x(7)); else, vm=1e-6; end
error=vref-vm;
raw=x(8)+b.avr(3)*error;
command=limit(raw,b.avr(6),b.avr(7));
dx(7)=((v.'*v)/(3*b.Vphase^2)-x(7))/b.avr(2);
dx(8)=b.avr(4)*error+b.avr(5)*(command-raw);
dx(9)=(command-x(9))/b.avr(8);
pm=limit(b.pm0+(1-x(6))/b.governor(1),b.governor(3),b.governor(4));
dx(10)=(pm-x(10))/b.governor(2);
grid=sqrt(2)*b.Vphase*cos(b.wb*t+b.phase);
dx(11:13)=(igen-x(14:16)-g*v)/b.network(1);
dx(14:16)=(v-grid-b.network(2)*x(14:16))/b.network(3);
end
function y=limit(x,lower,upper)
if real(x)<lower, y=lower;
elseif real(x)>upper, y=upper;
else, y=x;
end
end
